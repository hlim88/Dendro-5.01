/**
 * @file bssnCtx.cpp
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief BSSN contex file. 
 * @version 0.1
 * @date 2019-12-20
 * 
 * School of Computing, University of Utah. 
 * @copyright Copyright (c) 2019
 * 
 */

#include "bssnCtx.h"

namespace bssn
{
    BSSNCtx::BSSNCtx(ot::Mesh* pMesh) : Ctx()
    {   
        m_uiMesh = pMesh;
        // variable allocation for evolution variables
        m_uiEVar = this->create_vec(ts::CTXVType::EVOLUTION,true,false,false,BSSN_NUM_VARS);

        // variable allocation for constraint variables. 
        m_uiCVar = this->create_vec(ts::CTXVType::CONSTRAINT,true,false,false,BSSN_CONSTRAINT_NUM_VARS);

        // evars unzip in (0) evars unzip out (1) 
        m_uiEUnzip[0] = this->create_vec(ts::CTXVType::EVOLUTION,false,true,false,BSSN_NUM_VARS);
        m_uiEUnzip[1] = this->create_vec(ts::CTXVType::EVOLUTION,false,true,false,BSSN_NUM_VARS);

        // // constraint var unzip out (0)
        m_uiCUnzip[0] = this->create_vec(ts::CTXVType::CONSTRAINT,false,true,false,BSSN_CONSTRAINT_NUM_VARS); 

        m_uiTinfo._m_uiStep=0;
        m_uiTinfo._m_uiT = 0;
        m_uiTinfo._m_uiTb = bssn::BSSN_RK_TIME_BEGIN;
        m_uiTinfo._m_uiTe = bssn::BSSN_RK_TIME_END;
        m_uiTinfo._m_uiTh = bssn::BSSN_RK45_TIME_STEP_SIZE;     

        m_uiElementOrder = bssn::BSSN_ELE_ORDER;

        m_uiMinPt = Point(bssn::BSSN_GRID_MIN_X,bssn::BSSN_GRID_MIN_Y,bssn::BSSN_GRID_MIN_Z);
        m_uiMaxPt = Point(bssn::BSSN_GRID_MAX_X,bssn::BSSN_GRID_MAX_Y,bssn::BSSN_GRID_MAX_Z);
        
        return;

    }

    BSSNCtx::~BSSNCtx()
    {
        this->destroy_vec(m_uiEVar);
        this->destroy_vec(m_uiCVar);

        this->destroy_vec(m_uiEUnzip[0]);
        this->destroy_vec(m_uiEUnzip[1]);
        this->destroy_vec(m_uiCUnzip[0]);

    }

    int BSSNCtx::initialize()
    {
        if(bssn::BSSN_RESTORE_SOLVER)
        {
            this->restore_checkpt();
            return 0;
        }

        this->init_grid();

        bool isRefine=false;
        DendroIntL oldElements,oldElements_g;
        DendroIntL newElements,newElements_g;

        DendroIntL oldGridPoints,oldGridPoints_g;
        DendroIntL newGridPoints,newGridPoints_g;

        unsigned int iterCount=1;
        const unsigned int max_iter=bssn::BSSN_INIT_GRID_ITER;
        const unsigned int rank_global=m_uiMesh->getMPIRankGlobal();
        MPI_Comm gcomm = m_uiMesh->getMPIGlobalCommunicator();
        
        do
        {

            isRefine=this->is_remesh();
            if(isRefine)
            {
                ot::Mesh* newMesh = this->remesh(bssn::BSSN_DENDRO_GRAIN_SZ, bssn::BSSN_LOAD_IMB_TOL,bssn::BSSN_SPLIT_FIX);
                
                oldElements = m_uiMesh->getNumLocalMeshElements();
                newElements = newMesh->getNumLocalMeshElements();

                oldGridPoints =m_uiMesh->getNumLocalMeshNodes();
                newGridPoints =newMesh->getNumLocalMeshNodes();

                par::Mpi_Allreduce(&oldElements,&oldElements_g,1,MPI_SUM,gcomm);
                par::Mpi_Allreduce(&newElements,&newElements_g,1,MPI_SUM,gcomm);

                par::Mpi_Allreduce(&oldGridPoints,&oldGridPoints_g,1,MPI_SUM,m_uiMesh->getMPIGlobalCommunicator());
                par::Mpi_Allreduce(&newGridPoints,&newGridPoints_g,1,MPI_SUM,m_uiMesh->getMPIGlobalCommunicator());
                
                if(!rank_global)
                {
                    std::cout<<"[bssnCtx][initial grid convg. ] iter : "<<iterCount<<" (Remesh triggered) ->  old mesh : "<<oldElements_g<<" new mesh : "<<newElements_g<<std::endl;

                    std::cout<<"[bssnCtx][initial grid convg. ] iter : "<<iterCount<<" (Remesh triggered) ->  old mesh (zip nodes) : "<<oldGridPoints_g<<" new mesh (zip nodes) : "<<newGridPoints_g<<std::endl;
                }
                  

                this->grid_transfer(newMesh,true,false,false);
                this->update_app_vars();

                std::swap(m_uiMesh,newMesh);
                delete newMesh;
            }
            
            iterCount+=1;

        } while(isRefine && (newElements_g!=oldElements_g || newGridPoints_g!=oldGridPoints_g) && (iterCount<max_iter) );

        this->init_grid();

        unsigned int lmin, lmax;
        m_uiMesh->computeMinMaxLevel(lmin,lmax);
        bssn::BSSN_RK45_TIME_STEP_SIZE=bssn::BSSN_CFL_FACTOR*((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))));
        m_uiTinfo._m_uiTh=bssn::BSSN_RK45_TIME_STEP_SIZE;

        if(!m_uiMesh->getMPIRankGlobal())
        {
            std::cout<<"================= Grid Info (After init grid converge):======================================================="<<std::endl;
            std::cout<<"lmin: "<<lmin<<" lmax:"<<lmax<<std::endl;
            std::cout<<"dx: "<<((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))))<<std::endl;
            std::cout<<"dt: "<<bssn::BSSN_CFL_FACTOR*((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))))<<std::endl;
            std::cout<<"==============================================================================================================="<<std::endl;
        }
        
        return 0;
                    
    }

    int BSSNCtx::init_grid()
    {  

        unsigned int nodeLookUp_CG;
        unsigned int nodeLookUp_DG;
        DendroScalar x,y,z,len;
        const ot::TreeNode * pNodes=&(*(m_uiMesh->getAllElements().begin()));
        unsigned int ownerID,ii_x,jj_y,kk_z;
        unsigned int eleOrder=m_uiMesh->getElementOrder();
        const unsigned int * e2n_cg=&(*(m_uiMesh->getE2NMapping().begin()));
        const unsigned int * e2n_dg=&(*(m_uiMesh->getE2NMapping_DG().begin()));
        const unsigned int nPe=m_uiMesh->getNumNodesPerElement();
        const unsigned int nodeLocalBegin=m_uiMesh->getNodeLocalBegin();
        const unsigned int nodeLocalEnd=m_uiMesh->getNodeLocalEnd();


        DendroScalar* var=new double[bssn::BSSN_NUM_VARS];
        DendroScalar** zipIn = new DendroScalar*[bssn::BSSN_NUM_VARS];
        m_uiEVar.Get2DArray(zipIn,true);

        DendroScalar mp, mm, mp_adm, mm_adm, E, J1, J2, J3;
        // set the TP communicator. 
        if(bssn::BSSN_ID_TYPE==0)
        {
            TP_MPI_COMM=m_uiMesh->getMPIGlobalCommunicator();
            TwoPunctures((double)0,(double)0,(double)0,var,&mp, &mm, &mp_adm, &mm_adm, &E, &J1, &J2, &J3);
        }
        
        for(unsigned int elem=m_uiMesh->getElementLocalBegin(); elem<m_uiMesh->getElementLocalEnd(); elem++)
        {
            for(unsigned int k=0; k<(eleOrder+1); k++)
                for(unsigned int j=0; j<(eleOrder+1); j++ )
                    for(unsigned int i=0; i<(eleOrder+1); i++)
                    {
                        nodeLookUp_CG=e2n_cg[elem*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                        if(nodeLookUp_CG>=nodeLocalBegin && nodeLookUp_CG<nodeLocalEnd)
                        {
                            nodeLookUp_DG=e2n_dg[elem*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                            m_uiMesh->dg2eijk(nodeLookUp_DG,ownerID,ii_x,jj_y,kk_z);
                            len=(double)(1u<<(m_uiMaxDepth-pNodes[ownerID].getLevel()));
                            x=pNodes[ownerID].getX()+ ii_x*(len/(eleOrder));
                            y=pNodes[ownerID].getY()+ jj_y*(len/(eleOrder));
                            z=pNodes[ownerID].getZ()+ kk_z*(len/(eleOrder));
                            
                            if (bssn::BSSN_ID_TYPE == 0) {
                                x = GRIDX_TO_X(x); y = GRIDY_TO_Y(y); z = GRIDZ_TO_Z(z);
                                TwoPunctures((double)x,(double)y,(double)z,var,
                                            &mp, &mm, &mp_adm, &mm_adm, &E, &J1, &J2, &J3);
                            }
                            else if (bssn::BSSN_ID_TYPE == 1) {
                                bssn::punctureData((double)x,(double)y,(double)z,var);
                            }
                            else if (bssn::BSSN_ID_TYPE == 2) {
                                bssn::KerrSchildData((double)x,(double)y,(double)z,var);
                            }
                            else if (bssn::BSSN_ID_TYPE == 3) {
                                bssn::noiseData((double)x,(double)y,(double)z,var);
                            }
                            else if (bssn::BSSN_ID_TYPE == 4) {
                                bssn::fake_initial_data((double)x,(double)y,(double)z,var);
                            }
                            else {
                                std::cout<<"Unknown ID type"<<std::endl;
                            }
                            for(unsigned int v=0; v<bssn::BSSN_NUM_VARS; v++)
                                zipIn[v][nodeLookUp_CG]=var[v];
                        }

                    }
        }
        
        for(unsigned int node=m_uiMesh->getNodeLocalBegin(); node<m_uiMesh->getNodeLocalEnd(); node++)
            enforce_bssn_constraints(zipIn,node);


        delete [] var;
        delete [] zipIn;

        #ifdef BSSN_EXTRACT_BH_LOCATIONS
            m_uiBHLoc[0]=Point(bssn::BH1.getBHCoordX(),bssn::BH1.getBHCoordY(),bssn::BH1.getBHCoordZ());
            m_uiBHLoc[1]=Point(bssn::BH2.getBHCoordX(),bssn::BH2.getBHCoordY(),bssn::BH2.getBHCoordZ());
        #endif

        return 0;

    }
    
    int BSSNCtx::finalize()
    {
        return 0;    
    }

    int BSSNCtx::rhs(DVec* in , DVec* out, unsigned int sz , DendroScalar time)
    {
        // all the variables should be packed together. 
        assert(sz==1);
        DendroScalar ** sVar;
        in[0].Get2DArray(sVar,false);

        for(unsigned int node=m_uiMesh->getNodeLocalBegin(); node< m_uiMesh->getNodeLocalEnd(); node++)
            enforce_bssn_constraints(sVar, node);

        delete [] sVar;

        this->unzip(in[0],m_uiEUnzip[0] , bssn::BSSN_ASYNC_COMM_K);
        
        DendroScalar ** unzipIn;
        DendroScalar **  unzipOut; 
        
        m_uiEUnzip[0].Get2DArray(unzipIn , false);
        m_uiEUnzip[1].Get2DArray(unzipOut ,false);

        const ot::Block* blkList = m_uiMesh->getLocalBlockList().data();
        const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();
        
        bssnRHS(unzipOut,(const DendroScalar**)unzipIn,blkList,numBlocks);

        this->zip(m_uiEUnzip[1], out[0], bssn::BSSN_ASYNC_COMM_K);

        delete [] unzipIn;
        delete [] unzipOut;
        

        return 0;

    }

    int BSSNCtx::rhs_blkwise(DVec in, DVec out, const unsigned int* const blkIDs, unsigned int numIds, DendroScalar*  blk_time) const
    {

        DendroScalar ** unzipIn;
        DendroScalar **  unzipOut; 

        assert(in.GetDof() == out.GetDof());
        assert(in.IsUnzip() == out.IsUnzip());
        
        in.Get2DArray(unzipIn , false);
        out.Get2DArray(unzipOut ,false);

        const ot::Block* blkList = m_uiMesh->getLocalBlockList().data();
        const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();

        unsigned int offset;
        double ptmin[3], ptmax[3];
        unsigned int sz[3];
        unsigned int bflag;
        double dx,dy,dz;
        const Point pt_min(bssn::BSSN_COMPD_MIN[0],bssn::BSSN_COMPD_MIN[1],bssn::BSSN_COMPD_MIN[2]);
        const Point pt_max(bssn::BSSN_COMPD_MAX[0],bssn::BSSN_COMPD_MAX[1],bssn::BSSN_COMPD_MAX[2]);

        for(unsigned int i=0; i < numIds; i++)
        {
            const unsigned int blk = blkIDs[i];
            assert(blk<numBlocks);

            offset=blkList[blk].getOffset();
            sz[0]=blkList[blk].getAllocationSzX();
            sz[1]=blkList[blk].getAllocationSzY();
            sz[2]=blkList[blk].getAllocationSzZ();

            bflag=blkList[blk].getBlkNodeFlag();

            dx=blkList[blk].computeDx(pt_min,pt_max);
            dy=blkList[blk].computeDy(pt_min,pt_max);
            dz=blkList[blk].computeDz(pt_min,pt_max);

            ptmin[0]=GRIDX_TO_X(blkList[blk].getBlockNode().minX())-3*dx;
            ptmin[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().minY())-3*dy;
            ptmin[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ())-3*dz;

            ptmax[0]=GRIDX_TO_X(blkList[blk].getBlockNode().maxX())+3*dx;
            ptmax[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().maxY())+3*dy;
            ptmax[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ())+3*dz;

            #ifdef BSSN_RHS_STAGED_COMP
                bssnrhs_sep(unzipOut, (const double **)unzipIn, offset, ptmin, ptmax, sz, bflag);
            #else
                bssnrhs(unzipOut, (const double **)unzipIn, offset, ptmin, ptmax, sz, bflag);
            #endif

        }

        delete [] unzipIn;
        delete [] unzipOut;

        return 0;

        
    }

    int BSSNCtx::rhs_blk(const DendroScalar* in, DendroScalar* out, unsigned int dof, unsigned int local_blk_id, DendroScalar  blk_time) const 
    {
        //return 0;
        //std::cout<<"bssn_rhs"<<std::endl;
        DendroScalar **  unzipIn   = new DendroScalar*[dof];
        DendroScalar **  unzipOut  = new DendroScalar*[dof]; 

        const unsigned int blk = local_blk_id;

        const ot::Block* blkList = m_uiMesh->getLocalBlockList().data();
        const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();
        double ptmin[3], ptmax[3];
        unsigned int sz[3];
        unsigned int bflag;
        double dx,dy,dz;
        const Point pt_min(bssn::BSSN_COMPD_MIN[0],bssn::BSSN_COMPD_MIN[1],bssn::BSSN_COMPD_MIN[2]);
        const Point pt_max(bssn::BSSN_COMPD_MAX[0],bssn::BSSN_COMPD_MAX[1],bssn::BSSN_COMPD_MAX[2]);

        sz[0]=blkList[blk].getAllocationSzX();
        sz[1]=blkList[blk].getAllocationSzY();
        sz[2]=blkList[blk].getAllocationSzZ();

        const unsigned int NN  = sz[0]*sz[1]*sz[2];

        for(unsigned int v =0; v < dof; v++)
        {
            unzipIn[v] = (DendroScalar*) (in + v*NN);
            unzipOut[v] = (DendroScalar*) (out + v*NN);
        }

        
        bflag=blkList[blk].getBlkNodeFlag();
        const unsigned int pw = blkList[blk].get1DPadWidth();

        // if(!bflag)
        // {
        //     // for(unsigned int node=0; node < NN; node++)
        //     //     enforce_bssn_constraints(unzipIn, node);
        //     for(unsigned int k=pw; k < sz[2]-pw; k++)
        //     for(unsigned int j=pw; j < sz[1]-pw; j++)
        //     for(unsigned int i=pw; i < sz[0]-pw; i++)
        //     {
        //         const unsigned nid = k*sz[1]*sz[0] + j*sz[0] + i;
        //         enforce_bssn_constraints(unzipIn,nid);
        //     }
            
        // }else
        // {
        //     // note that we can apply enforce bssn constraints in the right padd, at the left boundary block,
        //     // currently we only apply internal parts of the boundary blocks. 
        //     for(unsigned int k=pw; k < sz[2]-pw; k++)
        //     for(unsigned int j=pw; j < sz[1]-pw; j++)
        //     for(unsigned int i=pw; i < sz[0]-pw; i++)
        //     {
        //         const unsigned nid = k*sz[1]*sz[0] + j*sz[0] + i;
        //         enforce_bssn_constraints(unzipIn,nid);
        //     }
                
        // }

        dx=blkList[blk].computeDx(pt_min,pt_max);
        dy=blkList[blk].computeDy(pt_min,pt_max);
        dz=blkList[blk].computeDz(pt_min,pt_max);

        ptmin[0]=GRIDX_TO_X(blkList[blk].getBlockNode().minX())-3*dx;
        ptmin[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().minY())-3*dy;
        ptmin[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ())-3*dz;

        ptmax[0]=GRIDX_TO_X(blkList[blk].getBlockNode().maxX())+3*dx;
        ptmax[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().maxY())+3*dy;
        ptmax[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ())+3*dz;

        #ifdef BSSN_RHS_STAGED_COMP
            bssnrhs_sep(unzipOut, (const DendroScalar **)unzipIn, 0, ptmin, ptmax, sz, bflag);
        #else
            bssnrhs(unzipOut, (const DendroScalar **)unzipIn, 0, ptmin, ptmax, sz, bflag);
        #endif

        delete [] unzipIn;
        delete [] unzipOut;

        return 0;

    }

    int BSSNCtx::pre_stage_blk(DendroScalar* in, unsigned int dof, unsigned int local_blk_id, DendroScalar  blk_time) const
    {

        DendroScalar **  unzipIn   = new DendroScalar*[dof];
        const unsigned int blk = local_blk_id;

        const ot::Block* blkList = m_uiMesh->getLocalBlockList().data();
        const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();
        double ptmin[3], ptmax[3];
        unsigned int sz[3];
        unsigned int bflag;
        double dx,dy,dz;
        const Point pt_min(bssn::BSSN_COMPD_MIN[0],bssn::BSSN_COMPD_MIN[1],bssn::BSSN_COMPD_MIN[2]);
        const Point pt_max(bssn::BSSN_COMPD_MAX[0],bssn::BSSN_COMPD_MAX[1],bssn::BSSN_COMPD_MAX[2]);

        sz[0]=blkList[blk].getAllocationSzX();
        sz[1]=blkList[blk].getAllocationSzY();
        sz[2]=blkList[blk].getAllocationSzZ();

        const unsigned int NN  = sz[0]*sz[1]*sz[2];

        for(unsigned int v =0; v < dof; v++)
        {
            unzipIn[v] = (DendroScalar*) (in + v*NN);
        }

        
        bflag=blkList[blk].getBlkNodeFlag();
        const unsigned int pw = blkList[blk].get1DPadWidth();
        
        if(!bflag)
        {
            for(unsigned int node=0; node < NN; node++)
                enforce_bssn_constraints(unzipIn, node);
            // for(unsigned int k=pw; k < sz[2]-pw; k++)
            // for(unsigned int j=pw; j < sz[1]-pw; j++)
            // for(unsigned int i=pw; i < sz[0]-pw; i++)
            // {
            //     const unsigned nid = k*sz[1]*sz[0] + j*sz[0] + i;
            //     enforce_bssn_constraints(unzipIn,nid);
            // }
            
        }else
        {
            // note that we can apply enforce bssn constraints in the right padd, at the left boundary block,
            // currently we only apply internal parts of the boundary blocks. 
            for(unsigned int k=pw; k < sz[2]-pw; k++)
            for(unsigned int j=pw; j < sz[1]-pw; j++)
            for(unsigned int i=pw; i < sz[0]-pw; i++)
            {
                const unsigned nid = k*sz[1]*sz[0] + j*sz[0] + i;
                enforce_bssn_constraints(unzipIn,nid);
            }
                
        }

        delete [] unzipIn;
        return 0;
    }
    
    int BSSNCtx::post_stage_blk(DendroScalar* in, unsigned int dof, unsigned int local_blk_id, DendroScalar  blk_time) const
    {
        return 0;
    }

    int BSSNCtx::pre_timestep_blk(DendroScalar* in, unsigned int dof, unsigned int local_blk_id, DendroScalar  blk_time) const
    {
        return 0;
    }
    
    int BSSNCtx::post_timestep_blk(DendroScalar* in, unsigned int dof, unsigned int local_blk_id, DendroScalar  blk_time) const
    {
        DendroScalar **  unzipIn   = new DendroScalar*[dof];
        const unsigned int blk = local_blk_id;

        const ot::Block* blkList = m_uiMesh->getLocalBlockList().data();
        const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();
        double ptmin[3], ptmax[3];
        unsigned int sz[3];
        unsigned int bflag;
        double dx,dy,dz;
        const Point pt_min(bssn::BSSN_COMPD_MIN[0],bssn::BSSN_COMPD_MIN[1],bssn::BSSN_COMPD_MIN[2]);
        const Point pt_max(bssn::BSSN_COMPD_MAX[0],bssn::BSSN_COMPD_MAX[1],bssn::BSSN_COMPD_MAX[2]);

        sz[0]=blkList[blk].getAllocationSzX();
        sz[1]=blkList[blk].getAllocationSzY();
        sz[2]=blkList[blk].getAllocationSzZ();

        const unsigned int NN  = sz[0]*sz[1]*sz[2];

        for(unsigned int v =0; v < dof; v++)
        {
            unzipIn[v] = (DendroScalar*) (in + v*NN);
        }

        
        bflag=blkList[blk].getBlkNodeFlag();
        const unsigned int pw = blkList[blk].get1DPadWidth();
        
        if(!bflag)
        {
            for(unsigned int node=0; node < NN; node++)
                enforce_bssn_constraints(unzipIn, node);
            // for(unsigned int k=pw; k < sz[2]-pw; k++)
            // for(unsigned int j=pw; j < sz[1]-pw; j++)
            // for(unsigned int i=pw; i < sz[0]-pw; i++)
            // {
            //     const unsigned nid = k*sz[1]*sz[0] + j*sz[0] + i;
            //     enforce_bssn_constraints(unzipIn,nid);
            // }
            
        }else
        {
            // note that we can apply enforce bssn constraints in the right padd, at the left boundary block,
            // currently we only apply internal parts of the boundary blocks. 
            for(unsigned int k=pw; k < sz[2]-pw; k++)
            for(unsigned int j=pw; j < sz[1]-pw; j++)
            for(unsigned int i=pw; i < sz[0]-pw; i++)
            {
                const unsigned nid = k*sz[1]*sz[0] + j*sz[0] + i;
                enforce_bssn_constraints(unzipIn,nid);
            }
                
        }

        delete [] unzipIn;
        return 0;
    }

    int BSSNCtx::write_vtu()
    {

        unzip(m_uiEVar,m_uiEUnzip[0],BSSN_ASYNC_COMM_K);
        
        DendroScalar ** evolUnzipVar = NULL;
        DendroScalar ** consUnzipVar = NULL;
        DendroScalar ** consVar = NULL;
        DendroScalar ** evolVar = NULL;

        m_uiEUnzip[0].Get2DArray(evolUnzipVar,false);
        m_uiCUnzip[0].Get2DArray(consUnzipVar,false);
       
        m_uiEVar.Get2DArray(evolVar,false);
        m_uiCVar.Get2DArray(consVar,false);

        #ifdef BSSN_COMPUTE_CONSTRAINTS

            const std::vector<ot::Block> blkList=m_uiMesh->getLocalBlockList();

            unsigned int offset;
            double ptmin[3], ptmax[3];
            unsigned int sz[3];
            unsigned int bflag;
            double dx,dy,dz;
            const Point pt_min(bssn::BSSN_COMPD_MIN[0],bssn::BSSN_COMPD_MIN[1],bssn::BSSN_COMPD_MIN[2]);
            const Point pt_max(bssn::BSSN_COMPD_MAX[0],bssn::BSSN_COMPD_MAX[1],bssn::BSSN_COMPD_MAX[2]);

            for(unsigned int blk=0; blk<blkList.size(); blk++)
            {
                offset=blkList[blk].getOffset();
                sz[0]=blkList[blk].getAllocationSzX();
                sz[1]=blkList[blk].getAllocationSzY();
                sz[2]=blkList[blk].getAllocationSzZ();

                bflag=blkList[blk].getBlkNodeFlag();

                dx=blkList[blk].computeDx(pt_min,pt_max);
                dy=blkList[blk].computeDy(pt_min,pt_max);
                dz=blkList[blk].computeDz(pt_min,pt_max);

                ptmin[0]=GRIDX_TO_X(blkList[blk].getBlockNode().minX())-3*dx;
                ptmin[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().minY())-3*dy;
                ptmin[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ())-3*dz;

                ptmax[0]=GRIDX_TO_X(blkList[blk].getBlockNode().maxX())+3*dx;
                ptmax[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().maxY())+3*dy;
                ptmax[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ())+3*dz;

                physical_constraints(consUnzipVar, (const DendroScalar **) evolUnzipVar, offset, ptmin, ptmax, sz, bflag);
            }

            /*double consVecMin[bssn::BSSN_CONSTRAINT_NUM_VARS];
            double consVecMax[bssn::BSSN_CONSTRAINT_NUM_VARS];*/
            double constraintMaskedL2[bssn::BSSN_CONSTRAINT_NUM_VARS];
            
            this->zip(m_uiCUnzip[0],m_uiCVar,bssn::BSSN_ASYNC_COMM_K);
            


            bssn::extractConstraints(m_uiMesh,(const DendroScalar **)consVar,evolVar[BHLOC::EXTRACTION_VAR_ID],BHLOC::EXTRACTION_TOL,m_uiTinfo._m_uiStep,m_uiTinfo._m_uiT);
            #ifndef BSSN_KERR_SCHILD_TEST
                #ifdef BSSN_EXTRACT_GRAVITATIONAL_WAVES
                GW::extractFarFieldPsi4(m_uiMesh,(const DendroScalar **)consVar,m_uiTinfo._m_uiStep,m_uiTinfo._m_uiT);
                #endif
            #endif

        #endif

        #ifdef BSSN_ENABLE_VTU_OUTPUT

            if((m_uiTinfo._m_uiStep % bssn::BSSN_IO_OUTPUT_FREQ)==0)
            {
                std::vector<std::string> pDataNames;
                const unsigned int numConstVars = bssn::BSSN_NUM_CONST_VARS_VTU_OUTPUT;
                const unsigned int numEvolVars  = bssn::BSSN_NUM_EVOL_VARS_VTU_OUTPUT;

                double *pData[(numConstVars+numEvolVars)];

                for(unsigned int i=0; i<numEvolVars; i++)
                {
                    pDataNames.push_back(std::string(bssn::BSSN_VAR_NAMES[BSSN_VTU_OUTPUT_EVOL_INDICES[i]]));
                    pData[i]=evolVar[BSSN_VTU_OUTPUT_EVOL_INDICES[i]];
                }


                for(unsigned int i=0; i<numConstVars; i++)
                {
                    pDataNames.push_back(std::string(bssn::BSSN_CONSTRAINT_VAR_NAMES[BSSN_VTU_OUTPUT_CONST_INDICES[i]]));
                    pData[numEvolVars+i]=consVar[BSSN_VTU_OUTPUT_CONST_INDICES[i]];
                }

                std::vector<char*> pDataNames_char;
                pDataNames_char.reserve(pDataNames.size());

                for(unsigned int  i = 0; i < pDataNames.size(); i++)
                    pDataNames_char.push_back(const_cast<char*>(pDataNames[i].c_str()));

                const char * fDataNames[]= {"Time","Cycle"};
                const double fData[]= {m_uiTinfo._m_uiT,(double)m_uiTinfo._m_uiStep};

                char fPrefix[256];
                sprintf(fPrefix,"%s_%d",bssn::BSSN_VTU_FILE_PREFIX.c_str(),m_uiTinfo._m_uiStep);

                if(bssn::BSSN_VTU_Z_SLICE_ONLY)
                {
                    unsigned int s_val[3]= {1u<<(m_uiMaxDepth-1), 1u<<(m_uiMaxDepth-1), 1u<<(m_uiMaxDepth-1)};
                    unsigned int s_norm[3] ={0,0,1};
                    io::vtk::mesh2vtu_slice(m_uiMesh,s_val, s_norm, fPrefix,2,fDataNames,fData,(numEvolVars+numConstVars),(const char **)&pDataNames_char[0],(const double **)pData);
                }
                else
                    io::vtk::mesh2vtuFine(m_uiMesh,fPrefix,2,fDataNames,fData,(numEvolVars+numConstVars),(const char **)&pDataNames_char[0],(const double **)pData);
                

            }

            
        #endif


        #ifdef BSSN_EXTRACT_BH_LOCATIONS
            bssn::writeBHCoordinates((const ot::Mesh *)m_uiMesh,(const Point *) m_uiBHLoc,2,m_uiTinfo._m_uiStep,m_uiTinfo._m_uiT);
        #endif

        delete [] evolUnzipVar;
        delete [] consUnzipVar;
        delete [] evolVar;
        delete [] consVar;

        return 0;
    }

    int BSSNCtx::write_checkpt()
    {   
        if(m_uiMesh->isActive())
        {
            unsigned int cpIndex;
            (m_uiTinfo._m_uiStep %(2*bssn::BSSN_CHECKPT_FREQ)==0) ? cpIndex=0 : cpIndex=1; // to support alternate file writing.

            const bool is_merged = ((bssn::BSSN_BH_LOC[0]-bssn::BSSN_BH_LOC[1]).abs() < 0.1);
            if(is_merged && !bssn::BSSN_MERGED_CHKPT_WRITTEN)
            {
                cpIndex=3;
                bssn::BSSN_MERGED_CHKPT_WRITTEN=true;
            }
            
            unsigned int rank=m_uiMesh->getMPIRank();
            unsigned int npes=m_uiMesh->getMPICommSize();

            DendroScalar** eVar = NULL;
            m_uiEVar.Get2DArray(eVar,false);
            

            char fName[256];
            const ot::TreeNode * pNodes=&(*(m_uiMesh->getAllElements().begin()+m_uiMesh->getElementLocalBegin()));
            sprintf(fName,"%s_octree_%d_%d.oct",bssn::BSSN_CHKPT_FILE_PREFIX.c_str(),cpIndex,rank);
            io::checkpoint::writeOctToFile(fName,pNodes,m_uiMesh->getNumLocalMeshElements());

            unsigned int numVars=bssn::BSSN_NUM_VARS;
            const char ** varNames=bssn::BSSN_VAR_NAMES;

            /*for(unsigned int i=0;i<numVars;i++)
            {
                sprintf(fName,"%s_%s_%d_%d.var",fNamePrefix,varNames[i],cpIndex,rank);
                io::checkpoint::writeVecToFile(fName,m_uiMesh,m_uiPrevVar[i]);
            }*/

            sprintf(fName,"%s_%d_%d.var",bssn::BSSN_CHKPT_FILE_PREFIX.c_str(),cpIndex,rank);
            io::checkpoint::writeVecToFile(fName,m_uiMesh,(const double **)eVar,bssn::BSSN_NUM_VARS);

            if(!rank)
            {
                sprintf(fName,"%s_step_%d.cp",bssn::BSSN_CHKPT_FILE_PREFIX.c_str(),cpIndex);
                std::cout<<"[BSSNCtx] \t writing checkpoint file : "<<fName<<std::endl;
                std::ofstream outfile(fName);
                if(!outfile) {
                    std::cout<<fName<<" file open failed "<<std::endl;
                    return 0;
                }

                json checkPoint;
                checkPoint["DENDRO_TS_TIME_BEGIN"]    = m_uiTinfo._m_uiTb;
                checkPoint["DENDRO_TS_TIME_END"]      = m_uiTinfo._m_uiTe;
                checkPoint["DENDRO_TS_ELEMENT_ORDER"] = m_uiElementOrder;

                checkPoint["DENDRO_TS_TIME_CURRENT"]    = m_uiTinfo._m_uiT;
                checkPoint["DENDRO_TS_STEP_CURRENT"]    = m_uiTinfo._m_uiStep;
                checkPoint["DENDRO_TS_TIME_STEP_SIZE"]  = m_uiTinfo._m_uiTh;
                checkPoint["DENDRO_TS_LAST_IO_TIME"]    = m_uiTinfo._m_uiT;

                checkPoint["DENDRO_TS_WAVELET_TOLERANCE"]=bssn::BSSN_WAVELET_TOL;
                checkPoint["DENDRO_TS_LOAD_IMB_TOLERANCE"]=bssn::BSSN_LOAD_IMB_TOL;
                checkPoint["DENDRO_TS_NUM_VARS"]=numVars; // number of variables to restore.
                checkPoint["DENDRO_TS_ACTIVE_COMM_SZ"]=m_uiMesh->getMPICommSize(); // (note that rank 0 is always active).
                
                checkPoint["DENDRO_BH1_X"]=m_uiBHLoc[0].x();
                checkPoint["DENDRO_BH1_Y"]=m_uiBHLoc[0].y();
                checkPoint["DENDRO_BH1_Z"]=m_uiBHLoc[0].z();
                
                
                checkPoint["DENDRO_BH2_X"]=m_uiBHLoc[1].x();
                checkPoint["DENDRO_BH2_Y"]=m_uiBHLoc[1].y();
                checkPoint["DENDRO_BH2_Z"]=m_uiBHLoc[1].z();
                

                outfile<<std::setw(4)<<checkPoint<<std::endl;
                outfile.close();

            }
            
            delete [] eVar;
        }
        return 0;
    }

    int BSSNCtx::restore_checkpt()
    {
        unsigned int numVars=0;
        std::vector<ot::TreeNode> octree;
        json checkPoint;

        int rank;
        int npes;
        MPI_Comm comm=m_uiMesh->getMPIGlobalCommunicator();
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        unsigned int activeCommSz;

        char fName[256];
        unsigned int restoreStatus=0;
        unsigned int restoreStatusGlobal=0; // 0 indicates successfully restorable.

        ot::Mesh* newMesh;
        unsigned int restoreStep[2];
        restoreStep[0]=0;
        restoreStep[1]=0;

        unsigned int restoreFileIndex=0;

        for(unsigned int cpIndex=0; cpIndex<2; cpIndex++) {

            restoreStatus=0;

            if(!rank)
            {
                sprintf(fName,"%s_step_%d.cp", bssn::BSSN_CHKPT_FILE_PREFIX.c_str() , cpIndex);
                std::ifstream infile(fName);
                if(!infile) {
                    std::cout<<fName<<" file open failed "<<std::endl;
                    restoreStatus=1;
                }


                if(restoreStatus==0)
                {
                    infile>>checkPoint;
                    m_uiTinfo._m_uiTb    = checkPoint["DENDRO_TS_TIME_BEGIN"];
                    m_uiTinfo._m_uiTe    = checkPoint["DENDRO_TS_TIME_END"];
                    m_uiTinfo._m_uiT     = checkPoint["DENDRO_TS_TIME_CURRENT"];
                    m_uiTinfo._m_uiStep  = checkPoint["DENDRO_TS_STEP_CURRENT"]; 
                    m_uiTinfo._m_uiTh    = checkPoint["DENDRO_TS_TIME_STEP_SIZE"];
                    m_uiElementOrder     = checkPoint["DENDRO_TS_ELEMENT_ORDER"];
                    
                    bssn::BSSN_WAVELET_TOL=checkPoint["DENDRO_TS_WAVELET_TOLERANCE"];
                    bssn::BSSN_LOAD_IMB_TOL=checkPoint["DENDRO_TS_LOAD_IMB_TOLERANCE"];
                    
                    numVars=checkPoint["DENDRO_TS_NUM_VARS"];
                    activeCommSz=checkPoint["DENDRO_TS_ACTIVE_COMM_SZ"];
                    
                    m_uiBHLoc[0]=Point((double)checkPoint["DENDRO_BH1_X"],(double)checkPoint["DENDRO_BH1_Y"],(double)checkPoint["DENDRO_BH1_Z"]);
                    m_uiBHLoc[1]=Point((double)checkPoint["DENDRO_BH2_X"],(double)checkPoint["DENDRO_BH2_Y"],(double)checkPoint["DENDRO_BH2_Z"]);
                    restoreStep[cpIndex]=m_uiTinfo._m_uiStep;

                    }
              }

            }

            if(!rank)
            {
                if(restoreStep[0]<restoreStep[1])
                    restoreFileIndex=1;
                else
                    restoreFileIndex=0;
            }

            par::Mpi_Bcast(&restoreFileIndex,1,0,comm);

            restoreStatus=0;
            octree.clear();
            if(!rank) std::cout<<"[BSSNCtx] :  Trying to restore from checkpoint index : "<<restoreFileIndex<<std::endl;
        
            if(!rank)
            {
                sprintf(fName,"%s_step_%d.cp", bssn::BSSN_CHKPT_FILE_PREFIX.c_str(), restoreFileIndex);
                std::ifstream infile(fName);
                if(!infile) {
                    std::cout<<fName<<" file open failed "<<std::endl;
                    restoreStatus=1;
                }


                if(restoreStatus==0)
                {
                    infile>>checkPoint;
                    m_uiTinfo._m_uiTb    = checkPoint["DENDRO_TS_TIME_BEGIN"];
                    m_uiTinfo._m_uiTe    = checkPoint["DENDRO_TS_TIME_END"];
                    m_uiTinfo._m_uiT     = checkPoint["DENDRO_TS_TIME_CURRENT"];
                    m_uiTinfo._m_uiStep  = checkPoint["DENDRO_TS_STEP_CURRENT"]; 
                    m_uiTinfo._m_uiTh    = checkPoint["DENDRO_TS_TIME_STEP_SIZE"];
                    m_uiElementOrder     = checkPoint["DENDRO_TS_ELEMENT_ORDER"];
                    
                    bssn::BSSN_WAVELET_TOL=checkPoint["DENDRO_TS_WAVELET_TOLERANCE"];
                    bssn::BSSN_LOAD_IMB_TOL=checkPoint["DENDRO_TS_LOAD_IMB_TOLERANCE"];
                    
                    numVars=checkPoint["DENDRO_TS_NUM_VARS"];
                    activeCommSz=checkPoint["DENDRO_TS_ACTIVE_COMM_SZ"];
                    
                    m_uiBHLoc[0]=Point((double)checkPoint["DENDRO_BH1_X"],(double)checkPoint["DENDRO_BH1_Y"],(double)checkPoint["DENDRO_BH1_Z"]);
                    m_uiBHLoc[1]=Point((double)checkPoint["DENDRO_BH2_X"],(double)checkPoint["DENDRO_BH2_Y"],(double)checkPoint["DENDRO_BH2_Z"]);
                    restoreStep[restoreFileIndex]=m_uiTinfo._m_uiStep;
                }


            }

            par::Mpi_Allreduce(&restoreStatus,&restoreStatusGlobal,1,MPI_MAX,comm);
            if(restoreStatusGlobal==1) 
            {
                if(!rank)
                    std::cout<<"[BSSNCtx] : Restore step failed, restore file corrupted. "<<std::endl;
                MPI_Abort(comm,0);
            }


            MPI_Bcast(&m_uiTinfo,sizeof(ts::TSInfo),MPI_BYTE,0,comm);
            par::Mpi_Bcast(&bssn::BSSN_WAVELET_TOL,1,0,comm);
            par::Mpi_Bcast(&bssn::BSSN_LOAD_IMB_TOL,1,0,comm);

            par::Mpi_Bcast(&numVars,1,0,comm);
            par::Mpi_Bcast(&m_uiElementOrder,1,0,comm);
            par::Mpi_Bcast(&activeCommSz,1,0,comm);
            
            par::Mpi_Bcast(m_uiBHLoc,2,0,comm);
            bssn::BSSN_BH_LOC[0]=m_uiBHLoc[0];
            bssn::BSSN_BH_LOC[1]=m_uiBHLoc[1];

            if(activeCommSz>npes)
            {
                if(!rank)
                    std::cout<<" [BSSNCtx] : checkpoint file written from  a larger communicator than the current global comm. (i.e. communicator shrinking not allowed in the restore step. )"<<std::endl;
                
                MPI_Abort(comm,0);
            }



            bool isActive=(rank<activeCommSz);

            MPI_Comm newComm;
            par::splitComm2way(isActive,&newComm,comm);

            if(isActive) {

                int activeRank;
                int activeNpes;

                MPI_Comm_rank(newComm, &activeRank);
                MPI_Comm_size(newComm, &activeNpes);
                assert(activeNpes == activeCommSz);

                sprintf(fName, "%s_octree_%d_%d.oct", bssn::BSSN_CHKPT_FILE_PREFIX.c_str(),restoreFileIndex,activeRank);
                restoreStatus=io::checkpoint::readOctFromFile(fName, octree);
                assert(par::test::isUniqueAndSorted(octree, newComm));

            }

            par::Mpi_Allreduce(&restoreStatus,&restoreStatusGlobal,1,MPI_MAX,comm);
            if(restoreStatusGlobal==1) {

                if(!rank) std::cout<<"[BSSNCtx]: octree (*.oct) restore file is corrupted "<<std::endl;
                MPI_Abort(comm,0);
            }

            newMesh=new ot::Mesh(octree,1,m_uiElementOrder,activeCommSz,comm);
            newMesh->setDomainBounds(Point(bssn::BSSN_GRID_MIN_X,bssn::BSSN_GRID_MIN_Y,bssn::BSSN_GRID_MIN_Z), Point(bssn::BSSN_GRID_MAX_X, bssn::BSSN_GRID_MAX_Y,bssn::BSSN_GRID_MAX_Z));
            // no need to transfer data only to resize the contex variables. 
            this->grid_transfer( newMesh , false, false, false);
            this->update_app_vars();
            
            // only reads the evolution variables. 
            if(isActive) {

                int activeRank;
                int activeNpes;

                DendroScalar** inVec=NULL;
                m_uiEVar.Get2DArray(inVec,false);

                MPI_Comm_rank(newComm, &activeRank);
                MPI_Comm_size(newComm, &activeNpes);
                assert(activeNpes == activeCommSz);

                sprintf(fName,"%s_%d_%d.var",bssn::BSSN_CHKPT_FILE_PREFIX.c_str(),restoreFileIndex,activeRank);
                restoreStatus=io::checkpoint::readVecFromFile(fName,newMesh,inVec,bssn::BSSN_NUM_VARS);

                delete [] inVec;
            }

            MPI_Comm_free(&newComm);
            par::Mpi_Allreduce(&restoreStatus,&restoreStatusGlobal,1,MPI_MAX,comm);
            if(restoreStatusGlobal==1) {

                if(!rank) std::cout<<"[BSSNCtx]: varible (*.var) restore file currupted "<<std::endl;
                MPI_Abort(comm,0);
            }

            std::swap(m_uiMesh,newMesh);
            delete newMesh;
            
        unsigned int localSz=m_uiMesh->getNumLocalMeshElements();
        unsigned int totalElems=0;
        par::Mpi_Allreduce(&localSz, &totalElems ,1,MPI_SUM,comm);
        
        if(!rank) std::cout<<" checkpoint at step : "<<m_uiTinfo._m_uiStep<<"active Comm. sz: "<<activeCommSz<<" restore successful: "<<" restored mesh size: "<<totalElems<<std::endl;

        m_uiIsETSSynced = false;
        return 0;

    }

    int BSSNCtx::pre_timestep(DVec sIn)
    {

        return 0;
    }

    int BSSNCtx::pre_stage(DVec  sIn)
    {
        
        return 0;
    }

    int BSSNCtx::post_stage(DVec sIn)
    { 
        return 0;
    }

    int BSSNCtx::post_timestep(DVec sIn)
    {
        // we need to enforce constraint before computing the HAM and MOM_i constraints. 
        DendroScalar ** sVar;
        sIn.Get2DArray(sVar,false);

        for(unsigned int node=m_uiMesh->getNodeLocalBegin(); node< m_uiMesh->getNodeLocalEnd(); node++)
            enforce_bssn_constraints(sVar, node);

        delete [] sVar;

        return 0;
    }

    bool BSSNCtx::is_remesh()
    {
        bool isRefine = false;
        if(bssn::BSSN_ENABLE_BLOCK_ADAPTIVITY)
            return false;

        
        MPI_Comm comm = m_uiMesh->getMPIGlobalCommunicator();
        
        this->unzip(m_uiEVar,m_uiEUnzip[0],bssn::BSSN_ASYNC_COMM_K);

        DendroScalar** unzipVar;
        m_uiEUnzip[0].Get2DArray(unzipVar, false);

        unsigned int refineVarIds[bssn::BSSN_NUM_REFINE_VARS];
        for(unsigned int vIndex=0; vIndex<bssn::BSSN_NUM_REFINE_VARS; vIndex++)
            refineVarIds[vIndex]=bssn::BSSN_REFINE_VARIABLE_INDICES[vIndex];

        double wTol=bssn::BSSN_WAVELET_TOL;
        std::function<double(double,double,double,double*hx)> waveletTolFunc =[](double x,double y, double z,double*hx) {
            return bssn::computeWTolDCoords(x,y,z,hx);
        };
        

        if(bssn::BSSN_REFINEMENT_MODE == bssn::RefinementMode::WAMR) {
            isRefine = bssn::isReMeshWAMR(m_uiMesh,(const double **)unzipVar,refineVarIds,bssn::BSSN_NUM_REFINE_VARS,waveletTolFunc,bssn::BSSN_DENDRO_AMR_FAC);

        }else if(bssn::BSSN_REFINEMENT_MODE == bssn::RefinementMode::EH)
        {
            isRefine = bssn::isRemeshEH(m_uiMesh,(const double **)unzipVar,bssn::VAR::U_ALPHA,bssn::BSSN_EH_REFINE_VAL,bssn::BSSN_EH_COARSEN_VAL,true);

        }else if(bssn::BSSN_REFINEMENT_MODE == bssn::RefinementMode::EH_WAMR)
        {
            const bool isR1 = bssn::isReMeshWAMR(m_uiMesh,(const double **)unzipVar,refineVarIds,bssn::BSSN_NUM_REFINE_VARS,waveletTolFunc,bssn::BSSN_DENDRO_AMR_FAC);
            const bool isR2 = bssn::isRemeshEH(m_uiMesh,(const double **)unzipVar,bssn::VAR::U_ALPHA,bssn::BSSN_EH_REFINE_VAL,bssn::BSSN_EH_COARSEN_VAL,false);

            isRefine = (isR1 || isR2);
        }else if( bssn::BSSN_REFINEMENT_MODE == bssn::RefinementMode::BH_LOC)
        {
            isRefine = bssn::isRemeshBH(m_uiMesh,m_uiBHLoc);
        }

        delete [] unzipVar;
        
        return isRefine;

    }

    int BSSNCtx::update_app_vars() 
    {
        m_uiEVar = m_uiEvolutionVar[0];
        m_uiCVar = m_uiConstrainedVar[0];

        m_uiEUnzip[0] = m_uiEvolutionUnzipVar[0];
        m_uiEUnzip[1] = m_uiEvolutionUnzipVar[1];

        m_uiCUnzip[0] = m_uiConstraintUnzipVar[0];

        return 0;
        
    }

    DVec BSSNCtx::get_evolution_vars()
    {
        return m_uiEVar;
    }
    
    DVec BSSNCtx::get_constraint_vars()
    {
        return m_uiCVar;
    }

    DVec BSSNCtx::get_primitive_vars()
    {
        return m_uiPVar;
    }

    int BSSNCtx::terminal_output()
    {
        DendroScalar min=0, max=0;
        m_uiEVar.VecMinMax(m_uiMesh,min,max,bssn::VAR::U_ALPHA);

        if(m_uiMesh->isActive())
        {
            if(!(m_uiMesh->getMPIRank()))
                std::cout<<"[BSSNCtx]:  "<<bssn::BSSN_VAR_NAMES[bssn::VAR::U_ALPHA]<<" (min,max) : \t ( "<<min<<", "<<max<<" ) "<<std::endl;
        }

        return 0; 
    }

    unsigned int BSSNCtx::compute_lts_ts_offset()
    {

        const ot::Block* blkList     = m_uiMesh->getLocalBlockList().data();
        const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();
        const ot::TreeNode* pNodes   = m_uiMesh->getAllElements().data();

        const unsigned int ldiff     = 0;

        unsigned int lmin, lmax;
        m_uiMesh->computeMinMaxLevel(lmin,lmax);
        BSSN_LTS_TS_OFFSET=0;

        const double dt_min        = bssn::BSSN_CFL_FACTOR*((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))));
        const double dt_eta_fac    = dt_min*ETA_CONST;
        
        for(unsigned int blk = 0; blk < numBlocks; blk++)
        {
            const unsigned int blev      = pNodes[blkList[blk].getLocalElementBegin()].getLevel(); 
            ot::TreeNode blkNode         = blkList[blk].getBlockNode(); 
            const unsigned int szby2     = (1u<<(m_uiMaxDepth-blkNode.getLevel()-1));

            Point pt_oct(blkNode.minX() + szby2 , blkNode.minY() + szby2,  blkNode.minZ() + szby2);
            Point pt_domain;
            m_uiMesh->octCoordToDomainCoord(pt_oct,pt_domain);
            const double r_coord = pt_domain.abs();
            
            double eta;
            eta=ETA_CONST;
            if (r_coord >= ETA_R0) {
                eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
            }

            const double dt_eta              = dt_eta_fac*(1/eta);
            const double dt_cfl              = (1u<<(lmax-blev))*dt_min;

            const double dt_feasible         = std::min(dt_eta,dt_cfl);

            if(dt_feasible > dt_min)
            {
                unsigned int lts_offset    =  lmax- blev - std::floor(std::log2(dt_feasible/dt_min)) ;
                if(BSSN_LTS_TS_OFFSET  < lts_offset)
                    BSSN_LTS_TS_OFFSET = lts_offset;
            }    
            
            
        }

        unsigned int lts_offset_max=0;
        par::Mpi_Allreduce(&BSSN_LTS_TS_OFFSET,&lts_offset_max,1,MPI_MAX,m_uiMesh->getMPIGlobalCommunicator());
        BSSN_LTS_TS_OFFSET=lts_offset_max;

        return BSSN_LTS_TS_OFFSET;

    }

    unsigned int BSSNCtx::getBlkTimestepFac(unsigned int blev, unsigned int lmin, unsigned int lmax)
    {   
        const unsigned int ldiff=BSSN_LTS_TS_OFFSET;
        if((lmax-blev) <=ldiff)
            return 1;
        else
        {
            return 1u<<(lmax-blev-ldiff);
        }

    }

    void BSSNCtx::evolve_bh_loc(DVec sIn, double dt)
    {
        #ifdef BSSN_EXTRACT_BH_LOCATIONS

            m_uiMesh->readFromGhostBegin(sIn.GetVecArray()+ VAR::U_BETA0 * m_uiMesh->getDegOfFreedom(),3);
            m_uiMesh->readFromGhostEnd(sIn.GetVecArray()  + VAR::U_BETA0 * m_uiMesh->getDegOfFreedom(),3);
            Point bhLoc[2];
            DendroScalar ** evar = new DendroScalar*[bssn::BSSN_NUM_VARS];
            sIn.Get2DArray(evar,true);
            bssn::computeBHLocations((const ot::Mesh *)m_uiMesh,m_uiBHLoc,bhLoc,evar,dt);
            // if(!m_uiMesh->getMPIRankGlobal())
            // {
            //     std::cout<<"bh0 "<<bhLoc[0]<<std::endl;
            //     std::cout<<"bh1 "<<bhLoc[1]<<std::endl;

            // }
            m_uiBHLoc[0] = bhLoc[0];
            m_uiBHLoc[1] = bhLoc[1];

            bssn::BSSN_BH_LOC[0]=m_uiBHLoc[0];
            bssn::BSSN_BH_LOC[1]=m_uiBHLoc[1];

            delete [] evar;

            // old bh location extractor. 
            #if 0
                DendroScalar ** evar = new DendroScalar*[bssn::BSSN_NUM_VARS];
                Point bhLoc[2];
                sIn.Get2DArray(evar,true);
                bssn::extractBHCoords((const ot::Mesh *)m_uiMesh,(const DendroScalar*)evar[BHLOC::EXTRACTION_VAR_ID],BHLOC::EXTRACTION_TOL,(const Point *) m_uiBHLoc,2,(Point*)bhLoc);
                
                m_uiBHLoc[0] = bhLoc[0];
                m_uiBHLoc[1] = bhLoc[1];

                delete [] evar;
            #endif
            
        #endif

        return;

    }
    
    void BSSNCtx::lts_smooth(DVec sIn, LTS_SMOOTH_MODE mode)
    {

        if(mode == LTS_SMOOTH_MODE::KO)
        {

            this->unzip(sIn, m_uiEUnzip[0] , bssn::BSSN_ASYNC_COMM_K);
        
            DendroScalar ** uZipVars;
            DendroScalar **  unzipVarsRHS; 
        
            m_uiEUnzip[0].Get2DArray(uZipVars , false);
            m_uiEUnzip[1].Get2DArray(unzipVarsRHS ,false);

            const ot::Block* blkList = m_uiMesh->getLocalBlockList().data();
            const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();

            unsigned int sz[3];
            const Point pt_min(bssn::BSSN_COMPD_MIN[0],bssn::BSSN_COMPD_MIN[1],bssn::BSSN_COMPD_MIN[2]);
            const Point pt_max(bssn::BSSN_COMPD_MAX[0],bssn::BSSN_COMPD_MAX[1],bssn::BSSN_COMPD_MAX[2]);
            double ptmin[3], ptmax[3];

            mem::memory_pool<double>* __mem_pool = &BSSN_MEM_POOL;

            for(unsigned int blk=0; blk<numBlocks; blk++)
            {
                const unsigned int offset=blkList[blk].getOffset();
                sz[0]=blkList[blk].getAllocationSzX();
                sz[1]=blkList[blk].getAllocationSzY();
                sz[2]=blkList[blk].getAllocationSzZ();

                const unsigned int bflag=blkList[blk].getBlkNodeFlag();

                const double dx=blkList[blk].computeDx(pt_min,pt_max);
                const double dy=blkList[blk].computeDy(pt_min,pt_max);
                const double dz=blkList[blk].computeDz(pt_min,pt_max);

                ptmin[0]=GRIDX_TO_X(blkList[blk].getBlockNode().minX())-3*dx;
                ptmin[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().minY())-3*dy;
                ptmin[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ())-3*dz;

                ptmax[0]=GRIDX_TO_X(blkList[blk].getBlockNode().maxX())+3*dx;
                ptmax[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().maxY())+3*dy;
                ptmax[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ())+3*dz;

                const double *alpha = &uZipVars[VAR::U_ALPHA][offset];
                const double *chi = &uZipVars[VAR::U_CHI][offset];
                const double *K = &uZipVars[VAR::U_K][offset];
                const double *gt0 = &uZipVars[VAR::U_SYMGT0][offset];
                const double *gt1 = &uZipVars[VAR::U_SYMGT1][offset];
                const double *gt2 = &uZipVars[VAR::U_SYMGT2][offset];
                const double *gt3 = &uZipVars[VAR::U_SYMGT3][offset];
                const double *gt4 = &uZipVars[VAR::U_SYMGT4][offset];
                const double *gt5 = &uZipVars[VAR::U_SYMGT5][offset];
                const double *beta0 = &uZipVars[VAR::U_BETA0][offset];
                const double *beta1 = &uZipVars[VAR::U_BETA1][offset];
                const double *beta2 = &uZipVars[VAR::U_BETA2][offset];
                const double *At0 = &uZipVars[VAR::U_SYMAT0][offset];
                const double *At1 = &uZipVars[VAR::U_SYMAT1][offset];
                const double *At2 = &uZipVars[VAR::U_SYMAT2][offset];
                const double *At3 = &uZipVars[VAR::U_SYMAT3][offset];
                const double *At4 = &uZipVars[VAR::U_SYMAT4][offset];
                const double *At5 = &uZipVars[VAR::U_SYMAT5][offset];
                const double *Gt0 = &uZipVars[VAR::U_GT0][offset];
                const double *Gt1 = &uZipVars[VAR::U_GT1][offset];
                const double *Gt2 = &uZipVars[VAR::U_GT2][offset];
                const double *B0 = &uZipVars[VAR::U_B0][offset];
                const double *B1 = &uZipVars[VAR::U_B1][offset];
                const double *B2 = &uZipVars[VAR::U_B2][offset];

                double *a_rhs = &unzipVarsRHS[VAR::U_ALPHA][offset];
                double *chi_rhs = &unzipVarsRHS[VAR::U_CHI][offset];
                double *K_rhs = &unzipVarsRHS[VAR::U_K][offset];
                double *gt_rhs00 = &unzipVarsRHS[VAR::U_SYMGT0][offset];
                double *gt_rhs01 = &unzipVarsRHS[VAR::U_SYMGT1][offset];
                double *gt_rhs02 = &unzipVarsRHS[VAR::U_SYMGT2][offset];
                double *gt_rhs11 = &unzipVarsRHS[VAR::U_SYMGT3][offset];
                double *gt_rhs12 = &unzipVarsRHS[VAR::U_SYMGT4][offset];
                double *gt_rhs22 = &unzipVarsRHS[VAR::U_SYMGT5][offset];
                double *b_rhs0 = &unzipVarsRHS[VAR::U_BETA0][offset];
                double *b_rhs1 = &unzipVarsRHS[VAR::U_BETA1][offset];
                double *b_rhs2 = &unzipVarsRHS[VAR::U_BETA2][offset];
                double *At_rhs00 = &unzipVarsRHS[VAR::U_SYMAT0][offset];
                double *At_rhs01 = &unzipVarsRHS[VAR::U_SYMAT1][offset];
                double *At_rhs02 = &unzipVarsRHS[VAR::U_SYMAT2][offset];
                double *At_rhs11 = &unzipVarsRHS[VAR::U_SYMAT3][offset];
                double *At_rhs12 = &unzipVarsRHS[VAR::U_SYMAT4][offset];
                double *At_rhs22 = &unzipVarsRHS[VAR::U_SYMAT5][offset];
                double *Gt_rhs0 = &unzipVarsRHS[VAR::U_GT0][offset];
                double *Gt_rhs1 = &unzipVarsRHS[VAR::U_GT1][offset];
                double *Gt_rhs2 = &unzipVarsRHS[VAR::U_GT2][offset];
                double *B_rhs0 = &unzipVarsRHS[VAR::U_B0][offset];
                double *B_rhs1 = &unzipVarsRHS[VAR::U_B1][offset];
                double *B_rhs2 = &unzipVarsRHS[VAR::U_B2][offset];

                const unsigned int nx = sz[0];
                const unsigned int ny = sz[1];
                const unsigned int nz = sz[2];
                
                double hx = (ptmax[0] - ptmin[0]) / (nx - 1);
                double hy = (ptmax[1] - ptmin[1]) / (ny - 1);
                double hz = (ptmax[2] - ptmin[2]) / (nz - 1);

                // const unsigned int n =sz[0]*sz[1]*sz[2];
                // #include "bssnrhs_lts_smooth_ko.h"

                // const  double sigma = KO_DISS_SIGMA;


                // for (unsigned int k = 3; k < nz-3; k++) {
                //     for (unsigned int j = 3; j < ny-3; j++) {
                //         for (unsigned int i = 3; i < nx-3; i++) {
                //             const unsigned int pp = i + nx*(j + ny*k);

                //             a_rhs[pp]  = alpha[pp] + sigma * (grad_0_alpha[pp] + grad_1_alpha[pp] + grad_2_alpha[pp]);
                //             b_rhs0[pp] = beta0[pp] + sigma * (grad_0_beta0[pp] + grad_1_beta0[pp] + grad_2_beta0[pp]);
                //             b_rhs1[pp] = beta1[pp] + sigma * (grad_0_beta1[pp] + grad_1_beta1[pp] + grad_2_beta1[pp]);
                //             b_rhs2[pp] = beta2[pp] + sigma * (grad_0_beta2[pp] + grad_1_beta2[pp] + grad_2_beta2[pp]);

                //             gt_rhs00[pp] = gt0[pp] + sigma * (grad_0_gt0[pp] + grad_1_gt0[pp] + grad_2_gt0[pp]);
                //             gt_rhs01[pp] = gt1[pp] + sigma * (grad_0_gt1[pp] + grad_1_gt1[pp] + grad_2_gt1[pp]);
                //             gt_rhs02[pp] = gt2[pp] + sigma * (grad_0_gt2[pp] + grad_1_gt2[pp] + grad_2_gt2[pp]);
                //             gt_rhs11[pp] = gt3[pp] + sigma * (grad_0_gt3[pp] + grad_1_gt3[pp] + grad_2_gt3[pp]);
                //             gt_rhs12[pp] = gt4[pp] + sigma * (grad_0_gt4[pp] + grad_1_gt4[pp] + grad_2_gt4[pp]);
                //             gt_rhs22[pp] = gt5[pp] + sigma * (grad_0_gt5[pp] + grad_1_gt5[pp] + grad_2_gt5[pp]);

                //             chi_rhs[pp]  = chi[pp] + sigma * (grad_0_chi[pp] + grad_1_chi[pp] + grad_2_chi[pp]);

                //             At_rhs00[pp] = At0[pp] + sigma * (grad_0_At0[pp] + grad_1_At0[pp] + grad_2_At0[pp]);
                //             At_rhs01[pp] = At1[pp] + sigma * (grad_0_At1[pp] + grad_1_At1[pp] + grad_2_At1[pp]);
                //             At_rhs02[pp] = At2[pp] + sigma * (grad_0_At2[pp] + grad_1_At2[pp] + grad_2_At2[pp]);
                //             At_rhs11[pp] = At3[pp] + sigma * (grad_0_At3[pp] + grad_1_At3[pp] + grad_2_At3[pp]);
                //             At_rhs12[pp] = At4[pp] + sigma * (grad_0_At4[pp] + grad_1_At4[pp] + grad_2_At4[pp]);
                //             At_rhs22[pp] = At5[pp] + sigma * (grad_0_At5[pp] + grad_1_At5[pp] + grad_2_At5[pp]);

                //             K_rhs[pp] = K[pp] + sigma * (grad_0_K[pp] + grad_1_K[pp] + grad_2_K[pp]);

                //             Gt_rhs0[pp] = Gt0[pp] + sigma * (grad_0_Gt0[pp] + grad_1_Gt0[pp] + grad_2_Gt0[pp]);
                //             Gt_rhs1[pp] = Gt1[pp] + sigma * (grad_0_Gt1[pp] + grad_1_Gt1[pp] + grad_2_Gt1[pp]);
                //             Gt_rhs2[pp] = Gt2[pp] + sigma * (grad_0_Gt2[pp] + grad_1_Gt2[pp] + grad_2_Gt2[pp]);

                //             B_rhs0[pp] = B0[pp] + sigma * (grad_0_B0[pp] + grad_1_B0[pp] + grad_2_B0[pp]);
                //             B_rhs1[pp] = B1[pp] + sigma * (grad_0_B1[pp] + grad_1_B1[pp] + grad_2_B1[pp]);
                //             B_rhs2[pp] = B2[pp] + sigma * (grad_0_B2[pp] + grad_1_B2[pp] + grad_2_B2[pp]);
                //         }
                //     }
                // }

                
            }

            this->zip(m_uiEUnzip[1], sIn, bssn::BSSN_ASYNC_COMM_K);
            delete [] uZipVars;
            delete [] unzipVarsRHS;


        }

        



    }


}// end of namespace bssn. 