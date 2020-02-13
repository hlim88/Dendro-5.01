// Created by milinda on 11/19/18.
/**
 *@file oda.cpp
 *@brief The class that manages the octree mesh that support FEM computations. Note that this file is a refactored version from Dendro4
 *which is changed to support Dendro-5.0  Currently we use mesh based FEM computations but in future we will move towards the mesh - free FEM computation methods.
 *Distributed Array (DA) methods currently uses ot::Mesh class to provide interface to write FEM loops.
  @author Milinda Fernando
  date: 11/19/2018
 * */

#include "oda.h"

namespace ot
{

    /**@brief: Constructor for the DA data structures
      * @param [in] in : input octree, need to be 2:1 balanced unique sorted octree.
      * @param [in] comm: MPI global communicator for mesh generation.
      * @param [in] order: order of the element.
      * @param [in] grainSz: Number of suggested elements per processor,
      * @param [in] sfc_tol: SFC partitioning tolerance,
     * */
    DA::DA(std::vector<ot::TreeNode> &balOct,MPI_Comm comm,unsigned int order, unsigned int grainSz,double sfc_tol){

        assert(par::test::isUniqueAndSorted(balOct,comm));

        // generates the mesh.
        m_uiMesh=new ot::Mesh(balOct,1,order,comm, false, grainSz,sfc_tol);
        balOct.clear();


#ifdef DEBUG_ODA
        std::vector<ot::TreeNode> independent;
        std::vector<ot::TreeNode> w_dependent;
        std::vector<ot::TreeNode> w_bdy;
#endif

        // perform element flags.
        ot::computeODAFlags(m_uiMesh,m_uiOctantFlags);
        // compute the local node to global node mapping.
        ot::computeLocalToGlobalNodalMap(m_uiMesh,m_uiLocalToGlobalNodalMap,m_uiGlobalNodeSz);



#ifdef DEBUG_ODA
        for(init<DA_FLAGS::INDEPENDENT>();curr()<end<DA_FLAGS::INDEPENDENT>();next<DA_FLAGS::INDEPENDENT>())
            independent.push_back(allElements[curr()]);

        for(init<DA_FLAGS::W_DEPENDENT>();curr()<end<DA_FLAGS::W_DEPENDENT>();next<DA_FLAGS::W_DEPENDENT>())
            w_dependent.push_back(allElements[curr()]);

        for(init<DA_FLAGS::W_BOUNDARY>();curr()<end<DA_FLAGS::W_BOUNDARY>();next<DA_FLAGS::W_BOUNDARY>())
            w_bdy.push_back(allElements[curr()]);//std::cout<<"loop ele: "<<curr()<<" : "<<allElements[curr()]<<std::endl;//


        io::vtk::oct2vtu(&(*(independent.begin())),independent.size(),"independent",comm);
        io::vtk::oct2vtu(&(*(w_dependent.begin())),w_dependent.size(),"w_dependent",comm);
        io::vtk::oct2vtu(&(*(w_bdy.begin())),w_bdy.size(),"w_bdy",comm);
#endif

        m_uiTotalElementSz=m_uiMesh->getAllElements().size();
        m_uiLocalElementSz=m_uiMesh->getNumLocalMeshElements();

        m_uiTotalNodalSz=m_uiMesh->getDegOfFreedom();
        m_uiLocalNodalSz=m_uiMesh->getNumLocalMeshNodes();


        m_uiLoopInfo.currentIndex=0;
        m_uiLoopInfo.indexBegin=0;
        m_uiLoopInfo.indexEnd=0;

        m_uiOctreeLowerBound[0]=0;
        m_uiOctreeLowerBound[1]=0;
        m_uiOctreeLowerBound[2]=0;

        m_uiOctreeUpperBound[0]=1u<<(m_uiMaxDepth);
        m_uiOctreeUpperBound[1]=1u<<(m_uiMaxDepth);
        m_uiOctreeUpperBound[2]=1u<<(m_uiMaxDepth);


        m_uiCommTag=0;
        m_uiMPIContexts.clear();


    }

    DA::DA(ot::Mesh *pMesh)
    {
        m_uiMesh=pMesh;

#ifdef DEBUG_ODA
        std::vector<ot::TreeNode> independent;
        std::vector<ot::TreeNode> w_dependent;
        std::vector<ot::TreeNode> w_bdy;
#endif

        // perform element flags.
        ot::computeODAFlags(m_uiMesh,m_uiOctantFlags);
        // compute the local node to global node mapping.
        ot::computeLocalToGlobalNodalMap(m_uiMesh,m_uiLocalToGlobalNodalMap,m_uiGlobalNodeSz);



#ifdef DEBUG_ODA
        for(init<DA_FLAGS::INDEPENDENT>();curr()<end<DA_FLAGS::INDEPENDENT>();next<DA_FLAGS::INDEPENDENT>())
            independent.push_back(allElements[curr()]);

        for(init<DA_FLAGS::W_DEPENDENT>();curr()<end<DA_FLAGS::W_DEPENDENT>();next<DA_FLAGS::W_DEPENDENT>())
            w_dependent.push_back(allElements[curr()]);

        for(init<DA_FLAGS::W_BOUNDARY>();curr()<end<DA_FLAGS::W_BOUNDARY>();next<DA_FLAGS::W_BOUNDARY>())
            w_bdy.push_back(allElements[curr()]);//std::cout<<"loop ele: "<<curr()<<" : "<<allElements[curr()]<<std::endl;//


        io::vtk::oct2vtu(&(*(independent.begin())),independent.size(),"independent",comm);
        io::vtk::oct2vtu(&(*(w_dependent.begin())),w_dependent.size(),"w_dependent",comm);
        io::vtk::oct2vtu(&(*(w_bdy.begin())),w_bdy.size(),"w_bdy",comm);
#endif

        m_uiTotalElementSz=m_uiMesh->getAllElements().size();
        m_uiLocalElementSz=m_uiMesh->getNumLocalMeshElements();

        m_uiTotalNodalSz=m_uiMesh->getDegOfFreedom();
        m_uiLocalNodalSz=m_uiMesh->getNumLocalMeshNodes();


        m_uiLoopInfo.currentIndex=0;
        m_uiLoopInfo.indexBegin=0;
        m_uiLoopInfo.indexEnd=0;

        m_uiOctreeLowerBound[0]=0;
        m_uiOctreeLowerBound[1]=0;
        m_uiOctreeLowerBound[2]=0;

        m_uiOctreeUpperBound[0]=1u<<(m_uiMaxDepth);
        m_uiOctreeUpperBound[1]=1u<<(m_uiMaxDepth);
        m_uiOctreeUpperBound[2]=1u<<(m_uiMaxDepth);


        m_uiCommTag=0;
        m_uiMPIContexts.clear();
    }


    DA::~DA()
    {
        m_uiMPIContexts.clear();
        m_uiOctantFlags.clear();
        delete m_uiMesh;
        m_uiMesh=NULL;
    }

    bool DA::isBoundaryOctant(unsigned int eleID) const
    {
        return (binOp::getBit(m_uiOctantFlags[eleID],ot::DA_FLAGS::W_BOUNDARY)==1) ? true : false;
    }

    void DA::getElementalCoords(unsigned int eleID, double* coords) const
    {
        m_uiMesh->getElementCoordinates(eleID,coords);
    }

    void DA::getOctreeBoundaryNodeIndices(std::vector<unsigned int >& bdyIndex,std::vector<double>& coords,bool isGhosted)
    {
        bdyIndex.clear();
        coords.clear();
        ot::TreeNode tmpOct;
        const ot::TreeNode* allElements=&(*(m_uiMesh->getAllElements().begin()));
        const unsigned int * e2n=&(*(m_uiMesh->getE2NMapping().begin()));
        const unsigned int * e2n_dg=&(*(m_uiMesh->getE2NMapping_DG().begin()));
        const unsigned int * cgToDg=&(*(m_uiMesh->getCG2DGMap().begin()));

        const unsigned int octreeMin=0;
        const unsigned int octreeMax=1u<<(m_uiMaxDepth);

        const unsigned int nx=m_uiMesh->getElementOrder()+1;
        const unsigned int ny=m_uiMesh->getElementOrder()+1;
        const unsigned int nz=m_uiMesh->getElementOrder()+1;

        unsigned int eleOrder=m_uiMesh->getElementOrder();

        const unsigned int nPe=m_uiMesh->getNumNodesPerElement();
        const unsigned int nodeLocalBegin=m_uiMesh->getNodeLocalBegin();
        const unsigned int nodeLocalEnd=m_uiMesh->getNodeLocalEnd();

        unsigned int ele,lookup;



        for(init<ot::DA_FLAGS::W_BOUNDARY>();curr()<end<ot::DA_FLAGS::W_BOUNDARY>();next<ot::DA_FLAGS::W_BOUNDARY>())
        {
            ele=curr();
            tmpOct=allElements[ele];
            //std::cout<<"bdy element : "<<tmpOct<<std::endl;

            if(tmpOct.minX()==octreeMin)
            { // x_min boundary

                for(unsigned int k=0;k<nz;k++)
                   for(unsigned int j=0;j<ny;j++)
                   {
                       lookup=e2n[ele*nPe + k*ny*nx + j*nx + 0];
                       if(lookup>=nodeLocalBegin && lookup<nodeLocalEnd)
                           bdyIndex.push_back(lookup);
                   }

            }

            if(tmpOct.minY()==octreeMin)
            { // y_min boundary

                for(unsigned int k=0;k<nz;k++)
                    for(unsigned int i=0;i<nx;i++)
                    {
                        lookup=e2n[ele*nPe + k*ny*nx + 0*nx + i];
                        if(lookup>=nodeLocalBegin && lookup<nodeLocalEnd)
                            bdyIndex.push_back(lookup);
                    }

            }

            if(tmpOct.minZ()==octreeMin)
            { // z_min boundary

                for(unsigned int j=0;j<ny;j++)
                    for(unsigned int i=0;i<nx;i++)
                    {
                        lookup=e2n[ele*nPe + 0*ny*nx + j*nx + i];
                        if(lookup>=nodeLocalBegin && lookup<nodeLocalEnd)
                            bdyIndex.push_back(lookup);
                    }

            }

            if(tmpOct.maxX()==octreeMax)
            { // x_max boundary

                for(unsigned int k=0;k<nz;k++)
                    for(unsigned int j=0;j<ny;j++)
                    {
                        lookup=e2n[ele*nPe + k*ny*nx + j*nx + eleOrder];
                        if(lookup>=nodeLocalBegin && lookup<nodeLocalEnd)
                            bdyIndex.push_back(lookup);
                    }

            }

            if(tmpOct.maxY()==octreeMax)
            { // y_max boundary

                for(unsigned int k=0;k<nz;k++)
                    for(unsigned int i=0;i<nx;i++)
                    {
                        lookup=e2n[ele*nPe + k*ny*nx + eleOrder*nx + i];
                        if(lookup>=nodeLocalBegin && lookup<nodeLocalEnd)
                            bdyIndex.push_back(lookup);
                    }

            }

            if(tmpOct.maxZ()==octreeMax)
            { // z_max boundary

                for(unsigned int j=0;j<ny;j++)
                    for(unsigned int i=0;i<nx;i++)
                    {
                        lookup=e2n[ele*nPe + eleOrder*ny*nx + j*nx + i];
                        if(lookup>=nodeLocalBegin && lookup<nodeLocalEnd)
                            bdyIndex.push_back(lookup);
                    }

            }

        }


        std::sort(bdyIndex.begin(),bdyIndex.end());
        bdyIndex.erase( std::unique( bdyIndex.begin(), bdyIndex.end() ), bdyIndex.end());

        coords.resize(m_uiDim*bdyIndex.size());

        unsigned int ownerID,ii_x,jj_y,kk_z;
        double hx;


        for(unsigned int i=0;i<bdyIndex.size();i++)
        {
            lookup=cgToDg[bdyIndex[i]];
            m_uiMesh->dg2eijk(lookup,ownerID,ii_x,jj_y,kk_z);
            tmpOct=allElements[ownerID];
            hx=(tmpOct.maxX()-tmpOct.minX())/((double)eleOrder);

            coords[i*m_uiDim+0]=tmpOct.minX() + ii_x*hx;
            coords[i*m_uiDim+1]=tmpOct.minY() + jj_y*hx;
            coords[i*m_uiDim+2]=tmpOct.minZ() + kk_z*hx;

            //std::cout<<"x: "<<coords[i*m_uiDim+0]<<" y "<<coords[i*m_uiDim+1]<<" z: "<<coords[i*m_uiDim+2]<<std::endl;

        }

        if(!isGhosted)
        {
            for(unsigned int i=0;i<bdyIndex.size();i++)
                bdyIndex[i]=bdyIndex[i]-m_uiMesh->getNodeLocalBegin();
        }



    }

    int DA::getNodeIndices(DendroIntL* nodeIdx,unsigned int ele,bool isGhosted) const
    {

        const unsigned int * e2n=&(*(m_uiMesh->getE2NMapping().begin()));
        const unsigned int nPe=m_uiMesh->getNumNodesPerElement();

        for(unsigned int node=0;node<nPe;node++)
            nodeIdx[node]=e2n[ele*nPe+node];

        if(!isGhosted)
            for(unsigned int node=0;node<nPe;node++)
                nodeIdx[node]=nodeIdx[node]-m_uiMesh->getNodeLocalBegin();

        return 0;
    }

    int DA::getGlobalNodeIndices(DendroIntL* nodeIdx, unsigned int ele) const 
    {
        this->getNodeIndices(nodeIdx,ele,true);
        const unsigned int nPe=m_uiMesh->getNumNodesPerElement();
        for(unsigned int node=0;node<nPe;node++)
            nodeIdx[node]=m_uiLocalToGlobalNodalMap[nodeIdx[node]];
        
        return 0;
    }


    ot::DA* DA::remesh(const DA_FLAGS::Refine * flags, unsigned int sz,unsigned int grainSz,double ld_bal, unsigned int sfK)
    {
        const unsigned int localElementBegin=m_uiMesh->getElementLocalBegin();
        const unsigned int localElementEnd= m_uiMesh->getElementLocalEnd();
        bool isRemesh=false;
        bool isRemesh_g=false;
        MPI_Comm commGlobal=m_uiMesh->getMPIGlobalCommunicator();

        std::vector<unsigned int> octflags;
        octflags.resize(sz,OCT_NO_CHANGE);

        if(m_uiMesh->isActive())
        {
            const ot::TreeNode* allElements=&(*(m_uiMesh->getAllElements().begin()));
            //1. check to see if we need a remesh.
            for(unsigned int i=0;i<sz;i++)
            {
                unsigned int ele=i+localElementBegin;
                if(flags[i]==DA_FLAGS::Refine::DA_REFINE) {
                    if((allElements[ele].getLevel()+MAXDEAPTH_LEVEL_DIFF+1)>=m_uiMaxDepth)
                    {
                        octflags[i]=OCT_NO_CHANGE;
                    }else{
                        octflags[i]=OCT_SPLIT;
                        isRemesh=true;
                    }

                }else if(flags[i]==DA_FLAGS::Refine::DA_COARSEN) {

                    if(((i+NUM_CHILDREN-1)<sz)  && (allElements[ele].getParent() == allElements[ele+NUM_CHILDREN-1].getParent()) && (allElements[ele].getLevel()>0))
                    {
                        for(unsigned int child=0;child<NUM_CHILDREN;child++)
                            octflags[i+child]=OCT_COARSE;

                        isRemesh= true;
                        i+=(NUM_CHILDREN-1);
                    }
                }else{
                    octflags[i]=OCT_NO_CHANGE;
                }

            }

        }


        MPI_Allreduce(&isRemesh,&isRemesh_g,1,MPI_CXX_BOOL,MPI_LOR,commGlobal);
        // return null da since no need to remesh
        if(!isRemesh_g)
            return NULL;

        m_uiMesh->setOctreeRefineFlags(&(*(octflags.begin())),octflags.size());
        ot::Mesh* newMesh=m_uiMesh->ReMesh(grainSz,ld_bal,sfK);

        ot::DA* newDA= new DA(newMesh);

        return newDA;

    }


    // all the petsc functionalities goes below.
#ifdef BUILD_WITH_PETSC


    PetscErrorCode DA::petscCreateVector(Vec &local, bool isElemental, bool isGhosted, unsigned int dof)
    {
        unsigned int sz=0;
        MPI_Comm globalComm=this->getGlobalComm();
        if(!(m_uiMesh->isActive()))
        {
            local=NULL;

        }else {

            if(isElemental)
            {
                if(isGhosted)
                    sz=dof*m_uiTotalElementSz;
                else
                    sz=dof*m_uiLocalElementSz;

            }else {

                if(isGhosted)
                    sz=dof*m_uiTotalNodalSz;
                else
                    sz=dof*m_uiLocalNodalSz;
            }

        }

        VecCreate(globalComm,&local);
        PetscErrorCode status=VecSetSizes(local,sz,PETSC_DECIDE);

        if (this->getNpesAll() > 1) {
            VecSetType(local,VECMPI);
        } else {
            VecSetType(local,VECSEQ);
        }


        return status;


    }

    PetscErrorCode DA::createMatrix(Mat &M, MatType mtype, unsigned int dof)
    {



        if(m_uiMesh->isActive())
        {

            const unsigned int npesAll=m_uiMesh->getMPICommSizeGlobal();
            const unsigned int eleOrder=m_uiMesh->getElementOrder();
            // in linear cases, 53 can be generated with 27 + 27 -1(self) node.
            const unsigned int preAllocFactor=dof*(53*(eleOrder+1));

            // first determine the size ...
            unsigned int lSz = dof*(m_uiLocalNodalSz);
            MPI_Comm activeComm=m_uiMesh->getMPICommunicator();

            PetscBool isAij, isAijSeq, isAijPrl, isSuperLU, isSuperLU_Dist;
            PetscStrcmp(mtype,MATAIJ,&isAij);
            PetscStrcmp(mtype,MATSEQAIJ,&isAijSeq);
            PetscStrcmp(mtype,MATMPIAIJ,&isAijPrl);
            isSuperLU = PETSC_FALSE; // PetscStrcmp(mtype,MATSUPERLU,&isSuperLU);
            isSuperLU_Dist = PETSC_FALSE; // PetscStrcmp(mtype,MATSUPERLU_DIST,&isSuperLU_Dist);

            MatCreate(activeComm, &M);
            MatSetSizes(M, lSz,lSz, PETSC_DETERMINE, PETSC_DETERMINE);
            MatSetType(M,mtype);


            if(isAij || isAijSeq || isAijPrl || isSuperLU || isSuperLU_Dist) {
                if(npesAll > 1) {
                    MatMPIAIJSetPreallocation(M, preAllocFactor , PETSC_NULL, preAllocFactor , PETSC_NULL);
                }else {
                    MatSeqAIJSetPreallocation(M, preAllocFactor, PETSC_NULL);
                }
            }

        }



        return 0;
    }



    PetscErrorCode DA::petscNodalVecToGhostedNodal(const Vec& in,Vec& out,bool isAllocated,unsigned int dof)
    {

        if(!(m_uiMesh->isActive()))
            return 0 ;

        unsigned int status=0;
        if(!isAllocated)
            status=petscCreateVector(out,false,true,dof);

        PetscScalar * inArry=NULL;
        PetscScalar * outArry=NULL;

        VecGetArray(in,&inArry);
        VecGetArray(out,&outArry);


        for(unsigned int var=0;var<dof;var++)
            std::memcpy((outArry+var*m_uiTotalNodalSz+m_uiMesh->getNodeLocalBegin()),(inArry+var*m_uiLocalNodalSz),sizeof(PetscScalar)*(m_uiLocalNodalSz));

        VecRestoreArray(in,&inArry);
        VecRestoreArray(out,&outArry);

        return status;

    }


    PetscErrorCode DA::petscGhostedNodalToNodalVec(const Vec& gVec,Vec& local,bool isAllocated,unsigned int dof)
    {
        if(!(m_uiMesh->isActive()))
            return 0;

        unsigned int status=0;
        if(!isAllocated)
            status=petscCreateVector(local,false,false,dof);

        PetscScalar * gVecArry=NULL;
        PetscScalar * localArry=NULL;

        VecGetArray(gVec,&gVecArry);
        VecGetArray(local,&localArry);

        for(unsigned int var=0;var<dof;var++)
            std::memcpy((localArry + var*m_uiLocalNodalSz ),(gVecArry+var*m_uiTotalNodalSz+m_uiMesh->getNodeLocalBegin()),sizeof(PetscScalar)*(m_uiLocalNodalSz));

        VecRestoreArray(gVec,&gVecArry);
        VecRestoreArray(local,&localArry);

        return status;

    }


    void DA::petscReadFromGhostBegin(PetscScalar* vecArry, unsigned int dof)
    {
        if(!m_uiMesh->isActive())
            return;

        readFromGhostBegin(vecArry,dof);

        return;

    }

    void DA::petscReadFromGhostEnd(PetscScalar* vecArry, unsigned int dof)
    {
        if(!m_uiMesh->isActive())
            return;

        readFromGhostEnd(vecArry,dof);

        return;

    }


    void DA::petscVecTopvtu(const Vec& local, const char * fPrefix,char** nodalVarNames,bool isElemental,bool isGhosted,unsigned int dof)
    {

        PetscScalar *arry=NULL;
        VecGetArray(local,&arry);

        vecTopvtu(arry,fPrefix,nodalVarNames,isElemental,isGhosted,dof);

        VecRestoreArray(local,&arry);

    }


    PetscErrorCode DA::petscSetValuesInMatrix(Mat mat, std::vector<ot::MatRecord>& records,unsigned int dof, InsertMode mode)
    {

        if(records.empty())
            return 0;

        // assembly code based on Dendro4
        std::vector<PetscScalar > values;
        std::vector<PetscInt> colIndices;
        
        const unsigned int nodeLocalBegin=m_uiMesh->getNodeLocalBegin();
        const unsigned int nodeLocalEnd=m_uiMesh->getNodeLocalEnd();
        
        //Can make it more efficient later.
        if (!records.empty())
        {
            ///TODO: change this to more efficient one later. 
            std::vector<ot::MatRecord> tmpRecords;
            
            for(unsigned int i=0;i<records.size();i++)
            {
                if(records[i].getRowID()>=nodeLocalBegin && records[i].getRowID()<nodeLocalEnd)
                    tmpRecords.push_back(records[i]);
            }
            
            std::swap(tmpRecords,records);
            tmpRecords.clear();
            
            //Sort Order: row first, col next, val last
            std::sort(records.begin(), records.end());
            
            unsigned int currRecord = 0;

            while (currRecord < (records.size() - 1))
            {
                
                values.push_back(records[currRecord].getMatVal());
                colIndices.push_back(static_cast<PetscInt>((records[currRecord].getColDim()) + dof * m_uiLocalToGlobalNodalMap[records[currRecord].getColID()]));

                if ((records[currRecord].getRowID() != records[currRecord + 1].getRowID()) ||  (records[currRecord].getRowDim() != records[currRecord + 1].getRowDim()))
                {
                    PetscInt rowId = static_cast<PetscInt>((records[currRecord].getRowDim()) + dof * m_uiLocalToGlobalNodalMap[records[currRecord].getRowID()] );
                    MatSetValues(mat, 1, &rowId, colIndices.size(), (&(*colIndices.begin())),  (&(*values.begin())), mode);

                    colIndices.clear();
                    values.clear();
                }
                currRecord++;
            } //end while

            PetscInt rowId = static_cast<PetscInt>((records[currRecord].getRowDim()) +  dof * m_uiLocalToGlobalNodalMap[records[currRecord].getRowID()]);
            if (values.empty())
            {
                //Last row is different from the previous row
                PetscInt colId = static_cast<PetscInt>((records[currRecord].getColDim()) + dof * m_uiLocalToGlobalNodalMap[records[currRecord].getColID()]);
                PetscScalar value = records[currRecord].getMatVal();
                MatSetValues(mat, 1, &rowId, 1, &colId, &value, mode);
            }
            else
            {
                //Last row is same as the previous row
                values.push_back(records[currRecord].getMatVal());
                colIndices.push_back(static_cast<PetscInt>(static_cast<PetscInt>((records[currRecord].getColDim()) + dof * m_uiLocalToGlobalNodalMap[records[currRecord].getColID()])));
                MatSetValues(mat, 1, &rowId, colIndices.size(), (&(*colIndices.begin())), (&(*values.begin())), mode);
                colIndices.clear();
                values.clear();
            }
            records.clear();
        } // records not empty

        return 0;

    }

    
    PetscErrorCode DA::petscChangeVecToMatBased(Vec& v1,bool isElemental,bool isGhosted, unsigned int dof)
    {
        Vec tmp;
        petscCreateVector(tmp,isElemental,isGhosted,dof);
        unsigned int sz;
        if(isElemental)
        {
            if(isGhosted)
                sz=m_uiTotalElementSz;
            else
                sz=m_uiLocalElementSz;

        }else {

            if(isGhosted)
                sz=m_uiTotalNodalSz;
            else
                sz=m_uiLocalNodalSz;
        }
        
        PetscScalar * tmpArry=NULL;
        PetscScalar * v1Arry=NULL;

        VecGetArray(tmp,&tmpArry);
        VecGetArray(v1,&v1Arry);
        
        for(unsigned int node=0;node<sz;node++)
        {
            for(unsigned int var=0;var<dof;var++)
            {
                tmpArry[dof*node+var]=v1Arry[var*sz+node];
            }
        }
        
        VecRestoreArray(tmp,&tmpArry);
        VecRestoreArray(v1,&v1Arry);
        
       
        std::swap(tmp,v1);
        VecDestroy(&tmp);
        
        return 0;
    }

    
    
    PetscErrorCode DA::petscChangeVecToMatFree(Vec& v1,bool isElemental,bool isGhosted,unsigned int dof)
    {
        Vec tmp;
        petscCreateVector(tmp,isElemental,isGhosted,dof);
        unsigned int sz;
        if(isElemental)
        {
            if(isGhosted)
                sz=m_uiTotalElementSz;
            else
                sz=m_uiLocalElementSz;

        }else {

            if(isGhosted)
                sz=m_uiTotalNodalSz;
            else
                sz=m_uiLocalNodalSz;
        }
        
        PetscScalar * tmpArry=NULL;
        PetscScalar * v1Arry=NULL;

        VecGetArray(tmp,&tmpArry);
        VecGetArray(v1,&v1Arry);
        
        for(unsigned int node=0;node<sz;node++)
        {
            for(unsigned int var=0;var<dof;var++)
            {
                tmpArry[var*sz+node]=v1Arry[dof*node+var];
            }
        }
        
        VecRestoreArray(tmp,&tmpArry);
        VecRestoreArray(v1,&v1Arry);
        
       
        std::swap(tmp,v1);
        VecDestroy(&tmp);
        
        return 0;
    }
    
    
    
#endif

}



