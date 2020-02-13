//
//Origianlly created by milinda on 7/25/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Header file for the GR simulation.
*/
// Modified by Hyun Lim on Mar.5.2018 

#include "gr.h"
#include "grUtils.h"
#include "mpi.h"
#include "TreeNode.h"
#include "mesh.h"
#include <vector>
#include <iostream>
#include "rk45CCZ4.h"
#include "rk4CCZ4.h"
#include "rk3CCZ4.h"
#include "octUtils.h"

int main (int argc, char** argv)
{


    if(argc<2)
        std::cout<<"Usage: "<<argv[0]<<" paramFile"<<std::endl;

    MPI_Init(&argc,&argv);
    MPI_Comm comm=MPI_COMM_WORLD;

    int rank,npes;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);

    ccz4::timer::initFlops();

    ccz4::timer::total_runtime.start();


    //1 . read the parameter file.
    if(!rank) std::cout<<" reading parameter file :"<<argv[1]<<std::endl;
    ccz4::readParamFile(argv[1],comm);



    if(rank==1|| npes==1)
    {
        std::cout<<"parameters read: "<<std::endl;

        std::cout<<YLW<<"\tnpes :"<<npes<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_DIM :"<<ccz4::CCZ4_DIM<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_IO_OUTPUT_FREQ :"<<ccz4::CCZ4_IO_OUTPUT_FREQ<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_REMESH_TEST_FREQ :"<<ccz4::CCZ4_REMESH_TEST_FREQ<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_CHECKPT_FREQ :"<<ccz4::CCZ4_CHECKPT_FREQ<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_RESTORE_SOLVER :"<<ccz4::CCZ4_RESTORE_SOLVER<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_ENABLE_BLOCK_ADAPTIVITY :"<<ccz4::CCZ4_ENABLE_BLOCK_ADAPTIVITY<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_VTU_FILE_PREFIX :"<<ccz4::CCZ4_VTU_FILE_PREFIX<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_CHKPT_FILE_PREFIX :"<<ccz4::CCZ4_CHKPT_FILE_PREFIX<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_PROFILE_FILE_PREFIX :"<<ccz4::CCZ4_PROFILE_FILE_PREFIX<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_IO_OUTPUT_GAP :"<<ccz4::CCZ4_IO_OUTPUT_GAP<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_DENDRO_GRAIN_SZ :"<<ccz4::CCZ4_DENDRO_GRAIN_SZ<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_ASYNC_COMM_K :"<<ccz4::CCZ4_ASYNC_COMM_K<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_DENDRO_AMR_FAC :"<<ccz4::CCZ4_DENDRO_AMR_FAC<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_WAVELET_TOL :"<<ccz4::CCZ4_WAVELET_TOL<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_LOAD_IMB_TOL :"<<ccz4::CCZ4_LOAD_IMB_TOL<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_RK45_TIME_BEGIN :"<<ccz4::CCZ4_RK45_TIME_BEGIN<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_RK45_TIME_END :"<<ccz4::CCZ4_RK45_TIME_END<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_RK45_TIME_STEP_SIZE :"<<ccz4::CCZ4_RK45_TIME_STEP_SIZE<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_RK45_DESIRED_TOL :"<<ccz4::CCZ4_RK45_DESIRED_TOL<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_COMPD_MIN : ( :"<<ccz4::CCZ4_COMPD_MIN[0]<<" ,"<<ccz4::CCZ4_COMPD_MIN[1]<<","<<ccz4::CCZ4_COMPD_MIN[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_COMPD_MAX : ( :"<<ccz4::CCZ4_COMPD_MAX[0]<<" ,"<<ccz4::CCZ4_COMPD_MAX[1]<<","<<ccz4::CCZ4_COMPD_MAX[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_BLK_MIN : ( :"<<ccz4::CCZ4_BLK_MIN_X<<" ,"<<ccz4::CCZ4_BLK_MIN_Y<<","<<ccz4::CCZ4_BLK_MIN_Z<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_BLK_MAX : ( :"<<ccz4::CCZ4_BLK_MAX_X<<" ,"<<ccz4::CCZ4_BLK_MAX_Y<<","<<ccz4::CCZ4_BLK_MAX_Z<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_OCTREE_MIN : ( :"<<ccz4::CCZ4_OCTREE_MIN[0]<<" ,"<<ccz4::CCZ4_OCTREE_MIN[1]<<","<<ccz4::CCZ4_OCTREE_MIN[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_OCTREE_MAX : ( :"<<ccz4::CCZ4_OCTREE_MAX[0]<<" ,"<<ccz4::CCZ4_OCTREE_MAX[1]<<","<<ccz4::CCZ4_OCTREE_MAX[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tETA_CONST :"<<ccz4::ETA_CONST<<NRM<<std::endl;
        std::cout<<YLW<<"\tETA_R0 :"<<ccz4::ETA_R0<<NRM<<std::endl;
        std::cout<<YLW<<"\tETA_DAMPING :"<<ccz4::ETA_DAMPING<<NRM<<std::endl;
        std::cout<<YLW<<"\tETA_DAMPING_EXP :"<<ccz4::ETA_DAMPING_EXP<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_LAMBDA : ("<<ccz4::CCZ4_LAMBDA[0]<<" ,"<<ccz4::CCZ4_LAMBDA[1]<<","<<ccz4::CCZ4_LAMBDA[2]<<ccz4::CCZ4_LAMBDA[3]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_LAMBDA_F : ("<<ccz4::CCZ4_LAMBDA_F[0]<<" ,"<<ccz4::CCZ4_LAMBDA_F[1]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_KAPPA : ("<<ccz4::CCZ4_KAPPA[0]<<" ,"<<ccz4::CCZ4_KAPPA[1]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tPSI_FLOOR :"<<ccz4::PSI_FLOOR<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_TRK0 :"<<ccz4::CCZ4_TRK0<<NRM<<std::endl;
        std::cout<<YLW<<"\tKO_DISS_SIGMA :"<<ccz4::KO_DISS_SIGMA<<NRM<<std::endl;
        std::cout<<YLW<<"\tP_EXPO :"<<ccz4::P_EXPO<<NRM<<std::endl;
        std::cout<<YLW<<"\tANG_PAR :"<<ccz4::ANG_PAR<<NRM<<std::endl;

        std::cout<<YLW<<"\tBH1 MASS :"<<ccz4::BH1.getBHMass()<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH1 POSITION (x,y,z) : ("<<ccz4::BH1.getBHCoordX()<<", "<<ccz4::BH1.getBHCoordY()<<", "<<ccz4::BH1.getBHCoordZ()<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH1 VELOCITY (x,y,z) : ("<<ccz4::BH1.getVx()<<", "<<ccz4::BH1.getVy()<<", "<<ccz4::BH1.getVz()<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH1 SPIN (||,theta,phi): ( "<<ccz4::BH1.getBHSpin()<<", "<<ccz4::BH1.getBHSpinTheta()<<", "<<ccz4::BH1.getBHSpinPhi()<<" )"<<NRM<<std::endl;

        std::cout<<YLW<<"\tBH2 MASS :"<<ccz4::BH2.getBHMass()<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH2 POSITION (x,y,z) : ("<<ccz4::BH2.getBHCoordX()<<", "<<ccz4::BH2.getBHCoordY()<<", "<<ccz4::BH2.getBHCoordZ()<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH2 VELOCITY (x,y,z) : ("<<ccz4::BH2.getVx()<<", "<<ccz4::BH2.getVy()<<", "<<ccz4::BH2.getVz()<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH2 SPIN (||,theta,phi): ( "<<ccz4::BH2.getBHSpin()<<", "<<ccz4::BH2.getBHSpinTheta()<<", "<<ccz4::BH2.getBHSpinPhi()<<" )"<<NRM<<std::endl;

        std::cout<<YLW<<"\tCCZ4_DIM :"<<ccz4::CCZ4_DIM<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_MAXDEPTH :"<<ccz4::CCZ4_MAXDEPTH<<NRM<<std::endl;
	
	std::cout<<YLW<<"\tCCZ4_NUM_REFINE_VARS :"<<ccz4::CCZ4_NUM_REFINE_VARS<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_REFINE_VARIABLE_INDICES :[";
        for(unsigned int i=0;i<ccz4::CCZ4_NUM_REFINE_VARS-1;i++)
            std::cout<<ccz4::CCZ4_REFINE_VARIABLE_INDICES[i]<<", ";
        std::cout<<ccz4::CCZ4_REFINE_VARIABLE_INDICES[ccz4::CCZ4_NUM_REFINE_VARS-1]<<"]"<<NRM<<std::endl;

        std::cout<<YLW<<"\tCCZ4_NUM_EVOL_VARS_VTU_OUTPUT :"<<ccz4::CCZ4_NUM_EVOL_VARS_VTU_OUTPUT<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_VTU_OUTPUT_EVOL_INDICES :[";
        for(unsigned int i=0;i<ccz4::CCZ4_NUM_EVOL_VARS_VTU_OUTPUT-1;i++)
            std::cout<<ccz4::CCZ4_VTU_OUTPUT_EVOL_INDICES[i]<<", ";
        std::cout<<ccz4::CCZ4_VTU_OUTPUT_EVOL_INDICES[ccz4::CCZ4_NUM_EVOL_VARS_VTU_OUTPUT-1]<<"]"<<NRM<<std::endl;

        std::cout<<YLW<<"\tCCZ4_NUM_CONST_VARS_VTU_OUTPUT :"<<ccz4::CCZ4_NUM_CONST_VARS_VTU_OUTPUT<<NRM<<std::endl;
        std::cout<<YLW<<"\tCCZ4_VTU_OUTPUT_CONST_INDICES :[";
        for(unsigned int i=0;i<ccz4::CCZ4_NUM_CONST_VARS_VTU_OUTPUT-1;i++)
            std::cout<<ccz4::CCZ4_VTU_OUTPUT_CONST_INDICES[i]<<", ";
        std::cout<<ccz4::CCZ4_VTU_OUTPUT_CONST_INDICES[ccz4::CCZ4_NUM_CONST_VARS_VTU_OUTPUT-1]<<"]"<<NRM<<std::endl;

	std::cout<<YLW<<"\tTPID_TARGET_M_PLUS :"<<TPID::target_M_plus<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_TARGET_M_MINUS :"<<TPID::target_M_minus<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_PAR_B :"<<TPID::par_b<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_PAR_P_PLUS :"<<TPID::par_P_plus<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_PAR_P_MINUS :"<<TPID::par_P_minus<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_PAR_S_PLUS :"<<TPID::par_S_plus<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_PAR_S_MINUS :"<<TPID::par_S_minus<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_CENTER_OFFSET :"<<TPID::center_offset<<NRM<<std::endl;

        std::cout<<YLW<<"\tTPID_INITIAL_LAPSE_PSI_EXPONENT :"<<TPID::initial_lapse_psi_exponent<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_NPOINTS_A :"<<TPID::npoints_A<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_NPOINTS_B :"<<TPID::npoints_B<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_NPOINTS_PHI :"<<TPID::npoints_phi<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_GIVE_BARE_MASS :"<<TPID::give_bare_mass<<NRM<<std::endl;
        std::cout<<YLW<<"\tINITIAL_LAPSE :"<<TPID::initial_lapse<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_SOLVE_MOMENTUM_CONSTRAINT :"<<TPID::solve_momentum_constraint<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_GRID_SETUP_METHOD :"<<TPID::grid_setup_method<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_VERBOSE :"<<TPID::verbose<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_ADM_TOL :"<<TPID::adm_tol<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_NEWTON_TOL :"<<TPID::Newton_tol<<NRM<<std::endl;

    }

    _InitializeHcurve(ccz4::CCZ4_DIM);
    m_uiMaxDepth=ccz4::CCZ4_MAXDEPTH;

    if(ccz4::CCZ4_NUM_VARS%ccz4::CCZ4_ASYNC_COMM_K!=0)
    {
        if(!rank) std::cout<<"[overlap communication error]: total CCZ4_NUM_VARS: "<<ccz4::CCZ4_NUM_VARS<<" is not divisable by CCZ4_ASYNC_COMM_K: "<<ccz4::CCZ4_ASYNC_COMM_K<<std::endl;
        exit(0);
    }

    ccz4::CCZ4_RK45_TIME_STEP_SIZE=ccz4::CCZ4_CFL_FACTOR*(ccz4::CCZ4_COMPD_MAX[0]-ccz4::CCZ4_COMPD_MIN[0])*(1.0/(double)(1u<<ccz4::CCZ4_MAXDEPTH));

    //2. generate the initial grid.
    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){ccz4::punctureData(x,y,z,var);};

    const unsigned int interpVars=ccz4::CCZ4_NUM_VARS;
    unsigned int varIndex[interpVars];
    for(unsigned int i=0;i<ccz4::CCZ4_NUM_VARS;i++)
        varIndex[i]=i;

    /*varIndex[0]=ccz4::VAR::U_ALPHA;
    varIndex[1]=ccz4::VAR::U_CHI;*/
    DendroIntL localSz,globalSz;
    double t_stat;
    double t_stat_g[3];

    ccz4::timer::t_f2o.start();

    if(ccz4::CCZ4_ENABLE_BLOCK_ADAPTIVITY)
    {
        if(!rank) std::cout<<YLW<<"Using block adaptive mesh. AMR disabled "<<NRM<<std::endl;
        const Point pt_min(ccz4::CCZ4_BLK_MIN_X,ccz4::CCZ4_BLK_MIN_Y,ccz4::CCZ4_BLK_MIN_Z);
        const Point pt_max(ccz4::CCZ4_BLK_MAX_X,ccz4::CCZ4_BLK_MAX_Y,ccz4::CCZ4_BLK_MAX_Z);

        ccz4::blockAdaptiveOctree(tmpNodes,pt_min,pt_max,m_uiMaxDepth-2,m_uiMaxDepth,comm);
    }else
    {

        if(!rank) std::cout<<YLW<<"Using function2Octree. AMR enabled "<<NRM<<std::endl;
        function2Octree(f_init,ccz4::CCZ4_NUM_VARS,varIndex,interpVars,tmpNodes,m_uiMaxDepth,ccz4::CCZ4_WAVELET_TOL,ccz4::CCZ4_ELE_ORDER,comm);

    }

    ccz4::timer::t_f2o.stop();

    t_stat=ccz4::timer::t_f2o.seconds;
    par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,comm);
    t_stat_g[1]=t_stat_g[1]/(double)npes;

    localSz=tmpNodes.size();
    par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);

    if(!rank) std::cout<<GRN<<" function to octree max (s): "<<t_stat_g[2]<<NRM<<std::endl;
    if(!rank) std::cout<<GRN<<" function to octree # octants : "<<globalSz<<NRM<<std::endl;

    par::Mpi_Bcast(&globalSz,1,0,comm);
    const unsigned int grainSz=ccz4::CCZ4_DENDRO_GRAIN_SZ;//DENDRO_DEFAULT_GRAIN_SZ;

    bool isActive;
    MPI_Comm commActive;

    if((globalSz/grainSz)>=npes)
    {
        MPI_Comm_dup(comm,&commActive);
        isActive=true;

    }else
    {
        isActive=(rank*grainSz<globalSz);
        par::splitComm2way(isActive,&commActive,comm);

    }

    shrinkOrExpandOctree(tmpNodes,ccz4::CCZ4_LOAD_IMB_TOL,DENDRO_DEFAULT_SF_K,isActive,commActive,comm);

    if(!isActive)
        if(tmpNodes.size()!=0)
            std::cout<<" rank_g: "<<rank<<" isActive: "<<isActive<<" f2O octants: "<<tmpNodes.size()<<std::endl;




    std::vector<ot::TreeNode> balOct;
    localSz=0;
    if(isActive)
    {

        int rank_active,npes_active;

        MPI_Comm_size(commActive,&npes_active);
        MPI_Comm_rank(commActive,&rank_active);

        if(!rank_active) std::cout<<"[MPI_COMM_SWITCH]: "<<npes_active<<std::endl;

        ot::TreeNode root(ccz4::CCZ4_DIM,ccz4::CCZ4_MAXDEPTH);
        std::vector<ot::TreeNode> tmpVec;
        ccz4::timer::t_cons.start();

        SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,ccz4::CCZ4_LOAD_IMB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_REMOVE_DUPLICATES,ccz4::CCZ4_SPLIT_FIX,commActive);
        std::swap(tmpNodes,tmpVec);
        tmpVec.clear();

        SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,ccz4::CCZ4_LOAD_IMB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_CONSTRUCT_OCTREE,ccz4::CCZ4_SPLIT_FIX,commActive);
        std::swap(tmpNodes,tmpVec);
        tmpVec.clear();

        ccz4::timer::t_cons.stop();
        t_stat=ccz4::timer::t_cons.seconds;

        par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,commActive);
        par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,commActive);
        par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,commActive);
        t_stat_g[1]=t_stat_g[1]/(double)rank_active;

        localSz=tmpNodes.size();
        par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,commActive);

        if(!rank_active) std::cout<<GRN<<"remove duplicates + octree construction (s): "<<t_stat_g[2]<<NRM<<std::endl;
        if(!rank_active) std::cout<<GRN<<" # const. octants: "<<globalSz<<NRM<<std::endl;


        ccz4::timer::t_bal.start();

        SFC::parSort::SFC_treeSort(tmpNodes,balOct,balOct,balOct,ccz4::CCZ4_LOAD_IMB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_BALANCE_OCTREE,ccz4::CCZ4_SPLIT_FIX,commActive);
        tmpNodes.clear();

        ccz4::timer::t_bal.stop();

        t_stat=ccz4::timer::t_bal.seconds;
        par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,commActive);
        par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,commActive);
        par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,commActive);
        t_stat_g[1]=t_stat_g[1]/(double)rank_active;

        if(!rank_active) std::cout<<GRN<<" 2:1 balancing max (s): "<<t_stat_g[2]<<NRM<<std::endl;
        localSz=balOct.size();


    }
    MPI_Comm_free(&commActive);
    // all reduce act as barrier to sync all procs.
    par::Mpi_Allreduce(&localSz,&globalSz,1,MPI_SUM,comm);
    if(!rank) std::cout<<GRN<<" balanced # octants : "<<globalSz<<NRM<<std::endl;

    ccz4::timer::t_mesh.start();

    ot::Mesh * mesh=new ot::Mesh(balOct,1,ccz4::CCZ4_ELE_ORDER,comm,true,ccz4::CCZ4_DENDRO_GRAIN_SZ,ccz4::CCZ4_LOAD_IMB_TOL,ccz4::CCZ4_SPLIT_FIX);

    ccz4::timer::t_mesh.stop();

    t_stat=ccz4::timer::t_mesh.seconds;
    par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,comm);
    t_stat_g[1]=t_stat_g[1]/(double)npes;

    localSz=mesh->getNumLocalMeshNodes();
    par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
    if(!rank) std::cout<<GRN<<" # of CG nodes (vertices) : "<<globalSz<<NRM<<std::endl;
    if(!rank)
    {
        std::cout<< GRN<<"Mesh generation time (max): "<<t_stat_g[2]<<NRM<<std::endl;
        std::cout<<"\t"<<GRN<<" e2e (min,mean,max): "<<"( "<<t_e2e_g[0]<<"\t"<<t_e2e_g[1]<<"\t"<<t_e2e_g[2]<<" )"<<NRM<<std::endl;
        std::cout<<"\t"<<GRN<<" e2n (min,mean,max): "<<"( "<<t_e2n_g[0]<<"\t"<<t_e2n_g[1]<<"\t"<<t_e2n_g[2]<<" )"<<NRM<<std::endl;
        std::cout<<"\t"<<GRN<<" sm (min,mean,max): "<<"( "<<t_sm_g[0]<<"\t"<<t_sm_g[1]<<"\t"<<t_sm_g[2]<<" )"<<NRM<<std::endl;
        std::cout<<"\t"<<GRN<<" blk (min,mean,max): "<<"( "<<t_blk_g[0]<<"\t"<<t_blk_g[1]<<"\t"<<t_blk_g[2]<<" )"<<NRM<<std::endl;
    }


    //ode::solver::RK45_CCZ4 rk_ccz4(mesh,ccz4::CCZ4_RK45_TIME_BEGIN,ccz4::CCZ4_RK45_TIME_END,ccz4::CCZ4_RK45_TIME_STEP_SIZE);
    ode::solver::RK3_CCZ4 rk_ccz4(mesh,ccz4::CCZ4_RK45_TIME_BEGIN,ccz4::CCZ4_RK45_TIME_END,ccz4::CCZ4_RK45_TIME_STEP_SIZE);

    if(ccz4::CCZ4_RESTORE_SOLVER==1)
        rk_ccz4.restoreCheckPoint(ccz4::CCZ4_CHKPT_FILE_PREFIX.c_str(),comm);

    ccz4::timer::t_rkSolve.start();
    rk_ccz4.rkSolve();
    ccz4::timer::t_rkSolve.stop();

    ccz4::timer::total_runtime.stop();

    ccz4::timer::profileInfo(ccz4::CCZ4_PROFILE_FILE_PREFIX.c_str(),mesh);


    delete mesh;

    MPI_Finalize();

    return 0;
}
