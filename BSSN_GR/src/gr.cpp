//
// Created by milinda on 7/25/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Header file for the GR simulation.
*/
//

#include "gr.h"
#include "grUtils.h"
#include "mpi.h"
#include "TreeNode.h"
#include "mesh.h"
#include <vector>
#include <iostream>
#include "rkBSSN.h"
#include "octUtils.h"
#include "meshUtils.h"

int main (int argc, char** argv)
{


    if(argc<2)
        std::cout<<"Usage: "<<argv[0]<<" paramFile"<<std::endl;

    MPI_Init(&argc,&argv);
    MPI_Comm comm=MPI_COMM_WORLD;

    int rank,npes;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);

    bssn::timer::initFlops();

    bssn::timer::total_runtime.start();

    // Print out CMAKE options
    if (!rank) {
        #ifdef BSSN_COMPUTE_CONSTRAINTS
          std::cout<<GRN<<"  Compiled with BSSN_COMPUTE_CONSTRAINTS"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without BSSN_COMPUTE_CONSTRAINTS"<<NRM<<std::endl;
        #endif
        #ifdef BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT
          std::cout<<GRN<<"  Compiled with BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT"<<NRM<<std::endl;
        #endif
        #ifdef BSSN_ENABLE_VTU_OUTPUT
          std::cout<<GRN<<"  Compiled with BSSN_ENABLE_VTU_OUTPUT"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without BSSN_ENABLE_VTU_OUTPUT"<<NRM<<std::endl;
        #endif
        #ifdef BSSN_ETA_FUNCTION 
          std::cout<<GRN<<"  Compiled with  BSSN_ETA_FUNCTION"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  BSSN_ETA_FUNCTION"<<NRM<<std::endl;
        #endif
        #ifdef BSSN_EXTRACT_BH_LOCATIONS 
          std::cout<<GRN<<"  Compiled with  BSSN_EXTRACT_BH_LOCATIONS"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  BSSN_EXTRACT_BH_LOCATIONS"<<NRM<<std::endl;
        #endif
        #ifdef BSSN_EXTRACT_GRAVITATIONAL_WAVES 
          std::cout<<GRN<<"  Compiled with  BSSN_EXTRACT_GRAVITATIONAL_WAVES"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  BSSN_EXTRACT_GRAVITATIONAL_WAVES"<<NRM<<std::endl;
        #endif
        #ifdef BSSN_EXTRACT_GRAVITATIONAL_WAVES 
          std::cout<<GRN<<"  Compiled with  BSSN_EXTRACT_GRAVITATIONAL_WAVES"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  BSSN_EXTRACT_GRAVITATIONAL_WAVES"<<NRM<<std::endl;
        #endif
        #ifdef BSSN_GAUGE_ROCHESTER 
          std::cout<<GRN<<"  Compiled with  BSSN_GAUGE_ROCHESTER"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  BSSN_GAUGE_ROCHESTER"<<NRM<<std::endl;
        #endif
        #ifdef BSSN_KERR_SCHILD_TEST 
          std::cout<<GRN<<"  Compiled with  BSSN_KERR_SCHILD_TEST"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  BSSN_KERR_SCHILD_TEST"<<NRM<<std::endl;
        #endif

        #ifdef BSSN_REFINE_BASE_EH 
          std::cout<<GRN<<"  Compiled with  BSSN_REFINE_BASE_EH"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  BSSN_REFINE_BASE_EH"<<NRM<<std::endl;
        #endif

        #ifdef USE_FD_INTERP_FOR_UNZIP 
          std::cout<<GRN<<"  Compiled with  USE_FD_INTERP_FOR_UNZIP"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  USE_FD_INTERP_FOR_UNZIP"<<NRM<<std::endl;
        #endif

    }


    //1 . read the parameter file.
    if(!rank) std::cout<<" reading parameter file :"<<argv[1]<<std::endl;
    bssn::readParamFile(argv[1],comm);

    if(rank==1|| npes==1)
    {
        std::cout<<"parameters read: "<<std::endl;

        std::cout<<YLW<<"\tnpes :"<<npes<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_DIM :"<<bssn::BSSN_DIM<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_ELE_ORDER :"<<bssn::BSSN_ELE_ORDER<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_CFL_FACTOR :"<<bssn::BSSN_CFL_FACTOR<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_IO_OUTPUT_FREQ :"<<bssn::BSSN_IO_OUTPUT_FREQ<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_GW_EXTRACT_FREQ :"<<bssn::BSSN_GW_EXTRACT_FREQ<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_REMESH_TEST_FREQ :"<<bssn::BSSN_REMESH_TEST_FREQ<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_CHECKPT_FREQ :"<<bssn::BSSN_CHECKPT_FREQ<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_RESTORE_SOLVER :"<<bssn::BSSN_RESTORE_SOLVER<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_ENABLE_BLOCK_ADAPTIVITY :"<<bssn::BSSN_ENABLE_BLOCK_ADAPTIVITY<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_VTU_FILE_PREFIX :"<<bssn::BSSN_VTU_FILE_PREFIX<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_CHKPT_FILE_PREFIX :"<<bssn::BSSN_CHKPT_FILE_PREFIX<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_PROFILE_FILE_PREFIX :"<<bssn::BSSN_PROFILE_FILE_PREFIX<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_VTU_Z_SLICE_ONLY :"<<bssn::BSSN_VTU_Z_SLICE_ONLY<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_IO_OUTPUT_GAP :"<<bssn::BSSN_IO_OUTPUT_GAP<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_DENDRO_GRAIN_SZ :"<<bssn::BSSN_DENDRO_GRAIN_SZ<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_ASYNC_COMM_K :"<<bssn::BSSN_ASYNC_COMM_K<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_DENDRO_AMR_FAC :"<<bssn::BSSN_DENDRO_AMR_FAC<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_USE_WAVELET_TOL_FUNCTION :"<<bssn::BSSN_USE_WAVELET_TOL_FUNCTION<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_USE_FD_GRID_TRANSFER :"<<bssn::BSSN_USE_FD_GRID_TRANSFER<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_WAVELET_TOL :"<<bssn::BSSN_WAVELET_TOL<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_WAVELET_TOL_MAX:"<<bssn::BSSN_WAVELET_TOL_MAX<<NRM<<std::endl;
        std::cout<<YLW<<"\t:BSSN_WAVELET_TOL_FUNCTION_R0: "<<bssn::BSSN_WAVELET_TOL_FUNCTION_R0<<NRM<<std::endl;
        std::cout<<YLW<<"\t:BSSN_WAVELET_TOL_FUNCTION_R1: "<<bssn::BSSN_WAVELET_TOL_FUNCTION_R1<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_LOAD_IMB_TOL :"<<bssn::BSSN_LOAD_IMB_TOL<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_RK_TIME_BEGIN :"<<bssn::BSSN_RK_TIME_BEGIN<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_RK_TIME_END :"<<bssn::BSSN_RK_TIME_END<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_RK_TYPE :"<<bssn::BSSN_RK_TYPE<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_RK45_TIME_STEP_SIZE :"<<bssn::BSSN_RK45_TIME_STEP_SIZE<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_RK45_DESIRED_TOL :"<<bssn::BSSN_RK45_DESIRED_TOL<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_COMPD_MIN : ( :"<<bssn::BSSN_COMPD_MIN[0]<<" ,"<<bssn::BSSN_COMPD_MIN[1]<<","<<bssn::BSSN_COMPD_MIN[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_COMPD_MAX : ( :"<<bssn::BSSN_COMPD_MAX[0]<<" ,"<<bssn::BSSN_COMPD_MAX[1]<<","<<bssn::BSSN_COMPD_MAX[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_BLK_MIN : ( :"<<bssn::BSSN_BLK_MIN_X<<" ,"<<bssn::BSSN_BLK_MIN_Y<<","<<bssn::BSSN_BLK_MIN_Z<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_BLK_MAX : ( :"<<bssn::BSSN_BLK_MAX_X<<" ,"<<bssn::BSSN_BLK_MAX_Y<<","<<bssn::BSSN_BLK_MAX_Z<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_OCTREE_MIN : ( :"<<bssn::BSSN_OCTREE_MIN[0]<<" ,"<<bssn::BSSN_OCTREE_MIN[1]<<","<<bssn::BSSN_OCTREE_MIN[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_OCTREE_MAX : ( :"<<bssn::BSSN_OCTREE_MAX[0]<<" ,"<<bssn::BSSN_OCTREE_MAX[1]<<","<<bssn::BSSN_OCTREE_MAX[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tETA_CONST :"<<bssn::ETA_CONST<<NRM<<std::endl;
        std::cout<<YLW<<"\tETA_R0 :"<<bssn::ETA_R0<<NRM<<std::endl;
        std::cout<<YLW<<"\tETA_DAMPING :"<<bssn::ETA_DAMPING<<NRM<<std::endl;
        std::cout<<YLW<<"\tETA_DAMPING_EXP :"<<bssn::ETA_DAMPING_EXP<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_ETA_R0 :"<<bssn::BSSN_ETA_R0<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_ETA_POWER : ("<<bssn::BSSN_ETA_POWER[0]<<" ,"<<bssn::BSSN_ETA_POWER[1]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_LAMBDA : ("<<bssn::BSSN_LAMBDA[0]<<" ,"<<bssn::BSSN_LAMBDA[1]<<","<<bssn::BSSN_LAMBDA[2]<<bssn::BSSN_LAMBDA[3]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_LAMBDA_F : ("<<bssn::BSSN_LAMBDA_F[0]<<" ,"<<bssn::BSSN_LAMBDA_F[1]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_XI : ("<<bssn::BSSN_XI[0]<<" ,"<<bssn::BSSN_XI[1]<<" ,"<<bssn::BSSN_XI[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tCHI_FLOOR :"<<bssn::CHI_FLOOR<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_TRK0 :"<<bssn::BSSN_TRK0<<NRM<<std::endl;
        std::cout<<YLW<<"\tDISSIPATION_TYPE :"<<bssn::DISSIPATION_TYPE<<NRM<<std::endl;
        std::cout<<YLW<<"\tKO_DISS_SIGMA :"<<bssn::KO_DISS_SIGMA<<NRM<<std::endl;

        std::cout<<YLW<<"\tBH1 MASS :"<<bssn::BH1.getBHMass()<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH1 POSITION (x,y,z) : ("<<bssn::BH1.getBHCoordX()<<", "<<bssn::BH1.getBHCoordY()<<", "<<bssn::BH1.getBHCoordZ()<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH1 VELOCITY (x,y,z) : ("<<bssn::BH1.getVx()<<", "<<bssn::BH1.getVy()<<", "<<bssn::BH1.getVz()<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH1 SPIN (||,theta,phi): ( "<<bssn::BH1.getBHSpin()<<", "<<bssn::BH1.getBHSpinTheta()<<", "<<bssn::BH1.getBHSpinPhi()<<" )"<<NRM<<std::endl;

        std::cout<<YLW<<"\tBH2 MASS :"<<bssn::BH2.getBHMass()<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH2 POSITION (x,y,z) : ("<<bssn::BH2.getBHCoordX()<<", "<<bssn::BH2.getBHCoordY()<<", "<<bssn::BH2.getBHCoordZ()<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH2 VELOCITY (x,y,z) : ("<<bssn::BH2.getVx()<<", "<<bssn::BH2.getVy()<<", "<<bssn::BH2.getVz()<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH2 SPIN (||,theta,phi): ( "<<bssn::BH2.getBHSpin()<<", "<<bssn::BH2.getBHSpinTheta()<<", "<<bssn::BH2.getBHSpinPhi()<<" )"<<NRM<<std::endl;

        std::cout<<YLW<<"\tBSSN_DIM :"<<bssn::BSSN_DIM<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_MAXDEPTH :"<<bssn::BSSN_MAXDEPTH<<NRM<<std::endl;

        std::cout<<YLW<<"\tBSSN_NUM_REFINE_VARS :"<<bssn::BSSN_NUM_REFINE_VARS<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_REFINE_VARIABLE_INDICES :[";
        for(unsigned int i=0;i<bssn::BSSN_NUM_REFINE_VARS-1;i++)
            std::cout<<bssn::BSSN_REFINE_VARIABLE_INDICES[i]<<", ";
        std::cout<<bssn::BSSN_REFINE_VARIABLE_INDICES[bssn::BSSN_NUM_REFINE_VARS-1]<<"]"<<NRM<<std::endl;

        std::cout<<YLW<<"\tBSSN_REFINEMENT_MODE :"<<bssn::BSSN_REFINEMENT_MODE<<NRM<<std::endl;

        #ifdef BSSN_REFINE_BASE_EH
                std::cout<<YLW<<"\tBSSN_EH_REFINE_VAL  : "<<bssn::BSSN_EH_REFINE_VAL<<NRM<<std::endl;
                std::cout<<YLW<<"\tBSSN_EH_COARSEN_VAL : "<<bssn::BSSN_EH_COARSEN_VAL<<NRM<<std::endl;
        #endif 

        std::cout<<YLW<<"\tBSSN_NUM_EVOL_VARS_VTU_OUTPUT :"<<bssn::BSSN_NUM_EVOL_VARS_VTU_OUTPUT<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_VTU_OUTPUT_EVOL_INDICES :[";
        for(unsigned int i=0;i<bssn::BSSN_NUM_EVOL_VARS_VTU_OUTPUT-1;i++)
            std::cout<<bssn::BSSN_VTU_OUTPUT_EVOL_INDICES[i]<<", ";
        std::cout<<bssn::BSSN_VTU_OUTPUT_EVOL_INDICES[bssn::BSSN_NUM_EVOL_VARS_VTU_OUTPUT-1]<<"]"<<NRM<<std::endl;

        std::cout<<YLW<<"\tBSSN_NUM_CONST_VARS_VTU_OUTPUT :"<<bssn::BSSN_NUM_CONST_VARS_VTU_OUTPUT<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_VTU_OUTPUT_CONST_INDICES :[";
        for(unsigned int i=0;i<bssn::BSSN_NUM_CONST_VARS_VTU_OUTPUT-1;i++)
            std::cout<<bssn::BSSN_VTU_OUTPUT_CONST_INDICES[i]<<", ";
        std::cout<<bssn::BSSN_VTU_OUTPUT_CONST_INDICES[bssn::BSSN_NUM_CONST_VARS_VTU_OUTPUT-1]<<"]"<<NRM<<std::endl;


        std::cout<<YLW<<"\tTPID_TARGET_M_PLUS :"<<TPID::target_M_plus<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_TARGET_M_MINUS :"<<TPID::target_M_minus<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_PAR_B :"<<TPID::par_b<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_PAR_P_PLUS : ( "<<TPID::par_P_plus[0]<<", "<<TPID::par_P_plus[1]<<", "<<TPID::par_P_plus[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_PAR_P_MINUS : ( "<<TPID::par_P_minus[0]<<", "<<TPID::par_P_minus[1]<<", "<<TPID::par_P_minus[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_PAR_S_PLUS : ( "<<TPID::par_S_plus[0]<<", "<<TPID::par_S_plus[1]<<", "<<TPID::par_S_plus[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_PAR_S_MINUS : ( "<<TPID::par_S_minus[0]<<", "<<TPID::par_S_minus[1]<<", "<<TPID::par_S_minus[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_CENTER_OFFSET : ( "<<TPID::center_offset[0]<<", "<<TPID::center_offset[1]<<", "<<TPID::center_offset[2]<<" )"<<NRM<<std::endl;

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
        
        
        std::cout<<YLW<<"\tEXTRACTION_VAR_ID :"<<BHLOC::EXTRACTION_VAR_ID<<NRM<<std::endl;
        std::cout<<YLW<<"\tEXTRACTION_TOL :"<<BHLOC::EXTRACTION_TOL<<NRM<<std::endl;


        std::cout<<YLW<<"\tBSSN_GW_NUM_RADAII: "<<GW::BSSN_GW_NUM_RADAII<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_GW_NUM_LMODES: "<<GW::BSSN_GW_NUM_LMODES<<NRM<<std::endl;

        std::cout<<YLW<<"\tBSSN_GW_RADAII: {";
        for(unsigned int i=0;i<GW::BSSN_GW_NUM_RADAII;i++)
            std::cout<<" ,"<<GW::BSSN_GW_RADAII[i];
        std::cout<<"}"<<NRM<<std::endl;

        std::cout<<YLW<<"\tBSSN_GW_L_MODES: {";
        for(unsigned int i=0;i<GW::BSSN_GW_NUM_LMODES;i++)
            std::cout<<" ,"<<GW::BSSN_GW_L_MODES[i];
        std::cout<<"}"<<NRM<<std::endl;

        

    }

    _InitializeHcurve(bssn::BSSN_DIM);
    m_uiMaxDepth=bssn::BSSN_MAXDEPTH;
    
    if(bssn::BSSN_NUM_VARS%bssn::BSSN_ASYNC_COMM_K!=0)
    {
        if(!rank) std::cout<<"[overlap communication error]: total BSSN_NUM_VARS: "<<bssn::BSSN_NUM_VARS<<" is not divisable by BSSN_ASYNC_COMM_K: "<<bssn::BSSN_ASYNC_COMM_K<<std::endl;
        MPI_Abort(comm,0);
    }

    if(bssn::BSSN_GW_EXTRACT_FREQ> bssn::BSSN_IO_OUTPUT_FREQ)
    {
      if(!rank) std::cout<<" BSSN_GW_EXTRACT_FREQ  should be less BSSN_IO_OUTPUT_FREQ "<<std::endl;
      MPI_Abort(comm,0);
    }


    //2. generate the initial grid.
    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){bssn::punctureData(x,y,z,var);};
    std::function<double(double,double,double)> f_init_alpha=[](double x,double y,double z){ double var[24]; bssn::punctureData(x,y,z,var); return var[0];};
    //std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){bssn::KerrSchildData(x,y,z,var);};

    const unsigned int interpVars=bssn::BSSN_NUM_VARS;
    unsigned int varIndex[interpVars];
    for(unsigned int i=0;i<bssn::BSSN_NUM_VARS;i++)
        varIndex[i]=i;

    /*varIndex[0]=bssn::VAR::U_ALPHA;
    varIndex[1]=bssn::VAR::U_CHI;*/
    DendroIntL localSz,globalSz;
    double t_stat;
    double t_stat_g[3];

    bssn::timer::t_f2o.start();

    if(bssn::BSSN_ENABLE_BLOCK_ADAPTIVITY)
    {
        if(!rank) std::cout<<YLW<<"Using block adaptive mesh. AMR disabled "<<NRM<<std::endl;
        const Point pt_min(bssn::BSSN_BLK_MIN_X,bssn::BSSN_BLK_MIN_Y,bssn::BSSN_BLK_MIN_Z);
        const Point pt_max(bssn::BSSN_BLK_MAX_X,bssn::BSSN_BLK_MAX_Y,bssn::BSSN_BLK_MAX_Z);

        bssn::blockAdaptiveOctree(tmpNodes,pt_min,pt_max,m_uiMaxDepth-2,m_uiMaxDepth,comm);
    }else
    {

        if(!rank) std::cout<<YLW<<"Using function2Octree. AMR enabled "<<NRM<<std::endl;
        function2Octree(f_init,bssn::BSSN_NUM_VARS,varIndex,interpVars,tmpNodes,m_uiMaxDepth,bssn::BSSN_WAVELET_TOL,bssn::BSSN_ELE_ORDER,comm);
    }

    ot::Mesh * mesh= ot::createMesh(tmpNodes.data(),tmpNodes.size(),bssn::BSSN_ELE_ORDER,comm,1,ot::SM_TYPE::FDM,bssn::BSSN_DENDRO_GRAIN_SZ,bssn::BSSN_LOAD_IMB_TOL,bssn::BSSN_SPLIT_FIX, bssn::getOctantWeight);
    unsigned int lmin, lmax;
    mesh->computeMinMaxLevel(lmin,lmax);    
    bssn::BSSN_RK45_TIME_STEP_SIZE=bssn::BSSN_CFL_FACTOR*((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))));
    par::Mpi_Bcast(&bssn::BSSN_RK45_TIME_STEP_SIZE,1,0,comm);


    ode::solver::RK_BSSN rk_bssn(mesh,bssn::BSSN_RK_TIME_BEGIN,bssn::BSSN_RK_TIME_END,bssn::BSSN_RK45_TIME_STEP_SIZE,(RKType)bssn::BSSN_RK_TYPE);
    
    if(bssn::BSSN_RESTORE_SOLVER==1)
        rk_bssn.restoreCheckPoint(bssn::BSSN_CHKPT_FILE_PREFIX.c_str(),comm);

    bssn::timer::t_rkSolve.start();
    rk_bssn.rkSolve();
    bssn::timer::t_rkSolve.stop();

    bssn::timer::total_runtime.stop();
    rk_bssn.freeMesh();

    
    MPI_Finalize();
    return 0;

}
