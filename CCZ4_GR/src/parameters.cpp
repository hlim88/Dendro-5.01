//
// Created by milinda on 8/23/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief
*/
//

#include "parameters.h"

namespace ccz4
{
    unsigned int CCZ4_IO_OUTPUT_FREQ=10;
    unsigned int CCZ4_TIME_STEP_OUTPUT_FREQ=10;
    unsigned int CCZ4_REMESH_TEST_FREQ=10;
    double CCZ4_IO_OUTPUT_GAP=1.0;
    double CCZ4_WAVELET_TOL=0.0001;
    double CCZ4_LOAD_IMB_TOL=0.1;
    unsigned int CCZ4_SPLIT_FIX=2;
    unsigned int CCZ4_ASYNC_COMM_K=4;
    double CCZ4_RK45_TIME_BEGIN=0;
    double CCZ4_RK45_TIME_END=10;

    double CCZ4_RK45_DESIRED_TOL=1e-6;


    unsigned int CCZ4_CHECKPT_FREQ=10;
    unsigned int CCZ4_RESTORE_SOLVER=0;
    unsigned int CCZ4_ENABLE_BLOCK_ADAPTIVITY=0;

    std::string CCZ4_VTU_FILE_PREFIX="ccz4_gr";
    std::string CCZ4_CHKPT_FILE_PREFIX="ccz4_cp";
    std::string CCZ4_PROFILE_FILE_PREFIX="ccz4_prof";


    BH BH1;
    BH BH2;
    unsigned int CCZ4_DIM=3;
    unsigned int CCZ4_MAXDEPTH=8;

    unsigned int CCZ4_ID_TYPE=0;


    double CCZ4_GRID_MIN_X=-50.0;
    double CCZ4_GRID_MAX_X=50.0;
    double CCZ4_GRID_MIN_Y=-50.0;
    double CCZ4_GRID_MAX_Y=50.0;
    double CCZ4_GRID_MIN_Z=-50.0;
    double CCZ4_GRID_MAX_Z=50.0;

    double CCZ4_BLK_MIN_X=-6.0;
    double CCZ4_BLK_MIN_Y=-6.0;
    double CCZ4_BLK_MIN_Z=-6.0;

    double CCZ4_BLK_MAX_X=6.0;
    double CCZ4_BLK_MAX_Y=6.0;
    double CCZ4_BLK_MAX_Z=6.0;

    double CCZ4_COMPD_MIN[3]={CCZ4_GRID_MIN_X,CCZ4_GRID_MIN_Y,CCZ4_GRID_MIN_Z};
    double CCZ4_COMPD_MAX[3]={CCZ4_GRID_MAX_X,CCZ4_GRID_MAX_Y,CCZ4_GRID_MAX_Z};

    double CCZ4_OCTREE_MIN[3]={0.0,0.0,0.0};
    double CCZ4_OCTREE_MAX[3]={(double)(1u<<CCZ4_MAXDEPTH),(double)(1u<<CCZ4_MAXDEPTH),(double)(1u<<CCZ4_MAXDEPTH)};

    //@note assumes the computational domain is a cube as well.
    double CCZ4_RK45_TIME_STEP_SIZE=CCZ4_CFL_FACTOR*(CCZ4_COMPD_MAX[0]-CCZ4_COMPD_MIN[0])*(1.0/(double)(1u<<CCZ4_MAXDEPTH));

    unsigned int CCZ4_LAMBDA[4]={1, 1, 1, 1};
    double CCZ4_LAMBDA_F[2]={1.0, 0.0};
    double CCZ4_KAPPA[2]={1.0, 1.0};
    double CCZ4_TRK0=0.0;
    double ETA_CONST=2.0;
    double ETA_R0=50.0;
    double ETA_DAMPING=1.0;
    double ETA_DAMPING_EXP=1.0;
    double PSI_FLOOR=0.1;
    double KO_DISS_SIGMA=0.01;
    double P_EXPO=-1.0;
    double ANG_PAR=0.01;


    unsigned int CCZ4_DENDRO_GRAIN_SZ=1000;

    double CCZ4_DENDRO_AMR_FAC=0.1;

    unsigned int CCZ4_NUM_REFINE_VARS=1;
    unsigned int CCZ4_REFINE_VARIABLE_INDICES[CCZ4_NUM_VARS]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};

    unsigned int CCZ4_NUM_EVOL_VARS_VTU_OUTPUT=1;
    unsigned int CCZ4_NUM_CONST_VARS_VTU_OUTPUT=1;
    unsigned int CCZ4_VTU_OUTPUT_EVOL_INDICES[CCZ4_NUM_VARS]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
    unsigned int CCZ4_VTU_OUTPUT_CONST_INDICES[CCZ4_CONSTRAINT_NUM_VARS]={0,1,2,3,4,5};

}

namespace TPID {
  double target_M_plus=1.0;
  double target_M_minus=1.0;
  double par_m_plus=1.0;
  double par_m_minus=1.0;
  double par_b=4.0;
  double par_P_plus[3]={0.0,0.0,0.0};
  double par_P_minus[3]={0.0,0.0,0.0};
  double par_S_plus[3]={0.0,0.0,0.0};
  double par_S_minus[3]={0.0,0.0,0.0};
  double center_offset[3]={0.0,0.0,0.00014142135623730951};
  double initial_lapse_psi_exponent=-2;
  int npoints_A=30;
  int npoints_B=30;
  int npoints_phi=16;
  int give_bare_mass=0;
  int initial_lapse=2;
  int solve_momentum_constraint=1;
  int grid_setup_method=1;
  int verbose = 1;
  double adm_tol=1.0e-10;
  double Newton_tol=1.0e-10;
}
