//
// Created by milinda on 7/25/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief This file contains all the parameters related to CCZ4 simulation.
*/
//

#ifndef SFCSORTBENCH_PARAMETERS_H
#define SFCSORTBENCH_PARAMETERS_H

#include "bh.h"
#include <string.h>
#include <iostream>




namespace ccz4
{

    /**@brief element order*/
    static const unsigned int CCZ4_ELE_ORDER=4;

    /**@brief number of variables*/
    // HL : For ccz4, turn on this for 25 (for Theta_z4)
    static const unsigned int CCZ4_NUM_VARS=25;

    /**@brief number of constraints variables*/
    static const unsigned int CCZ4_CONSTRAINT_NUM_VARS=6;

    /***@brief number of RK45 stages*/
    static const unsigned int CCZ4_RK45_STAGES=6;

    /***@brief number of RK4 stages*/
    static const unsigned int CCZ4_RK4_STAGES=4;

    /**@brief number of rk4 stages*/
    static const unsigned int CCZ4_RK3_STAGES=3;

    /**@brief CFL stability number number (specifies how dt=CCZ4_CFL_FACTOR*dx)*/
    static const double CCZ4_CFL_FACTOR=0.1;

    /**@brief: parameter used for adaptive time step update. */
    static const double CCZ4_SAFETY_FAC=0.8;

    /**@brief number of internal variables*/
    static const unsigned int CCZ4_NUM_VARS_INTENL=(CCZ4_RK45_STAGES+1)*CCZ4_NUM_VARS;

    /**@brief min bh domain add these to the parameter file.*/
    extern double CCZ4_COMPD_MIN[3];
    /**@brief min bh domain @todo add these to the parameter file. */
    extern double CCZ4_COMPD_MAX[3];

    /**@brief min coords of the OCTREE */
    extern double CCZ4_OCTREE_MIN[3];
    /**@brief max coords of the OCTREE */
    extern double CCZ4_OCTREE_MAX[3];

    /**@brief solution output frequency*/
    extern unsigned int CCZ4_IO_OUTPUT_FREQ;

    /*brief timestep norms output frequency*/
    extern unsigned int CCZ4_TIME_STEP_OUTPUT_FREQ;

    /*brief remesh test frequency*/
    extern unsigned int CCZ4_REMESH_TEST_FREQ;

    /**@brief checkpoint store frequency*/
    extern unsigned int CCZ4_CHECKPT_FREQ;

    /**@brief restore the solver from check point if set to 1. */
    extern unsigned int CCZ4_RESTORE_SOLVER;

    /**@brief use the block adaptivity and disable the AMR*/
    extern unsigned int CCZ4_ENABLE_BLOCK_ADAPTIVITY;

    /**@brief file prefix for VTU*/
    extern std::string CCZ4_VTU_FILE_PREFIX;

    /**@brief file prefix for write check point*/
    extern std::string CCZ4_CHKPT_FILE_PREFIX;

    /**@brief file prefix to write profile info.*/
    extern std::string CCZ4_PROFILE_FILE_PREFIX;

    /**@brief number of refine variables */
    extern unsigned int CCZ4_NUM_REFINE_VARS;

    /**@brief indices of refine var ids*/
    extern unsigned int CCZ4_REFINE_VARIABLE_INDICES[CCZ4_NUM_VARS];

    /**@brief number of evolution variables written to vtu files*/
    extern unsigned int CCZ4_NUM_EVOL_VARS_VTU_OUTPUT;

    /**@brief number of constrint variables written to vtu files*/
    extern unsigned int CCZ4_NUM_CONST_VARS_VTU_OUTPUT;

    /**@brief evolution variable IDs written to vtu files*/
    extern unsigned int CCZ4_VTU_OUTPUT_EVOL_INDICES[CCZ4_NUM_VARS];

    /**@brief constraint variable IDs written to vtu files*/
    extern unsigned int CCZ4_VTU_OUTPUT_CONST_INDICES[CCZ4_CONSTRAINT_NUM_VARS];

    /**@brief solution output gap (instead of freq. we can use to output the solution if currentTime > lastIOOutputTime + CCZ4_IO_OUTPUT_GAP)*/
    extern  double CCZ4_IO_OUTPUT_GAP;

    /**@brief prefered grain sz to use when selecting active npes*/
    extern unsigned int CCZ4_DENDRO_GRAIN_SZ;

    /**@brief AMR coarsening factor (we coarsen if tol<CCZ4_DENDRO_AMR_FAC*CCZ4_WAVELET_TOL)*/
    extern double CCZ4_DENDRO_AMR_FAC;

    /**@brief wavelet tolerance value. */
    extern  double CCZ4_WAVELET_TOL;
    /**@brief load-imbalance tolerance value. */
    extern  double CCZ4_LOAD_IMB_TOL;
    /**@brief: Splitter fix value*/
    extern unsigned int CCZ4_SPLIT_FIX;

    /**@brief: async. communication at a time. (upper bound shoud be CCZ4_NUM_VARS) */
    extern unsigned int CCZ4_ASYNC_COMM_K;


    /**@brief simulation begin time. */
    extern double CCZ4_RK45_TIME_BEGIN;
    /**@brief simulation end time*/
    extern double CCZ4_RK45_TIME_END;
    /**@brief rk time step size. */
    extern double CCZ4_RK45_TIME_STEP_SIZE;

    /** desired tolerance value for the rk45 method (adaptive time stepping. )*/
    extern double CCZ4_RK45_DESIRED_TOL;

    /**@brief Black hole 1 */
    extern BH BH1;
    /**@brief Black hole 2 */
    extern BH BH2;

    /**@brief BBH initial data type */
    extern unsigned int CCZ4_ID_TYPE;

    /**@brief physical coordinates for grid, x_min */
    extern double CCZ4_GRID_MIN_X;

    /**@brief physical coordinates for grid, x_max */
    extern double CCZ4_GRID_MAX_X;

    /**@brief physical coordinates for grid, y_min */
    extern double CCZ4_GRID_MIN_Y;

    /**@brief physical coordinates for grid, y_max */
    extern double CCZ4_GRID_MAX_Y;

    /**@brief physical coordinates for grid, z_min */
    extern double CCZ4_GRID_MIN_Z;

    /**@brief physical coordinates for grid, z_max */
    extern double CCZ4_GRID_MAX_Z;

    /**@brief physical coordinates for the blk adaptive x_min*/
    extern double CCZ4_BLK_MIN_X;

    /**@brief physical coordinates for the blk adaptive x_min*/
    extern double CCZ4_BLK_MIN_Y;

    /**@brief physical coordinates for the blk adaptive x_min*/
    extern double CCZ4_BLK_MIN_Z;

    /**@brief physical coordinates for the blk adaptive x_min*/
    extern double CCZ4_BLK_MAX_X;

    /**@brief physical coordinates for the blk adaptive x_min*/
    extern double CCZ4_BLK_MAX_Y;

    /**@brief physical coordinates for the blk adaptive x_min*/
    extern double CCZ4_BLK_MAX_Z;

    /**@brief: dimension of the grid*/
    extern unsigned int CCZ4_DIM;

    /**@brief: max refinement level*/
    extern unsigned int CCZ4_MAXDEPTH;

    /**@brief: lambda values for evolution */
    extern unsigned int CCZ4_LAMBDA[4];

    /**@brief: lambda values for evolution */
    extern double CCZ4_LAMBDA_F[2];

    /**@brief: kappa values for ccz4 evolution */
    extern double CCZ4_KAPPA[2];

    /**@brief: lambda values for evolution */
    extern double CCZ4_TRK0;

    /**@brief: base value for eta in evolution */
    extern double ETA_CONST;

    /**@brief: eta_R0, radius where eta is damped for evolution */
    extern double ETA_R0;

    /**@brief: eta damping for evolution */
    extern double ETA_DAMPING;

    /**@brief: eta damping exponent for evolution */
    extern double ETA_DAMPING_EXP;

    /**@brief: psi floor value */
    extern double PSI_FLOOR;

    /**@brief: Kreiss-Oliger dissipation */
    extern double KO_DISS_SIGMA;

    /**@brief: Exponential for conformal factor */
    extern double P_EXPO;

    /**@brief: Angular parameter for Kerr-Schild */
    extern double ANG_PAR;
}

namespace TPID {
  static const double TP_epsilon = 1.0e-6;
  static const int swap_xz = 0;
  static const int use_sources = 0;
  static const int rescale_sources = 0;
  static const int use_external_initial_guess = 0;
  static const int do_residuum_debug_output = 1;
  static const int do_initial_debug_output = 1;
  static const int multiply_old_lapse = 0;
  static const double TP_Tiny = 1.0e-15;
  static const double TP_Extend_Radius = 0.0;
  static const int Newton_maxit = 5;

  extern double target_M_plus;
  extern double target_M_minus;
  extern double par_m_plus;
  extern double par_m_minus;
  extern double par_b;
  extern double par_P_plus[3];
  extern double par_P_minus[3];
  extern double par_S_plus[3];
  extern double par_S_minus[3];
  extern double center_offset[3];
  extern double initial_lapse_psi_exponent;
  extern int npoints_A;
  extern int npoints_B;
  extern int npoints_phi;
  extern int give_bare_mass;
  extern int initial_lapse;
  extern int solve_momentum_constraint;
  extern int grid_setup_method;
  extern int verbose;
  extern double adm_tol;
  extern double Newton_tol;
}

#endif //SFCSORTBENCH_PARAMETERS_H
