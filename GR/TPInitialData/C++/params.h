#ifndef _HAVE_PARAMS_H
#define _HAVE_PARAMS_H

#include "cctk.h"

namespace TPID {

#if ID_PARS == 0
  // data for GW150914
  static const int nx0 = 65;
  static const int ny0 = 65;
  static const int nz0 = 65;

  static const double xmin = -20.0;
  static const double xmax =  20.0;
  static const double ymin = -20.0;
  static const double ymax =  20.0;
  static const double zmin = -20.0;
  static const double zmax =  20.0;

  static const double target_M_plus =  0.5538461538461539;
  static const double target_M_minus =  0.4461538461538462;

  extern double par_m_plus;
  extern double par_m_minus;

  static const double par_b =  5.0;
  static const double center_offset[3] =  {-0.5384615384615383, 0, 0};
  static const double par_P_plus[3] =  {-0.00084541526517121, 0.09530152296974252, 0.0};
  static const double par_P_minus[3] =  {0.00084541526517121, -0.09530152296974252, -0.0};
  static const double par_S_plus[3] =  {0.0, 0.0, 0.09509112426035504};
  static const double par_S_minus[3] =  {0.0, 0.0, -0.09156449704142013};

  static const int initial_lapse = LAPSE_AVERAGED;
  static const int grid_setup_method = EVALUATION;
  static const int give_bare_mass = 0;
  static const double TP_epsilon = 1.0e-6;

  static const int swap_xz = 0;
  static const int verbose = 1;

  static const int use_sources = 0;
  static const int rescale_sources = 0;
  static const int use_external_initial_guess = 0;
  static const int do_residuum_debug_output = 1;
  static const int do_initial_debug_output = 1;
  static const int multiply_old_lapse = 0;
  static const int solve_momentum_constraint = 1;

  static const double adm_tol = 1.0e-10;
//  static const double adm_tol = 1.0e-6;

  static const double Newton_tol = 1.0e-10;
  static const double TP_Tiny = 1.0e-15;
  static const double TP_Extend_Radius = 0.0;

  static const double initial_lapse_psi_exponent = -2.0;

  static const int Newton_maxit = 5;

  static const int npoints_A = 30;
  static const int npoints_B = 30;
  static const int npoints_phi = 16;

#elif ID_PARS == 1

  static const int nx0 = 35;
  static const int ny0 = 35;
  static const int nz0 = 35;

  static const double xmin = -4.0;
  static const double xmax =  4.0;
  static const double ymin = -4.0;
  static const double ymax =  4.0;
  static const double zmin = -4.0;
  static const double zmax =  4.0;

  static const int initial_lapse = LAPSE_AVERAGED;
  static const int give_bare_mass = 1;
  static const int grid_setup_method = 1;

  static const int swap_xz = 0;

  static const int verbose = 1;

  static const int use_sources = 0;
  static const int rescale_sources = 0;
  static const int use_external_initial_guess = 0;
  static const int do_residuum_debug_output = 1;
  static const int do_initial_debug_output = 1;
  static const int multiply_old_lapse = 0;
  static const int solve_momentum_constraint = 0;

  static const double adm_tol = 1.0e-10;
  static const double Newton_tol = 1.0e-10;
  static const double TP_epsilon = 0.0;
  static const double TP_Tiny = 0.0;
  static const double TP_Extend_Radius = 0.0;

  static const double par_b = 1.5;

  extern double par_m_plus;
  extern double par_m_minus;

  static const double target_M_plus = 1.5;
  static const double target_M_minus = 1.0;

  static const double par_P_plus[3] = { 0.0, 2.0, 0.0 };
  static const double par_P_minus[3] = { 0.0, -2.0, 0.0 };
  static const double par_S_plus[3] = { 0.0, 0.5, -0.5 };
  static const double par_S_minus[3] = { -1.0, 0.0, -1.0 };
  static const double center_offset[3] = { 0.0, 0.0, 0.0 };

  static const double initial_lapse_psi_exponent = -2.0;

  static const int Newton_maxit = 5;

  static const int npoints_A = 20;
  static const int npoints_B = 20;
  static const int npoints_phi = 12;

#elif ID_PARS == 2

 /* Data for comparison to approximate solution in HAD. */

  static const int nx0 = 65;
  static const int ny0 = 65;
  static const int nz0 = 65;

  static const double xmin = -10.0;
  static const double xmax =  10.0;
  static const double ymin = -10.0;
  static const double ymax =  10.0;
  static const double zmin = -10.0;
  static const double zmax =  10.0;


  static const int initial_lapse = LAPSE_PSIN;
  static const int give_bare_mass = 1;
  static const int grid_setup_method = 1;

  static const int swap_xz = 0;

  static const int verbose = 1;

  static const int use_sources = 0;
  static const int rescale_sources = 0;
  static const int use_external_initial_guess = 0;
  static const int do_residuum_debug_output = 1;
  static const int do_initial_debug_output = 1;
  static const int multiply_old_lapse = 0;
  static const int solve_momentum_constraint = 0;

  static const double adm_tol = 1.0e-10;
  static const double Newton_tol = 1.0e-10;
  static const double TP_epsilon = 0.0;
  static const double TP_Tiny = 0.0;
  static const double TP_Extend_Radius = 0.0;

  static const double par_b = 4.0;

  extern double par_m_plus;
  extern double par_m_minus;

  static const double target_M_plus = 0.4824;
  static const double target_M_minus = 0.4824;

  static const double par_P_plus[3] = { 0.0, 0.114, 0.0 };
  static const double par_P_minus[3] = { 0.0, -0.114, 0.0 };
  static const double par_S_plus[3] = { 0.0, 0.0, 0.0 };
  static const double par_S_minus[3] = { 0.0, 0.0, 0.0 };
  static const double center_offset[3] = { 0.0, 0.0, 0.00123 };

  static const double initial_lapse_psi_exponent = -2.0;

  static const int Newton_maxit = 5;

  static const int npoints_A = 20;
  static const int npoints_B = 20;
  static const int npoints_phi = 12;

#endif

}
#endif
