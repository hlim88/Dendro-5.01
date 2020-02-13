#ifndef _HAVE_PARAMS_H
#define _HAVE_PARAMS_H

  extern const int     verbose;

  extern const int     initial_lapse;
  extern const int     give_bare_mass;
  extern const int     grid_setup_method;
  extern const int     swap_xz;

  extern const int     use_sources;
  extern const int     rescale_sources;
  extern const int     use_external_initial_guess;
  extern const int     do_residuum_debug_output;
  extern const int     do_initial_debug_output;
  extern const int     multiply_old_lapse;
  extern const int     solve_momentum_constraint;

  extern const double  adm_tol;
  extern const double  Newton_tol;
  extern const double  TP_epsilon;
  extern const double  TP_Tiny;
  extern const double  TP_Extend_Radius;
  extern const double  par_b;
  extern double  par_m_plus;
  extern double  par_m_minus;
  extern const double  target_M_plus;
  extern const double  target_M_minus;
  extern const double  par_P_plus[3];
  extern const double  par_P_minus[3];
  extern const double  par_S_plus[3];
  extern const double  par_S_minus[3];
  extern const double  center_offset[3];
  extern const double  initial_lapse_psi_exponent;

  extern const int     Newton_maxit;

  extern const int     npoints_A;
  extern const int     npoints_B;
  extern const int     npoints_phi;

  typedef struct {
  
  } Grid;

#endif
