/* TwoPunctures:  File  "TwoPunctures.c"*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <ctype.h>
#include "cctk.h"
#include "TP_utilities.h"
#include "TwoPunctures.h"
#include "params.h"



/* Swap two variables */
static inline
void swap (CCTK_REAL *  const a, CCTK_REAL *  const b)
{
  CCTK_REAL const t = *a; *a=*b; *b=t;
}
#undef SWAP
#define SWAP(a,b) (swap(&(a),&(b)))

inline double EXTEND(double M, double r) {
  return ( M * (3./8 * pow(r, 4) / pow(TPID::TP_Extend_Radius, 5) -
                 5./4 * pow(r, 2) / pow(TPID::TP_Extend_Radius, 3) +
                 15./8 / TPID::TP_Extend_Radius));
}

static
void set_initial_guess(CCTKGH *cctkGH, derivs v)
{

  int nvar = 1, n1 = TPID::npoints_A, n2 = TPID::npoints_B, n3 = TPID::npoints_phi;

  CCTK_REAL *s_x, *s_y, *s_z;
  CCTK_REAL al, A, Am1, be, B, phi, R, r, X;
  CCTK_INT ivar, i, j, k, i3D, indx;
  derivs U;
  FILE *debug_file;

  if (TPID::solve_momentum_constraint)
    nvar = 4;

  s_x    = (CCTK_REAL *) calloc(n1*n2*n3, sizeof(CCTK_REAL));
  s_y    = (CCTK_REAL *) calloc(n1*n2*n3, sizeof(CCTK_REAL));
  s_z    = (CCTK_REAL *) calloc(n1*n2*n3, sizeof(CCTK_REAL));
  allocate_derivs (&U, nvar);
  for (ivar = 0; ivar < nvar; ivar++)
    for (i = 0; i < n1; i++)
      for (j = 0; j < n2; j++)
        for (k = 0; k < n3; k++)
        {
          i3D = Index(ivar,i,j,k,1,n1,n2,n3);

          al = Pih * (2 * i + 1) / n1;
          A = -cos (al);
          be = Pih * (2 * j + 1) / n2;
          B = -cos (be);
          phi = 2. * Pi * k / n3;

          /* Calculation of (X,R)*/
          AB_To_XR (nvar, A, B, &X, &R, U);
          /* Calculation of (x,r)*/
          C_To_c (nvar, X, R, &(s_x[i3D]), &r, U);
          /* Calculation of (y,z)*/
          rx3_To_xyz (nvar, s_x[i3D], r, phi, &(s_y[i3D]), &(s_z[i3D]), U);
        }
  printf("routine Set_Initial_Guess_for_u commented out. no sources.\n");
  //Set_Initial_Guess_for_u(cctkGH, n1*n2*n3, v.d0, s_x, s_y, s_z);
  for (ivar = 0; ivar < nvar; ivar++)
    for (i = 0; i < n1; i++)
      for (j = 0; j < n2; j++)
        for (k = 0; k < n3; k++)
        {
          indx = Index(ivar,i,j,k,1,n1,n2,n3);
          v.d0[indx]/=(-cos(Pih * (2 * i + 1) / n1)-1.0);
        }
  Derivatives_AB3 (nvar, n1, n2, n3, v);
//if (do_initial_debug_output && CCTK_MyProc(cctkGH) == 0)
  if (TPID::do_initial_debug_output == 0)
  {
    debug_file=fopen("initial.dat", "w");
    assert(debug_file);
    for (ivar = 0; ivar < nvar; ivar++)
      for (i = 0; i < n1; i++)
        for (j = 0; j < n2; j++)
        {
          al = Pih * (2 * i + 1) / n1;
          A = -cos (al);
          Am1 = A -1.0;
          be = Pih * (2 * j + 1) / n2;
          B = -cos (be);
          phi = 0.0;
          indx = Index(ivar,i,j,0,1,n1,n2,n3);
          U.d0[0] = Am1 * v.d0[indx];        /* U*/
          U.d1[0] = v.d0[indx] + Am1 * v.d1[indx];        /* U_A*/
          U.d2[0] = Am1 * v.d2[indx];        /* U_B*/
          U.d3[0] = Am1 * v.d3[indx];        /* U_3*/
          U.d11[0] = 2 * v.d1[indx] + Am1 * v.d11[indx];        /* U_AA*/
          U.d12[0] = v.d2[indx] + Am1 * v.d12[indx];        /* U_AB*/
          U.d13[0] = v.d3[indx] + Am1 * v.d13[indx];        /* U_AB*/
          U.d22[0] = Am1 * v.d22[indx];        /* U_BB*/
          U.d23[0] = Am1 * v.d23[indx];        /* U_B3*/
          U.d33[0] = Am1 * v.d33[indx];        /* U_33*/
        /* Calculation of (X,R)*/
        AB_To_XR (nvar, A, B, &X, &R, U);
        /* Calculation of (x,r)*/
        C_To_c (nvar, X, R, &(s_x[indx]), &r, U);
        /* Calculation of (y,z)*/
        rx3_To_xyz (nvar, s_x[i3D], r, phi, &(s_y[indx]), &(s_z[indx]), U);
        fprintf(debug_file,
                "%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g "
                "%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
                (double)s_x[indx], (double)s_y[indx],
                (double)A,(double)B,
                (double)U.d0[0],
                (double)(-cos(Pih * (2 * i + 1) / n1)-1.0),
                (double)U.d1[0],
                (double)U.d2[0],
                (double)U.d3[0],
                (double)U.d11[0],
                (double)U.d22[0],
                (double)U.d33[0],
                (double)v.d0[indx],
                (double)v.d1[indx],
                (double)v.d2[indx],
                (double)v.d3[indx],
                (double)v.d11[indx],
                (double)v.d22[indx],
                (double)v.d33[indx]
                );
        }
    fprintf(debug_file, "\n\n");
    for (i=n2-10; i<n2; i++)
    {
      CCTK_REAL d;
      indx = Index(0,0,i,0,1,n1,n2,n3);
      d = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v,
              s_x[indx], 0.0, 0.0);
      fprintf(debug_file, "%.16g %.16g\n",
                (double)s_x[indx], (double)d);
    }
    fprintf(debug_file, "\n\n");
    for (i=n2-10; i<n2-1; i++)
    {
      CCTK_REAL d;
      int ip= Index(0,0,i+1,0,1,n1,n2,n3);
      indx = Index(0,0,i,0,1,n1,n2,n3);
      for (j=-10; j<10; j++)
      {
        d = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v,
                s_x[indx]+(s_x[ip]-s_x[indx])*j/10,
                0.0, 0.0);
        fprintf(debug_file, "%.16g %.16g\n",
                (double)(s_x[indx]+(s_x[ip]-s_x[indx])*j/10), (double)d);
      }
    }
    fprintf(debug_file, "\n\n");
    for (i = 0; i < n1; i++)
      for (j = 0; j < n2; j++)
      {
        X = 2*(2.0*i/n1-1.0);
        R = 2*(1.0*j/n2);
        if (X*X+R*R > 1.0)
        {
          C_To_c (nvar, X, R, &(s_x[indx]), &r, U);
          rx3_To_xyz (nvar, s_x[i3D], r, phi, &(s_y[indx]), &(s_z[indx]), U);
          *U.d0  = s_x[indx]*s_x[indx];
          *U.d1  = 2*s_x[indx];
          *U.d2  = 0.0;
          *U.d3 = 0.0;
          *U.d11 = 2.0;
          *U.d22 = 0.0;
          *U.d33 = *U.d12 = *U.d23 = *U.d13 = 0.0;
          C_To_c (nvar, X, R, &(s_x[indx]), &r, U);
          fprintf(debug_file,
                  "%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
                  (double)s_x[indx], (double)r, (double)X, (double)R, (double)U.d0[0],
                  (double)U.d1[0],
                  (double)U.d2[0],
                  (double)U.d3[0],
                  (double)U.d11[0],
                  (double)U.d22[0],
                  (double)U.d33[0]);
        }
      }
    fclose(debug_file);
  }
  free(s_z);
  free(s_y);
  free(s_x);
  free_derivs (&U, nvar);
  /*exit(0);*/
}

/* -------------------------------------------------------------------*/
void TwoPunctures(CCTKGH *cctkGH,
                  double *mp, double *mm, double *mp_adm, double *mm_adm,
                  double *E, double *J1, double *J2, double *J3)
{

  * mp = TPID::par_m_plus;
  * mm = TPID::par_m_minus;

  enum GRID_SETUP_METHOD { GSM_Taylor_expansion, GSM_evaluation };
  enum GRID_SETUP_METHOD gsm;

  int antisymmetric_lapse, averaged_lapse, pmn_lapse, brownsville_lapse;

  int const nvar = 1, n1 = TPID::npoints_A, n2 = TPID::npoints_B, n3 = TPID::npoints_phi;

  int const ntotal = n1 * n2 * n3 * nvar;
#if 1
  int percent10 = 0;
#endif
  static CCTK_REAL *F = NULL;
  static derivs u, v, cf_v;
  CCTK_REAL admMass;

  if (! F) {
    CCTK_REAL up, um;
    /* Solve only when called for the first time */
    F = dvector (0, ntotal - 1);
    allocate_derivs (&u, ntotal);
    allocate_derivs (&v, ntotal);
    allocate_derivs (&cf_v, ntotal);

    printf("INFO: b = %g\n", TPID::par_b);

    /* initialise to 0 */
    for (int j = 0; j < ntotal; j++)
    {
      cf_v.d0[j] = 0.0;
      cf_v.d1[j] = 0.0;
      cf_v.d2[j] = 0.0;
      cf_v.d3[j] = 0.0;
      cf_v.d11[j] = 0.0;
      cf_v.d12[j] = 0.0;
      cf_v.d13[j] = 0.0;
      cf_v.d22[j] = 0.0;
      cf_v.d23[j] = 0.0;
      cf_v.d33[j] = 0.0;
      v.d0[j] = 0.0;
      v.d1[j] = 0.0;
      v.d2[j] = 0.0;
      v.d3[j] = 0.0;
      v.d11[j] = 0.0;
      v.d12[j] = 0.0;
      v.d13[j] = 0.0;
      v.d22[j] = 0.0;
      v.d23[j] = 0.0;
      v.d33[j] = 0.0;
    }
    /* call for external initial guess */
    if (TPID::use_external_initial_guess)
    {
      set_initial_guess(cctkGH, v);
    }

    /* If bare masses are not given, iteratively solve for them given the
       target ADM masses target_M_plus and target_M_minus and with initial
       guesses given by par_m_plus and par_m_minus. */
    if(!(TPID::give_bare_mass)) {
      CCTK_REAL tmp, mp_adm_err, mm_adm_err;
      char valbuf[100];

      CCTK_REAL M_p = TPID::target_M_plus;
      CCTK_REAL M_m = TPID::target_M_minus;

      printf("Attempting to find bare masses.\n");
      printf("Target ADM masses: M_p=%g and M_m=%g\n",
                  (double) M_p, (double) M_m);
      printf("ADM mass tolerance: %g\n", (double) TPID::adm_tol);

      /* Loop until both ADM masses are within adm_tol of their target */
      do {
        printf("Bare masses: mp=%.15g, mm=%.15g\n",
                    (double)*mp, (double)*mm);
        Newton (cctkGH, nvar, n1, n2, n3, v, TPID::Newton_tol, 1);

        F_of_v (cctkGH, nvar, n1, n2, n3, v, F, u);

        up = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v, TPID::par_b, 0., 0.);
        um = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v,-TPID::par_b, 0., 0.);

        /* Calculate the ADM masses from the current bare mass guess */
        *mp_adm = (1 + up) * *mp + *mp * *mm / (4. * TPID::par_b);
        *mm_adm = (1 + um) * *mm + *mp * *mm / (4. * TPID::par_b);

        /* Check how far the current ADM masses are from the target */
        mp_adm_err = fabs(M_p-*mp_adm);
        mm_adm_err = fabs(M_m-*mm_adm);
        printf("ADM mass error: M_p_err=%.15g, M_m_err=%.15g\n",
                    mp_adm_err, mm_adm_err);

        /* Invert the ADM mass equation and update the bare mass guess so that
           it gives the correct target ADM masses */
        tmp = -4*TPID::par_b*( 1 + um + up + um*up ) +
                sqrt(16*TPID::par_b*M_m*(1 + um)*(1 + up) +
                  pow(-M_m + M_p + 4*TPID::par_b*(1 + um)*(1 + up),2));
        *mp = (tmp + M_p - M_m)/(2.*(1 + up));
        *mm = (tmp - M_p + M_m)/(2.*(1 + um));

        /* Set the par_m_plus and par_m_minus parameters */
/*
        sprintf (valbuf,"%.17g", (double) *mp);
        CCTK_ParameterSet ("par_m_plus", "TwoPunctures", valbuf);
        sprintf (valbuf,"%.17g", (double) *mm);
        CCTK_ParameterSet ("par_m_minus", "TwoPunctures", valbuf);
*/

        TPID::par_m_plus = *mp;
        TPID::par_m_minus = *mm;

      } while ( (mp_adm_err > TPID::adm_tol) ||
                (mm_adm_err > TPID::adm_tol) );

      printf("Found bare masses.");
    }

    Newton (cctkGH, nvar, n1, n2, n3, v, TPID::Newton_tol, TPID::Newton_maxit);

    F_of_v (cctkGH, nvar, n1, n2, n3, v, F, u);

    SpecCoef(n1, n2, n3, 0, v.d0, cf_v.d0);

    printf("The two puncture masses are mp=%.17g and mm=%.17g\n", *mp, *mm);

    up = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v, TPID::par_b, 0., 0.);
    um = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v,-TPID::par_b, 0., 0.);

    /* Calculate the ADM masses from the current bare mass guess */
    *mp_adm = (1 + up) * *mp + *mp * *mm / (4. * TPID::par_b);
    *mm_adm = (1 + um) * *mm + *mp * *mm / (4. * TPID::par_b);

    printf("Puncture 1 ADM mass is %g\n", *mp_adm);
    printf("Puncture 2 ADM mass is %g\n", *mm_adm);

    /* print out ADM mass, eq.: \Delta M_ADM=2*r*u=4*b*V for A=1,B=0,phi=0 */
    admMass = (*mp + *mm
               - 4*TPID::par_b*PunctEvalAtArbitPosition(v.d0, 0, 1, 0, 0, nvar, n1, n2, n3));
    printf("The total ADM mass is %g\n", admMass);
    *E = admMass;

    /*
      Run this in Mathematica (version 8 or later) with
        math -script <file>

      Needs["SymbolicC`"];
      co = Table["center_offset[" <> ToString[i] <> "]", {i, 0, 2}];
      r1 = co + {"par_b", 0, 0};
      r2 = co + {-"par_b", 0, 0};
      {p1, p2} = Table["par_P_" <> bh <> "[" <> ToString[i] <> "]", {bh, {"plus", "minus"}}, {i, 0, 2}];
      {s1, s2} = Table["par_S_" <> bh <> "[" <> ToString[i] <> "]", {bh, {"plus", "minus"}}, {i, 0, 2}];

      J = Cross[r1, p1] + Cross[r2, p2] + s1 + s2;

      JVar = Table["*J" <> ToString[i], {i, 1, 3}];
      Print[OutputForm@StringReplace[
        ToCCodeString@MapThread[CAssign[#1, CExpression[#2]] &, {JVar, J}],
        "\"" -> ""]];
     */

    *J1 = -(TPID::center_offset[2]*TPID::par_P_minus[1]) + TPID::center_offset[1]*TPID::par_P_minus[2] - TPID::center_offset[2]*TPID::par_P_plus[1] + TPID::center_offset[1]*TPID::par_P_plus[2] + TPID::par_S_minus[0] + TPID::par_S_plus[0];
    *J2 = TPID::center_offset[2]*TPID::par_P_minus[0] - TPID::center_offset[0]*TPID::par_P_minus[2] + TPID::par_b*TPID::par_P_minus[2] + TPID::center_offset[2]*TPID::par_P_plus[0] - TPID::center_offset[0]*TPID::par_P_plus[2] - TPID::par_b*TPID::par_P_plus[2] + TPID::par_S_minus[1] + TPID::par_S_plus[1];
    *J3 = -(TPID::center_offset[1]*TPID::par_P_minus[0]) + TPID::center_offset[0]*TPID::par_P_minus[1] - TPID::par_b*TPID::par_P_minus[1] - TPID::center_offset[1]*TPID::par_P_plus[0] + TPID::center_offset[0]*TPID::par_P_plus[1] + TPID::par_b*TPID::par_P_plus[1] + TPID::par_S_minus[2] + TPID::par_S_plus[2];
  }

  if (TPID::grid_setup_method == TAYLOR_EXPANSION)
  {
    gsm = GSM_Taylor_expansion;
  }
  else if (TPID::grid_setup_method == EVALUATION)
  {
    gsm = GSM_evaluation;
  }
  else
  {
    printf("internal error. unknown grid_setup_method = %d\n",
            TPID::grid_setup_method);
    exit(-1);
  }

  antisymmetric_lapse = (TPID::initial_lapse == LAPSE_ANTISYMMETRIC) ? 1 : 0;
  averaged_lapse = (TPID::initial_lapse == LAPSE_AVERAGED) ? 1 : 0;
	pmn_lapse = (TPID::initial_lapse == LAPSE_PSIN) ? 1 : 0 ;
  if (pmn_lapse)
		printf("Setting initial lapse to psi^%f profile.\n",
               TPID::initial_lapse_psi_exponent);
  brownsville_lapse = (TPID::initial_lapse == LAPSE_BROWNSVILLE) ? 1 : 0;
  if (brownsville_lapse)
    printf( "Setting initial lapse to a Brownsville-style profile "
            "with exp %f.", TPID::initial_lapse_psi_exponent);

  printf("Interpolating result.\n");

 /* Apparently in Cactus, you can choose to save one of the following
  * to 3D arrays:
  *    psi
  *    psi + first derivatives
  *    psi + first derivatives + second derivatives.
  *
  * For now I only provide the option to save psi.
  */

  int metric_type = MT_STANDARD;
  int conformal_storage = CS_FACTOR;
  int conformal_state;

  if (metric_type == MT_STATIC_CONFORMAL) {
    if (conformal_storage == CS_FACTOR) {
      conformal_state = 1;
    } else if (conformal_storage == CS_FACTOR_DERIVS) {
      conformal_state = 2;
    } else if (conformal_storage == CS_FACTOR_SECOND_DERIVS) {
      conformal_state = 3;
    }
  } else {
    conformal_state = 0;
  }

  int nx = cctkGH->cctk_lsh[0];
  int ny = cctkGH->cctk_lsh[1];
  int nz = cctkGH->cctk_lsh[2];

#pragma omp parallel for
  for (int k = 0; k < nz; k++)
  {
    for (int j = 0; j < ny; j++)
    {
      for (int i = 0; i < nx; i++)
      {
#if 1
        /* We can't output this when running in parallel */
        if (percent10 != 10*(i+j*cctkGH->cctk_lsh[0]+k*cctkGH->cctk_lsh[0]*cctkGH->cctk_lsh[1]) /
                            (cctkGH->cctk_lsh[0] * cctkGH->cctk_lsh[1] * cctkGH->cctk_lsh[2]))
        {
            percent10 = 10*(i+j*cctkGH->cctk_lsh[0]+k*cctkGH->cctk_lsh[0]*cctkGH->cctk_lsh[1]) /
                           (cctkGH->cctk_lsh[0] * cctkGH->cctk_lsh[1] * cctkGH->cctk_lsh[2]);
            printf("... %3d%% done", percent10*10);
        }
#endif

        const int ind = i + nx*(j + ny*k);

        CCTK_REAL xx, yy, zz;
        xx = cctkGH->x[ind] - TPID::center_offset[0];
        yy = cctkGH->y[ind] - TPID::center_offset[1];
        zz = cctkGH->z[ind] - TPID::center_offset[2];

        /* We implement swapping the x and z coordinates as follows.
           The bulk of the code that performs the actual calculations
           is unchanged.  This code looks only at local variables.
           Before the bulk --i.e., here-- we swap all x and z tensor
           components, and after the code --i.e., at the end of this
           main loop-- we swap everything back.  */
        if (TPID::swap_xz) {
          /* Swap the x and z coordinates */
          SWAP (xx, zz);
        }

        CCTK_REAL r_plus
          = sqrt(pow(xx - TPID::par_b, 2) + pow(yy, 2) + pow(zz, 2));
        CCTK_REAL r_minus
          = sqrt(pow(xx + TPID::par_b, 2) + pow(yy, 2) + pow(zz, 2));

        CCTK_REAL U;
        switch (gsm)
        {
        case GSM_Taylor_expansion:
          U = PunctTaylorExpandAtArbitPosition
            (0, nvar, n1, n2, n3, v, xx, yy, zz);
          break;
        case GSM_evaluation:
          U = PunctIntPolAtArbitPositionFast(0, nvar, n1, n2, n3, cf_v, xx, yy, zz);
          break;
        default:
          assert (0);
        }
        r_plus = pow (pow (r_plus, 4) + pow (TPID::TP_epsilon, 4), 0.25);
        r_minus = pow (pow (r_minus, 4) + pow (TPID::TP_epsilon, 4), 0.25);
        if (r_plus < TPID::TP_Tiny)
            r_plus = TPID::TP_Tiny;
        if (r_minus < TPID::TP_Tiny)
            r_minus = TPID::TP_Tiny;
        CCTK_REAL psi1 = 1
          + 0.5 * *mp / r_plus
          + 0.5 * *mm / r_minus + U;
/*
#define EXTEND(M,r) \
          ( M * (3./8 * pow(r, 4) / pow(TPID::TP_Extend_Radius, 5) - \
                 5./4 * pow(r, 2) / pow(TPID::TP_Extend_Radius, 3) + \
                 15./8 / TPID::TP_Extend_Radius))
*/
        if (r_plus < TPID::TP_Extend_Radius) {
          psi1 = 1
             + 0.5 * EXTEND(*mp,r_plus)
             + 0.5 * *mm / r_minus + U;
        }
        if (r_minus < TPID::TP_Extend_Radius) {
          psi1 = 1
             + 0.5 * EXTEND(*mm,r_minus)
             + 0.5 * *mp / r_plus + U;
        }
        CCTK_REAL static_psi = 1;

        CCTK_REAL Aij[3][3];
        BY_Aijofxyz (xx, yy, zz, Aij);

        CCTK_REAL old_alp=1.0;
        if (TPID::multiply_old_lapse)
            old_alp = cctkGH->alp[ind];

        if ((conformal_state > 0) || (pmn_lapse) || (brownsville_lapse)) {

          CCTK_REAL xp, yp, zp, rp, ir;
          CCTK_REAL s1, s3, s5;
          CCTK_REAL p, px, py, pz, pxx, pxy, pxz, pyy, pyz, pzz;
          p = 1.0;
          px = py = pz = 0.0;
          pxx = pxy = pxz = 0.0;
          pyy = pyz = pzz = 0.0;

          /* first puncture */
          xp = xx - TPID::par_b;
          yp = yy;
          zp = zz;
          rp = sqrt (xp*xp + yp*yp + zp*zp);
          rp = pow (pow (rp, 4) + pow (TPID::TP_epsilon, 4), 0.25);
          if (rp < TPID::TP_Tiny)
              rp = TPID::TP_Tiny;
          ir = 1.0/rp;

          if (rp < TPID::TP_Extend_Radius) {
            ir = EXTEND(1., rp);
          }

          s1 = 0.5* *mp *ir;
          s3 = -s1*ir*ir;
          s5 = -3.0*s3*ir*ir;

          p += s1;

          px += xp*s3;
          py += yp*s3;
          pz += zp*s3;

          pxx += xp*xp*s5 + s3;
          pxy += xp*yp*s5;
          pxz += xp*zp*s5;
          pyy += yp*yp*s5 + s3;
          pyz += yp*zp*s5;
          pzz += zp*zp*s5 + s3;

          /* second puncture */
          xp = xx + TPID::par_b;
          yp = yy;
          zp = zz;
          rp = sqrt (xp*xp + yp*yp + zp*zp);
          rp = pow (pow (rp, 4) + pow (TPID::TP_epsilon, 4), 0.25);
          if (rp < TPID::TP_Tiny)
              rp = TPID::TP_Tiny;
          ir = 1.0/rp;

          if (rp < TPID::TP_Extend_Radius) {
            ir = EXTEND(1., rp);
          }

          s1 = 0.5* *mm *ir;
          s3 = -s1*ir*ir;
          s5 = -3.0*s3*ir*ir;

          p += s1;

          px += xp*s3;
          py += yp*s3;
          pz += zp*s3;

          pxx += xp*xp*s5 + s3;
          pxy += xp*yp*s5;
          pxz += xp*zp*s5;
          pyy += yp*yp*s5 + s3;
          pyz += yp*zp*s5;
          pzz += zp*zp*s5 + s3;

          if (conformal_state >= 1) {
            static_psi = p;
            cctkGH->psi[ind] = static_psi;
          }
          if (conformal_state >= 2) {
            printf("Code doesn't yet work for conformal_state == 2.\n");
/*          psix[ind] = px / static_psi;
            psiy[ind] = py / static_psi;
            psiz[ind] = pz / static_psi;
*/
          }
          if (conformal_state >= 3) {
            printf("Code doesn't yet work for conformal_state == 3.\n");
/*
            psixx[ind] = pxx / static_psi;
            psixy[ind] = pxy / static_psi;
            psixz[ind] = pxz / static_psi;
            psiyy[ind] = pyy / static_psi;
            psiyz[ind] = pyz / static_psi;
            psizz[ind] = pzz / static_psi;
*/
          }

          if (pmn_lapse)
            cctkGH->alp[ind] = pow(p, TPID::initial_lapse_psi_exponent);
          if (brownsville_lapse)
            cctkGH->alp[ind] = 2.0/(1.0+pow(p, TPID::initial_lapse_psi_exponent));

        } /* if conformal-state > 0 */

        cctkGH->puncture_u[ind] = U;

        cctkGH->gxx[ind] = pow (psi1 / static_psi, 4);
        cctkGH->gxy[ind] = 0;
        cctkGH->gxz[ind] = 0;
        cctkGH->gyy[ind] = pow (psi1 / static_psi, 4);
        cctkGH->gyz[ind] = 0;
        cctkGH->gzz[ind] = pow (psi1 / static_psi, 4);

        cctkGH->kxx[ind] = Aij[0][0] / pow(psi1, 2);
        cctkGH->kxy[ind] = Aij[0][1] / pow(psi1, 2);
        cctkGH->kxz[ind] = Aij[0][2] / pow(psi1, 2);
        cctkGH->kyy[ind] = Aij[1][1] / pow(psi1, 2);
        cctkGH->kyz[ind] = Aij[1][2] / pow(psi1, 2);
        cctkGH->kzz[ind] = Aij[2][2] / pow(psi1, 2);

        if (antisymmetric_lapse || averaged_lapse) {
          cctkGH->alp[ind] =
            ((1.0 -0.5* *mp /r_plus -0.5* *mm/r_minus)
            /(1.0 +0.5* *mp /r_plus +0.5* *mm/r_minus));

          if (r_plus < TPID::TP_Extend_Radius) {
            cctkGH->alp[ind] =
              ((1.0 -0.5*EXTEND(*mp, r_plus) -0.5* *mm/r_minus)
              /(1.0 +0.5*EXTEND(*mp, r_plus) +0.5* *mm/r_minus));
          }
          if (r_minus < TPID::TP_Extend_Radius) {
            cctkGH->alp[ind] =
              ((1.0 -0.5*EXTEND(*mm, r_minus) -0.5* *mp/r_plus)
              /(1.0 +0.5*EXTEND(*mp, r_minus) +0.5* *mp/r_plus));
          }

          if (averaged_lapse) {
            cctkGH->alp[ind] = 0.5 * (1.0 + cctkGH->alp[ind]);
          }
        }

        double gd[3][3];
        gd[0][0] = pow (psi1 / static_psi, 4);
        gd[0][1] = 0.0;
        gd[0][2] = 0.0;
        gd[1][0] = 0.0;
        gd[1][1] = pow (psi1 / static_psi, 4);
        gd[1][2] = 0.0;
        gd[2][0] = 0.0;
        gd[2][1] = 0.0;
        gd[2][2] = pow (psi1 / static_psi, 4);

        double Kd[3][3];
        Kd[0][0] = Aij[0][0] / pow(psi1, 2);
        Kd[0][1] = Aij[0][1] / pow(psi1, 2);
        Kd[0][2] = Aij[0][2] / pow(psi1, 2);
        Kd[1][0] = Kd[0][1];
        Kd[1][1] = Aij[1][1] / pow(psi1, 2);
        Kd[1][2] = Aij[1][2] / pow(psi1, 2);
        Kd[2][0] = Kd[0][2];
        Kd[2][1] = Kd[1][2];
        Kd[2][2] = Aij[2][2] / pow(psi1, 2);

        double gu[3][3], gtd[0][0], Atd[3][3];
        double detgd, idetgd, trK;
        double chi;
        double t1, t2, t4, t6, t7, t9, t10, t12, t16;

#include "adm2bssn.h"

        cctkGH->gtxx[ind] = gtd[0][0];
        cctkGH->gtxy[ind] = gtd[0][1];
        cctkGH->gtxz[ind] = gtd[0][2];
        cctkGH->gtyy[ind] = gtd[1][1];
        cctkGH->gtyz[ind] = gtd[1][2];
        cctkGH->gtzz[ind] = gtd[2][2];

        cctkGH->Atxx[ind] = Atd[0][0];
        cctkGH->Atxy[ind] = Atd[0][1];
        cctkGH->Atxz[ind] = Atd[0][2];
        cctkGH->Atyy[ind] = Atd[1][1];
        cctkGH->Atyz[ind] = Atd[1][2];
        cctkGH->Atzz[ind] = Atd[2][2];

        cctkGH->trK[ind] = trK;
        cctkGH->chi[ind] = chi;

        if (TPID::multiply_old_lapse)
          cctkGH->alp[ind] *= old_alp;

        if (TPID::swap_xz) {
          /* Swap the x and z components of all tensors */
/*
          if (*conformal_state >= 2) {
            SWAP (psix[ind], psiz[ind]);
          }
          if (*conformal_state >= 3) {
            SWAP (psixx[ind], psizz[ind]);
            SWAP (psixy[ind], psiyz[ind]);
          }
*/
          SWAP (cctkGH->gxx[ind], cctkGH->gzz[ind]);
          SWAP (cctkGH->gxy[ind], cctkGH->gyz[ind]);
          SWAP (cctkGH->kxx[ind], cctkGH->kzz[ind]);
          SWAP (cctkGH->kxy[ind], cctkGH->kyz[ind]);
        } /* if swap_xz */

      } /* for i */
    }   /* for j */
  }     /* for k */

  if (0) {
    /* Keep the result around for the next time */
    free_dvector (F, 0, ntotal - 1);
    free_derivs (&u, ntotal);
    free_derivs (&v, ntotal);
    free_derivs (&cf_v, ntotal);
  }
}
