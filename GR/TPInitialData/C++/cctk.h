#ifndef _HAVE_CCTK_H
#define _HAVE_CCTK_H

#define ID_PARS 1

#define CCTK_REAL double
#define CCTK_INT int
#define CCTK_ERROR code_exit


#define TAYLOR_EXPANSION 0
#define EVALUATION 1

#define LAPSE_ANTISYMMETRIC 0
#define LAPSE_AVERAGED 1
#define LAPSE_PSIN 2
#define LAPSE_BROWNSVILLE 3

#define CS_FACTOR 0
#define CS_FACTOR_DERIVS 1
#define CS_FACTOR_SECOND_DERIVS 2

#define MT_STATIC_CONFORMAL 0
#define MT_STANDARD 1


typedef struct {

  int cctk_lsh[3];

  double *puncture_u;
  double *gxx;
  double *gxy;
  double *gxz;
  double *gyy;
  double *gyz;
  double *gzz;
  double *kxx;
  double *kxy;
  double *kxz;
  double *kyy;
  double *kyz;
  double *kzz;
  double *alp;
  double *psi;

  double *gtxx;
  double *gtxy;
  double *gtxz;
  double *gtyy;
  double *gtyz;
  double *gtzz;

  double *Atxx;
  double *Atxy;
  double *Atxz;
  double *Atyy;
  double *Atyz;
  double *Atzz;

  double *chi;
  double *trK;

  double *hAtxx;
  double *hAtxy;
  double *hAtxz;
  double *hAtyy;
  double *hAtyz;
  double *hAtzz;

  double *hchi;
  double *hAlpha;


  double *x, *y, *z;

} CCTKGH;

#define CCTK_POINTER_TO_CONST CCTKGH *

void code_exit(const char *);
void TwoPunctures(CCTKGH *cctkGH, 
                  double *mp, double *mm, double *mp_adm, double *mm_adm,
                  double *E, double *J1, double *J2, double *J3);
void init_data_puncture(const double *coord, double *u);


#define U_ALPHA 0
#define U_SHIFTX 1
#define U_SHIFTY 2
#define U_SHIFTZ 3
#define U_CHI 4
#define U_TRK 5
#define U_GAMTX 6
#define U_GAMTY 7
#define U_GAMTZ 8
#define U_GBX 9
#define U_GBY 10
#define U_GBZ 11
#define U_GTXX 12
#define U_GTXY 13
#define U_GTXZ 14
#define U_GTYY 15
#define U_GTYZ 16
#define U_GTZZ 17
#define U_ATXX 18
#define U_ATXY 19
#define U_ATXZ 20
#define U_ATYY 21
#define U_ATYZ 22
#define U_ATZZ 23

#endif
