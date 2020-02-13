#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cctk.h"
#include "params.h"


/* The initial data parameters are hard-coded. The options are
       0 =  GW150914 data
       1 =  Data for ET TPID sample solution
       2 =  Data for comparison to HAD's approximate solution.
            (no spin currently)
*/



void output_line_x(const char* fname, double *f, double *x1d, int *shp);
void output_line_y(const char* fname, double *f, double *y1d, int *shp);
void approx_sol(CCTKGH *cctkGH);


/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
int main(int argc, char **argv)
{

 /* create Cartesian grid for output */
  const int ntot = TPID::nx0 * TPID::ny0 * TPID::nz0;
  double dx =  (TPID::xmax - TPID::xmin)/(TPID::nx0-1);
  double dy =  (TPID::ymax - TPID::ymin)/(TPID::ny0-1);
  double dz =  (TPID::zmax - TPID::zmin)/(TPID::nz0-1);

  double *x1d = (double *) malloc(sizeof(double) * TPID::nx0);
  double *y1d = (double *) malloc(sizeof(double) * TPID::ny0);
  double *z1d = (double *) malloc(sizeof(double) * TPID::nz0);

  for (int i = 0; i < TPID::nx0; i++) {
    x1d[i] = TPID::xmin + i*dx;
  }
  for (int i = 0; i < TPID::ny0; i++) {
    y1d[i] = TPID::ymin + i*dy;
  }
  for (int i = 0; i < TPID::nz0; i++) {
    z1d[i] = TPID::zmin + i*dz;
  }

  double *x = (double *) malloc(sizeof(double) * ntot);
  double *y = (double *) malloc(sizeof(double) * ntot);
  double *z = (double *) malloc(sizeof(double) * ntot);

 /* Cactus uses 3D coordinate arrays */
  for (int k = 0; k < TPID::nz0; k++) {
    for (int j = 0; j < TPID::ny0; j++) {
      for (int i = 0; i < TPID::nx0; i++) {
        int pp = i + (k*TPID::ny0 + j)*TPID::nx0;
        x[pp] = x1d[i];
        y[pp] = y1d[j];
        z[pp] = z1d[k];
      }
    }
  }


  /* Storage for GR variables */
  CCTKGH cctkGH;

  double *g11 = (double *) malloc(sizeof(double) * ntot);
  double *g12 = (double *) malloc(sizeof(double) * ntot);
  double *g13 = (double *) malloc(sizeof(double) * ntot);
  double *g22 = (double *) malloc(sizeof(double) * ntot);
  double *g23 = (double *) malloc(sizeof(double) * ntot);
  double *g33 = (double *) malloc(sizeof(double) * ntot);

  double *psi = (double *) malloc(sizeof(double) * ntot);
  double *u= (double *) malloc(sizeof(double) * ntot);

  double *K11 = (double *) malloc(sizeof(double) * ntot);
  double *K12 = (double *) malloc(sizeof(double) * ntot);
  double *K13 = (double *) malloc(sizeof(double) * ntot);
  double *K22 = (double *) malloc(sizeof(double) * ntot);
  double *K23 = (double *) malloc(sizeof(double) * ntot);
  double *K33 = (double *) malloc(sizeof(double) * ntot);

  double *Alpha = (double *) malloc(sizeof(double) * ntot);

  double *gt11 = (double *) malloc(sizeof(double) * ntot);
  double *gt12 = (double *) malloc(sizeof(double) * ntot);
  double *gt13 = (double *) malloc(sizeof(double) * ntot);
  double *gt22 = (double *) malloc(sizeof(double) * ntot);
  double *gt23 = (double *) malloc(sizeof(double) * ntot);
  double *gt33 = (double *) malloc(sizeof(double) * ntot);

  double *At11 = (double *) malloc(sizeof(double) * ntot);
  double *At12 = (double *) malloc(sizeof(double) * ntot);
  double *At13 = (double *) malloc(sizeof(double) * ntot);
  double *At22 = (double *) malloc(sizeof(double) * ntot);
  double *At23 = (double *) malloc(sizeof(double) * ntot);
  double *At33 = (double *) malloc(sizeof(double) * ntot);

  double *trK = (double *) malloc(sizeof(double) * ntot);
  double *chi = (double *) malloc(sizeof(double) * ntot);

  double *hAlpha = (double *) malloc(sizeof(double) * ntot);
  double *hAt11 = (double *) malloc(sizeof(double) * ntot);
  double *hAt12 = (double *) malloc(sizeof(double) * ntot);
  double *hAt13 = (double *) malloc(sizeof(double) * ntot);
  double *hAt22 = (double *) malloc(sizeof(double) * ntot);
  double *hAt23 = (double *) malloc(sizeof(double) * ntot);
  double *hAt33 = (double *) malloc(sizeof(double) * ntot);
  double *hchi = (double *) malloc(sizeof(double) * ntot);


  cctkGH.gxx = g11;
  cctkGH.gxy = g12;
  cctkGH.gxz = g13;
  cctkGH.gyy = g22;
  cctkGH.gyz = g23;
  cctkGH.gzz = g33;

  cctkGH.kxx = K11;
  cctkGH.kxy = K12;
  cctkGH.kxz = K13;
  cctkGH.kyy = K22;
  cctkGH.kyz = K23;
  cctkGH.kzz = K33;

  cctkGH.alp = Alpha;

  cctkGH.psi = psi;
  cctkGH.puncture_u = u;

  cctkGH.gtxx = gt11;
  cctkGH.gtxy = gt12;
  cctkGH.gtxz = gt13;
  cctkGH.gtyy = gt22;
  cctkGH.gtyz = gt23;
  cctkGH.gtzz = gt33;

  cctkGH.Atxx = At11;
  cctkGH.Atxy = At12;
  cctkGH.Atxz = At13;
  cctkGH.Atyy = At22;
  cctkGH.Atyz = At23;
  cctkGH.Atzz = At33;

  cctkGH.chi = chi;
  cctkGH.trK = trK;

  cctkGH.hAtxx = hAt11;
  cctkGH.hAtxy = hAt12;
  cctkGH.hAtxz = hAt13;
  cctkGH.hAtyy = hAt22;
  cctkGH.hAtyz = hAt23;
  cctkGH.hAtzz = hAt33;

  cctkGH.hchi = hchi;
  cctkGH.hAlpha = hAlpha;


  cctkGH.cctk_lsh[0] = TPID::nx0;
  cctkGH.cctk_lsh[1] = TPID::ny0;
  cctkGH.cctk_lsh[2] = TPID::nz0;

  cctkGH.x = x;
  cctkGH.y = y;
  cctkGH.z = z;

 /* The solution from the TwoPuncture thorn */
  double massp, massm, massp_adm, massm_adm, E, J1, J2, J3;
  TwoPunctures(&cctkGH, &massp, &massm, &massp_adm, &massm_adm,
               &E, &J1, &J2, &J3);

 /* The approximate solution from the HAD code. */
  approx_sol(&cctkGH);

  output_line_x("u", u, x1d, cctkGH.cctk_lsh);
  output_line_x("g11", g11, x1d, cctkGH.cctk_lsh);
  output_line_x("K11", K11, x1d, cctkGH.cctk_lsh);
  output_line_x("Alpha", Alpha, x1d, cctkGH.cctk_lsh);

  output_line_y("u", u, y1d, cctkGH.cctk_lsh);
  output_line_y("g11", g11, y1d, cctkGH.cctk_lsh);
  output_line_y("K11", K11, y1d, cctkGH.cctk_lsh);
  output_line_y("Alpha", Alpha, y1d, cctkGH.cctk_lsh);

 /* BSSN vars */
  output_line_x("chi", chi, x1d, cctkGH.cctk_lsh);
  output_line_x("trK", trK, x1d, cctkGH.cctk_lsh);
  output_line_x("At11", At11, x1d, cctkGH.cctk_lsh);

  output_line_y("chi", chi, y1d, cctkGH.cctk_lsh);
  output_line_y("trK", trK, y1d, cctkGH.cctk_lsh);
  output_line_y("At11", At11, y1d, cctkGH.cctk_lsh);

 /* approx solutions */
  output_line_x("hchi", hchi, x1d, cctkGH.cctk_lsh);
  output_line_x("hAt11", hAt11, x1d, cctkGH.cctk_lsh);
  output_line_x("hAlpha", hAlpha, x1d, cctkGH.cctk_lsh);

  output_line_y("hchi", hchi, y1d, cctkGH.cctk_lsh);
  output_line_y("hAt11", hAt11, y1d, cctkGH.cctk_lsh);
  output_line_y("hAlpha", hAlpha, y1d, cctkGH.cctk_lsh);



  /* free memory */
  free(x1d);
  free(y1d);
  free(z1d);
  free(x);
  free(y);
  free(z);
  free(g11);
  free(g12);
  free(g13);
  free(g22);
  free(g23);
  free(g33);
  free(K11);
  free(K12);
  free(K13);
  free(K22);
  free(K23);
  free(K33);

  free(Alpha);

  free(psi);
  free(u);

  free(gt11);
  free(gt12);
  free(gt13);
  free(gt22);
  free(gt23);
  free(gt33);

  free(At11);
  free(At12);
  free(At13);
  free(At22);
  free(At23);
  free(At33);

  free(chi);
  free(trK);

 /* had funcs */
  free(hAt11);
  free(hAt12);
  free(hAt13);
  free(hAt22);
  free(hAt23);
  free(hAt33);
  free(hchi);
  free(hAlpha);


  printf("...finis...\n");
  return 0;
}

// code_exit {{{
/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void code_exit(const char *s)
{
  printf("%s\n",s);
  exit(-1);
}
// }}}

// output_line_x {{{
/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void output_line_x(const char* fname, double *f, double *x1d, int *shp)
{

  const int nx = shp[0];
  const int ny = shp[1];
  const int nz = shp[2];

  char filename[64];
  sprintf(filename, "%s_x.dat", fname);

  FILE *fp = fopen(filename, "w");

  int jj = (ny-1)/2;
  int kk = (nz-1)/2;

  printf("output_line_x: Writing %s to file. jj=%d, kk=%d\n",fname,jj,kk);

  for (int i = 0; i < nx; i++) {
    int pp = i + nx*(jj + ny*kk);
    fprintf(fp, "%f    %f\n",x1d[i],f[pp]);
  }

  fclose(fp);

}
// }}}

// output_line_y {{{
/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void output_line_y(const char* fname, double *f, double *y1d, int *shp)
{

  const int nx = shp[0];
  const int ny = shp[1];
  const int nz = shp[2];

  char filename[64];
  sprintf(filename, "%s_y.dat", fname);

  FILE *fp = fopen(filename, "w");

  int ii = (nx-1)/2;
  int kk = (nz-1)/2;

  printf("output_line_y: Writing %s to file. ii=%d, kk=%d\n",fname,ii,kk);

  for (int j = 0; j < ny; j++) {
    int pp = ii + nx*(j + ny*kk);
    fprintf(fp, "%f    %f\n",y1d[j],f[pp]);
  }

  fclose(fp);

}
// }}}

/*----------------------------------------------------------------------
 *
 * The approximate solution from the HAD code
 *
 *----------------------------------------------------------------------*/
void approx_sol(CCTKGH *cctkGH)
{
  int nx = cctkGH->cctk_lsh[0];
  int ny = cctkGH->cctk_lsh[1];
  int nz = cctkGH->cctk_lsh[2];

  double pos[3];
  double upt[24];

  double *hAlpha = cctkGH->hAlpha;
  double *hchi = cctkGH->hchi;
  double *hAtxx = cctkGH->hAtxx;
  double *hAtxy = cctkGH->hAtxy;
  double *hAtxz = cctkGH->hAtxz;
  double *hAtyy = cctkGH->hAtyy;
  double *hAtyz = cctkGH->hAtyz;
  double *hAtzz = cctkGH->hAtzz;

  double *x = cctkGH->x;
  double *y = cctkGH->y;
  double *z = cctkGH->z;


  for (int k = 0; k < nz; k++) {
  for (int j = 0; j < ny; j++) {
  for (int i = 0; i < nx; i++) {
    int pp = i + nx*(j + k*ny);
    pos[0] = x[pp];
    pos[1] = y[pp];
    pos[2] = z[pp];

    init_data_puncture(pos, upt);

    hAlpha[pp] = upt[U_ALPHA];
    hchi[pp] = upt[U_CHI];
    hAtxx[pp] = upt[U_ATXX];
    hAtxy[pp] = upt[U_ATXY];
    hAtxz[pp] = upt[U_ATXZ];
    hAtyy[pp] = upt[U_ATYY];
    hAtyz[pp] = upt[U_ATYZ];
    hAtzz[pp] = upt[U_ATZZ];

  }
  }
  }
}

/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
