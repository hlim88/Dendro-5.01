#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "cctk.h"
#include "params.h"
// #include <assert.h>

/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void init_data_puncture(const double *coord, double *u) 
{
  double epijk[3][3][3];
  int i,j,k;
  for (k=0;k<3;k++) {
    for (j=0;j<3;j++) {
      for (i=0;i<3;i++) {
        epijk[k][j][i] = 0.0;
      }
    }
  }
  epijk[0][1][2] = 1.0;epijk[1][2][0] = 1.0;epijk[2][0][1] = 1.0;
  epijk[0][2][1] = -1.0;epijk[2][1][0] = -1.0;epijk[1][0][2] = -1.0;

  double deltaij[3][3];
  for (j=0;j<3;j++) {
    for (i=0;i<3;i++) {
      deltaij[j][i] = 0.0;
    }
  }
  deltaij[0][0] = 1.0;deltaij[1][1] = 1.0;deltaij[2][2] = 1.0;

  double mass1 = par_m_plus;
  double bh1x = par_b;
  double bh1y =  0.0;
  double bh1z =  0.0;
  double vp1[3];
  vp1[0] = par_P_plus[0];
  vp1[1] = par_P_plus[1];
  vp1[2] = par_P_plus[2];
  double vp1tot = sqrt( vp1[0]*vp1[0] + vp1[1]*vp1[1] + vp1[2]*vp1[2] );
  double vs1[3];
  vs1[0] = par_S_plus[0];
  vs1[1] = par_S_plus[1];
  vs1[2] = par_S_plus[2];
  double spin1 = sqrt( vs1[0]*vs1[0] + vs1[1]*vs1[1] + vs1[2]*vs1[2] );

  // bh 2
  double mass2 = par_m_minus;
  double bh2x =  -par_b;
  double bh2y =  0.0;
  double bh2z =  0.0;
  double vp2[3];
  vp2[0] = par_P_minus[0];
  vp2[1] = par_P_minus[1];
  vp2[2] = par_P_minus[2];
  double vp2tot = sqrt( vp2[0]*vp2[0] + vp2[1]*vp2[1] + vp2[2]*vp2[2] );
  double vs2[3];
  vs2[0] = par_S_minus[0];
  vs2[1] = par_S_minus[1];
  vs2[2] = par_S_minus[2];
  double spin2 = sqrt( vs2[0]*vs2[0] + vs2[1]*vs2[1] + vs2[2]*vs2[2] );

  double xx = coord[0];
  double yy = coord[1];
  double zz = coord[2];

  double x1,y1,z1,rv1;
  double x2,y2,z2,rv2;
  double vn1[3],vn2[3];
  double vpsibl;
  double v_u_corr,amp_capj,amp_capr,l_r,u0_j,u2_j,mu_j,p2_mu_j,v_u_j1;
  double v1,v2,v3,v4,vt1,vt2;
  int i1,i2,i3,i4;
  double amp_capp,u0_p,u2_p,mu_p,p2_mu_p;
  double v_u_p1,v_u_c1,v_u_j2,v_u_p2;
  double v_u_c2,vpsibl_u,vpsibl_u2;

  //allocating the location of BH1

  x1 = xx - bh1x;
  y1 = yy - bh1y;
  z1 = zz - bh1z;

  //locating as a radial form
  rv1 = sqrt(x1*x1 + y1*y1 + z1*z1);
  vn1[0] = x1/rv1;
  vn1[1] = y1/rv1;
  vn1[2] = z1/rv1;

  //same as BH2
  x2 = xx - bh2x;
  y2 = yy - bh2y;
  z2 = zz - bh2z;

  rv2 = sqrt(x2*x2 + y2*y2 + z2*z2);
  vn2[0] = x2/rv2;
  vn2[1] = y2/rv2;
  vn2[2] = z2/rv2;

  //Initial data is related with the paper: http://arxiv.org/abs/0711.1165
  //Brill-Lindquist conformal factor
  vpsibl = 1.0 + mass1/(2.0*rv1);
  vpsibl = vpsibl + mass2/(2.0*rv2);

  v_u_corr = 0.0;
  // bh 1
  //For spinning puncture
  if ( fabs(spin1) > 1.e-6 ) {
      amp_capj = 4.0*spin1/(mass1*mass1);
      amp_capr = 2.0*rv1/mass1;
      l_r = 1.0/(1.0+amp_capr);
      u0_j = (l_r + l_r*l_r + l_r*l_r*l_r - 4.0*l_r*l_r*l_r*l_r + 2.0*l_r*l_r*l_r*l_r*l_r)/40.0;
      u2_j = -pow(l_r,5)/20.0;
      mu_j = vn1[0]*vs1[0];
      mu_j = mu_j + vn1[1]*vs1[1];
      mu_j = (mu_j + vn1[2]*vs1[2])/fabs(spin1);
      p2_mu_j = (3.0*mu_j*mu_j - 1.0)/2.0;
      v_u_j1 = amp_capj*amp_capj*(u0_j+u2_j*amp_capr*amp_capr*p2_mu_j);
      v_u_corr = v_u_corr + v_u_j1;
  }
  //For boosting puncture
  if (vp1tot > 1.e-6) {
      amp_capp = 2.0*vp1tot/mass1;
      amp_capr = 2.0*rv1/mass1;
      l_r = 1.0/(1.0 + amp_capr);
      u0_p = l_r - 2.0*l_r*l_r + 2.0*pow(l_r,3);
      u0_p = (u0_p - pow(l_r,4) + 0.20*pow(l_r,5))*(5.0/32.0);
      u2_p = 15.0*l_r + 132.0*l_r*l_r + 53.0*pow(l_r,3);
      u2_p = u2_p + 96.0*pow(l_r,4) + 82.0*pow(l_r,5);
      u2_p = u2_p + (84.0/amp_capr)*(pow(l_r,5)+log(l_r)/amp_capr);
      u2_p = (u2_p)/(80.0*amp_capr);
      mu_p =        vn1[0]*vp1[0]/vp1tot;
      mu_p = mu_p + vn1[1]*vp1[1]/vp1tot;
      mu_p = mu_p + vn1[2]*vp1[2]/vp1tot;
      p2_mu_p = (3.0*pow(mu_p,2) - 1.0)/2.0;
      v_u_p1 = pow(amp_capp,2)*(u0_p+u2_p*p2_mu_p);
      v_u_corr = v_u_corr + v_u_p1;
  }
  //For spinning boosted pucture
  if ( vp1tot > 1.e-6 && fabs(spin1) > 1.e-6 ) {
        v1 =      (vp1[1]*vs1[2]-vp1[2]*vs1[1])*vn1[0];
        v1 = v1 + (vp1[2]*vs1[0]-vp1[0]*vs1[2])*vn1[1];
        v1 = v1 + (vp1[0]*vs1[1]-vp1[1]*vs1[0])*vn1[2];
        v1 = v1*(16.0/pow(mass1,4))*rv1;
        
        amp_capr = 2.0*rv1/mass1;
        l_r = 1.0/(1.0 + amp_capr);
        
        v2 = 1.0 + 5.0*amp_capr + 10.0*pow(amp_capr,2);
        
        v_u_c1 = (v1*v2*pow(l_r,5))/80.0;
        v_u_corr = v_u_corr + v_u_c1;
  }
  // bh 2 same puncture as bh 1
  if ( fabs(spin2) > 1.e-6 ) {
        amp_capj = 4.0*spin2/(mass2*mass2);
        amp_capr = 2.0*rv2/mass2;
        l_r = 1.0/(1.0+amp_capr);
        u0_j = (l_r + l_r*l_r + l_r*l_r*l_r - 4.0*l_r*l_r*l_r*l_r + 2.0*l_r*l_r*l_r*l_r*l_r)/40.0;
        u2_j = -pow(l_r,5)/20.0;
        mu_j = vn2[0]*vs2[0];
        mu_j = mu_j + vn2[1]*vs2[1];
        mu_j = (mu_j + vn2[2]*vs2[2])/fabs(spin2);
        p2_mu_j = (3.0*mu_j*mu_j - 1.0)/2.0;
        v_u_j2 = amp_capj*amp_capj*(u0_j+u2_j*amp_capr*amp_capr*p2_mu_j);
        v_u_corr = v_u_corr + v_u_j2;
  }
    
  if (vp2tot > 1.e-6) {
        amp_capp = 2.0*vp2tot/mass2;
        amp_capr = 2.0*rv2/mass2;
        l_r = 1.0/(1.0 + amp_capr);
        u0_p = l_r - 2.0*l_r*l_r + 2.0*pow(l_r,3);
        u0_p = (u0_p - pow(l_r,4) + 0.20*pow(l_r,5))*(5.0/32.0);
        u2_p = 15.0*l_r + 132.0*l_r*l_r + 53.0*pow(l_r,3);
        u2_p = u2_p + 96.0*pow(l_r,4) + 82.0*pow(l_r,5);
        u2_p = u2_p + (84.0/amp_capr)*(pow(l_r,5)+log(l_r)/amp_capr);
        u2_p = (u2_p)/(80.0*amp_capr);
        mu_p =        vn2[0]*vp2[0]/vp2tot;
        mu_p = mu_p + vn2[1]*vp2[1]/vp2tot;
        mu_p = mu_p + vn2[2]*vp2[2]/vp2tot;
        p2_mu_p = (3.0*pow(mu_p,2) - 1.0)/2.0;
        v_u_p2 = pow(amp_capp,2)*(u0_p+u2_p*p2_mu_p);
        v_u_corr = v_u_corr + v_u_p2;
  }
    
  if ( vp2tot > 1.e-6 && fabs(spin2) > 1.e-6 ) {
        v1 =      (vp2[1]*vs2[2]-vp2[2]*vs2[1])*vn2[0];
        v1 = v1 + (vp2[2]*vs2[0]-vp2[0]*vs2[2])*vn2[1];
        v1 = v1 + (vp2[0]*vs2[1]-vp2[1]*vs2[0])*vn2[2];
        v1 = v1*(16.0/pow(mass2,4))*rv2;
        
        amp_capr = 2.0*rv2/mass2;
        l_r = 1.0/(1.0 + amp_capr);
        
        v2 = 1.0 + 5.0*amp_capr + 10.0*pow(amp_capr,2);
        
        v_u_c2 = (v1*v2*pow(l_r,5))/80.0;
        v_u_corr = v_u_corr + v_u_c2;
  }
 
  // vpsibl_u will be used for the conformal factor,
  vpsibl_u  = vpsibl + v_u_corr;
  // vpsibl_u2 is for the Aij terms...
  // ! since the corrections are first order...
  // ! adding half of the correction seems to give the best results...
  // ! update - do a fit for spin = 0.6...
  vpsibl_u2 = vpsibl + v_u_corr;

  u[U_ALPHA] = 1.0/(vpsibl_u*vpsibl_u);
  v2 = 1.0/pow(vpsibl_u,4);
  u[U_CHI] = v2;
  u[U_TRK] = 0.0;

  u[U_SHIFTX] = 0.0;
  u[U_SHIFTY] = 0.0;
  u[U_SHIFTZ] = 0.0;

  u[U_GAMTX] = 0.0;
  u[U_GAMTY] = 0.0;
  u[U_GAMTZ] = 0.0;

  u[U_GBX] = 0.0;
  u[U_GBY] = 0.0;
  u[U_GBZ] = 0.0;

  u[U_GTXX] = 1.0;
  u[U_GTXY] = 0.0;
  u[U_GTXZ] = 0.0;
  u[U_GTYY] = 1.0;
  u[U_GTYZ] = 0.0;
  u[U_GTZZ] = 1.0;

  for (i1=0;i1<3;i1++) {
    for (i2=0;i2<3;i2++) {
      // first BH
      v2 = 0.0;
      for (i3=0;i3<3;i3++) {
        for (i4=0;i4<3;i4++) {
          vt1 = epijk[i1][i3][i4]*vs1[i3]*vn1[i4]*vn1[i2];
          vt2 = epijk[i2][i3][i4]*vs1[i3]*vn1[i4]*vn1[i1];
          v2 = v2 + vt1 + vt2;
        }
      }

      v3 = vp1[i1]*vn1[i2] + vp1[i2]*vn1[i1];
      vt1 = 0.0;
      for (i3=0;i3<3;i3++) {
        vt1 = vt1 + vp1[i3]*vn1[i3];
      }
      vt1 = vt1*(vn1[i1]*vn1[i2] - deltaij[i1][i2]);
      v3 = v3 + vt1;

      v1 = 3.0/(pow(vpsibl_u2,6)*pow(rv1,3));
      v4 = v1*(v2+(rv1/2.0)*v3);

      // second BH
      v2 = 0.0;
      for (i3=0;i3<3;i3++) {
        for (i4=0;i4<3;i4++) {
          vt1 = epijk[i1][i3][i4]*vs2[i3]*vn2[i4]*vn2[i2];
          vt2 = epijk[i2][i3][i4]*vs2[i3]*vn2[i4]*vn2[i1];
          v2 = v2 + vt1 + vt2;
        }
      }

      v3 = vp2[i1]*vn2[i2] + vp2[i2]*vn2[i1];
      vt1 = 0.0;
      for (i3=0;i3<3;i3++) {
        vt1 = vt1 + vp2[i3]*vn2[i3];
      }
      vt1 = vt1*(vn2[i1]*vn2[i2] - deltaij[i1][i2]);
      v3 = v3 + vt1;

      v1 = 3.0/(pow(vpsibl_u2,6)*pow(rv2,3));
      v4 = v4 + v1*(v2+(rv2/2.0)*v3);

      if ( i1 == 0 && i2 == 0 ) {
        u[U_ATXX] = v4;
      } else if ( i1 == 0 && i2 == 1 ) {
        u[U_ATXY] = v4;
      } else if ( i1 == 0 && i2 == 2 ) {
        u[U_ATXZ] = v4;
      } else if ( i1 == 1 && i2 == 1 ) {
        u[U_ATYY] = v4;
      } else if ( i1 == 1 && i2 == 2 ) {
        u[U_ATYZ] = v4;
      } else if ( i1 == 2 && i2 == 2 ) {
        u[U_ATZZ] = v4;
      }

    }
  }

}


