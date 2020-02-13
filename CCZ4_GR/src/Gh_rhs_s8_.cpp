    ccz4::timer::t_rhs.start();
for (unsigned int k = 3; k < nz-3; k++) { 
    z = pmin[2] + k*hz;
for (unsigned int j = 3; j < ny-3; j++) { 
    y = pmin[1] + j*hy; 
for (unsigned int i = 3; i < nx-3; i++) {
    x = pmin[0] + i*hx;
    pp = i + nx*(j + ny*k);
    r_coord = sqrt(x*x + y*y + z*z);
    eta=ETA_CONST;
    if (r_coord >= ETA_R0) {
    eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
    }
// Dendro: {{{ 
// Dendro: original ops:  1980
// Dendro: printing temp variables
double DENDRO_0 = gt1[pp]*gt5[pp] - gt2[pp]*gt4[pp];
double DENDRO_1 = pow(gt4[pp], 2);
double DENDRO_2 = pow(gt1[pp], 2);
double DENDRO_3 = pow(gt2[pp], 2);
double DENDRO_4 = gt0[pp]*gt3[pp];
double DENDRO_5 = gt1[pp]*gt2[pp];
double DENDRO_6 = DENDRO_1*gt0[pp] + DENDRO_2*gt5[pp] + DENDRO_3*gt3[pp] - DENDRO_4*gt5[pp] - 2*DENDRO_5*gt4[pp];
double DENDRO_7 = 1.0/DENDRO_6;
double DENDRO_8 = 2*DENDRO_7*(alpha[pp]*grad_1_theta_z4[pp] - theta_z4[pp]*grad_1_alpha[pp]);
double DENDRO_9 = gt1[pp]*gt4[pp] - gt2[pp]*gt3[pp];
double DENDRO_10 = 2*DENDRO_7*(alpha[pp]*grad_2_theta_z4[pp] - theta_z4[pp]*grad_2_alpha[pp]);
double DENDRO_11 = -DENDRO_1 + gt3[pp]*gt5[pp];
double DENDRO_12 = 2*DENDRO_7*(alpha[pp]*grad_0_theta_z4[pp] - theta_z4[pp]*grad_0_alpha[pp]);
double DENDRO_13 = 0.333333333333333*alpha[pp]*(2*K[pp] + 3*kappa[0]);
double DENDRO_14 = -DENDRO_5 + gt0[pp]*gt4[pp];
double DENDRO_15 = pow(DENDRO_6, -2);
double DENDRO_16 = DENDRO_14*DENDRO_15;
double DENDRO_17 = grad_2_gt3[pp];
double DENDRO_18 = grad_1_gt5[pp];
double DENDRO_19 = grad_1_gt2[pp];
double DENDRO_20 = grad_2_gt1[pp];
double DENDRO_21 = grad_0_gt4[pp];
double DENDRO_22 = DENDRO_19 + DENDRO_20 - DENDRO_21;
double DENDRO_23 = DENDRO_0*DENDRO_15;
double DENDRO_24 = grad_0_gt3[pp];
double DENDRO_25 = grad_1_gt0[pp];
double DENDRO_26 = DENDRO_19 - DENDRO_20 + DENDRO_21;
double DENDRO_27 = 1.0*DENDRO_15*DENDRO_9;
double DENDRO_28 = grad_0_gt5[pp];
double DENDRO_29 = grad_2_gt0[pp];
double DENDRO_30 = -DENDRO_19 + DENDRO_20 + DENDRO_21;
double DENDRO_31 = -DENDRO_2 + DENDRO_4;
double DENDRO_32 = DENDRO_15*DENDRO_31;
double DENDRO_33 = grad_2_gt5[pp];
double DENDRO_34 = 0.5*gt1[pp]*gt4[pp] - 0.5*gt2[pp]*gt3[pp];
double DENDRO_35 = 0.5*DENDRO_18 - 1.0*grad_2_gt4[pp];
double DENDRO_36 = 0.5*DENDRO_28 - 1.0*grad_2_gt2[pp];
double DENDRO_37 = -DENDRO_3 + gt0[pp]*gt5[pp];
double DENDRO_38 = DENDRO_15*DENDRO_37;
double DENDRO_39 = grad_1_gt3[pp];
double DENDRO_40 = 0.5*gt1[pp]*gt5[pp] - 0.5*gt2[pp]*gt4[pp];
double DENDRO_41 = -0.5*DENDRO_17 + 1.0*grad_1_gt4[pp];
double DENDRO_42 = 0.5*DENDRO_24 - 1.0*grad_1_gt1[pp];
double DENDRO_43 = DENDRO_11*DENDRO_15;
double DENDRO_44 = grad_0_gt0[pp];
double DENDRO_45 = -0.5*DENDRO_25 + 1.0*grad_0_gt1[pp];
double DENDRO_46 = -0.5*DENDRO_29 + 1.0*grad_0_gt2[pp];
double DENDRO_47 = 0.5*gt0[pp]*gt4[pp] - 0.5*gt1[pp]*gt2[pp];
// Dendro: printing variables
//--
Gh_rhs_s8_0[pp] = DENDRO_0*DENDRO_8 - DENDRO_10*DENDRO_9 - DENDRO_11*DENDRO_12 + DENDRO_13*(DENDRO_16*(DENDRO_0*DENDRO_17 - DENDRO_11*DENDRO_22 - DENDRO_18*DENDRO_9) + DENDRO_23*(DENDRO_0*DENDRO_24 - DENDRO_11*DENDRO_25 - DENDRO_26*DENDRO_9) - DENDRO_27*(DENDRO_0*DENDRO_30 - DENDRO_11*DENDRO_29 - DENDRO_28*DENDRO_9) - DENDRO_32*(-DENDRO_0*DENDRO_35 + DENDRO_11*DENDRO_36 - DENDRO_33*DENDRO_34) - DENDRO_38*(DENDRO_11*DENDRO_42 + DENDRO_39*DENDRO_40 - DENDRO_41*DENDRO_9) - DENDRO_43*(DENDRO_0*DENDRO_45 - 0.5*DENDRO_11*DENDRO_44 - DENDRO_46*DENDRO_9) - Gh0[pp]);
//--
Gh_rhs_s8_1[pp] = DENDRO_0*DENDRO_12 + DENDRO_10*DENDRO_14 + DENDRO_13*(DENDRO_16*(DENDRO_0*DENDRO_22 + DENDRO_14*DENDRO_18 - DENDRO_17*DENDRO_37) + DENDRO_23*(DENDRO_0*DENDRO_25 + DENDRO_14*DENDRO_26 - DENDRO_24*DENDRO_37) - DENDRO_27*(DENDRO_0*DENDRO_29 + DENDRO_14*DENDRO_28 - DENDRO_30*DENDRO_37) - DENDRO_32*(-DENDRO_0*DENDRO_36 + DENDRO_33*DENDRO_47 + DENDRO_35*DENDRO_37) - DENDRO_38*(-DENDRO_0*DENDRO_42 + DENDRO_14*DENDRO_41 - 0.5*DENDRO_37*DENDRO_39) - DENDRO_43*(DENDRO_14*DENDRO_46 - DENDRO_37*DENDRO_45 + DENDRO_40*DENDRO_44) - Gh1[pp]) - DENDRO_37*DENDRO_8;
//--
Gh_rhs_s8_2[pp] = -DENDRO_10*DENDRO_31 - DENDRO_12*DENDRO_9 + DENDRO_13*(DENDRO_16*(DENDRO_14*DENDRO_17 - DENDRO_18*DENDRO_31 - DENDRO_22*DENDRO_9) + DENDRO_23*(DENDRO_14*DENDRO_24 - DENDRO_25*DENDRO_9 - DENDRO_26*DENDRO_31) - DENDRO_27*(DENDRO_14*DENDRO_30 - DENDRO_28*DENDRO_31 - DENDRO_29*DENDRO_9) - DENDRO_32*(-DENDRO_14*DENDRO_35 - 0.5*DENDRO_31*DENDRO_33 + DENDRO_36*DENDRO_9) - DENDRO_38*(-DENDRO_31*DENDRO_41 + DENDRO_39*DENDRO_47 + DENDRO_42*DENDRO_9) - DENDRO_43*(DENDRO_14*DENDRO_45 - DENDRO_31*DENDRO_46 - DENDRO_34*DENDRO_44) - Gh2[pp]) + DENDRO_14*DENDRO_8;
// Dendro: reduced ops:  271
// Dendro: }}} 
     /* debugging */
     /*unsigned int qi = 46 - 1;
     unsigned int qj = 10 - 1;
     unsigned int qk = 60 - 1;
     unsigned int qidx = qi + nx*(qj + ny*qk);
     if (0 && qidx == pp) {
     std::cout << ".... end OPTIMIZED debug stuff..." << std::endl;
     }*/
  }
 }
}
     ccz4::timer::t_rhs.stop();
