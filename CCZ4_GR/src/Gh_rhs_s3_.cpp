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
// Dendro: original ops:  18
// Dendro: printing temp variables
double DENDRO_0 = grad_0_beta0[pp] + grad_1_beta1[pp] + grad_2_beta2[pp];
// Dendro: printing variables
//--
Gh_rhs_s3_0[pp] = DENDRO_0*Gh0[pp];
//--
Gh_rhs_s3_1[pp] = DENDRO_0*Gh1[pp];
//--
Gh_rhs_s3_2[pp] = DENDRO_0*Gh2[pp];
// Dendro: reduced ops:  8
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
