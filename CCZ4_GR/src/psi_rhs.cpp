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
// Dendro: original ops:  20
// Dendro: printing temp variables
// Dendro: printing variables
//--
psi_rhs[pp] = beta0[pp]*agrad_0_psi[pp] + beta1[pp]*agrad_1_psi[pp] + beta2[pp]*agrad_2_psi[pp] + (2.0L/3.0L)*psi[pp]*(-K[pp]*alpha[pp] + grad_0_beta0[pp] + grad_1_beta1[pp] + grad_2_beta2[pp])/p_expo;
// Dendro: reduced ops:  20
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
