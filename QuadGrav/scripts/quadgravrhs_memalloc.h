  const unsigned int bytes = n * sizeof(double);
  double *grad_0_alpha = (double *) malloc(bytes);
  double *grad_1_alpha = (double *) malloc(bytes);
  double *grad_2_alpha = (double *) malloc(bytes);
  double *grad_0_beta0 = (double *) malloc(bytes);
  double *grad_1_beta0 = (double *) malloc(bytes);
  double *grad_2_beta0 = (double *) malloc(bytes);
  double *grad_0_beta1 = (double *) malloc(bytes);
  double *grad_1_beta1 = (double *) malloc(bytes);
  double *grad_2_beta1 = (double *) malloc(bytes);
  double *grad_0_beta2 = (double *) malloc(bytes);
  double *grad_1_beta2 = (double *) malloc(bytes);
  double *grad_2_beta2 = (double *) malloc(bytes);
  double *grad_0_B0 = (double *) malloc(bytes);
  double *grad_1_B0 = (double *) malloc(bytes);
  double *grad_2_B0 = (double *) malloc(bytes);
  double *grad_0_B1 = (double *) malloc(bytes);
  double *grad_1_B1 = (double *) malloc(bytes);
  double *grad_2_B1 = (double *) malloc(bytes);
  double *grad_0_B2 = (double *) malloc(bytes);
  double *grad_1_B2 = (double *) malloc(bytes);
  double *grad_2_B2 = (double *) malloc(bytes);
  double *grad_0_chi = (double *) malloc(bytes);
  double *grad_1_chi = (double *) malloc(bytes);
  double *grad_2_chi = (double *) malloc(bytes);
  double *grad_0_Gt0 = (double *) malloc(bytes);
  double *grad_1_Gt0 = (double *) malloc(bytes);
  double *grad_2_Gt0 = (double *) malloc(bytes);
  double *grad_0_Gt1 = (double *) malloc(bytes);
  double *grad_1_Gt1 = (double *) malloc(bytes);
  double *grad_2_Gt1 = (double *) malloc(bytes);
  double *grad_0_Gt2 = (double *) malloc(bytes);
  double *grad_1_Gt2 = (double *) malloc(bytes);
  double *grad_2_Gt2 = (double *) malloc(bytes);
  double *grad_0_K = (double *) malloc(bytes);
  double *grad_1_K = (double *) malloc(bytes);
  double *grad_2_K = (double *) malloc(bytes);
  double *grad_0_gt0 = (double *) malloc(bytes);
  double *grad_1_gt0 = (double *) malloc(bytes);
  double *grad_2_gt0 = (double *) malloc(bytes);
  double *grad_0_gt1 = (double *) malloc(bytes);
  double *grad_1_gt1 = (double *) malloc(bytes);
  double *grad_2_gt1 = (double *) malloc(bytes);
  double *grad_0_gt2 = (double *) malloc(bytes);
  double *grad_1_gt2 = (double *) malloc(bytes);
  double *grad_2_gt2 = (double *) malloc(bytes);
  double *grad_0_gt3 = (double *) malloc(bytes);
  double *grad_1_gt3 = (double *) malloc(bytes);
  double *grad_2_gt3 = (double *) malloc(bytes);
  double *grad_0_gt4 = (double *) malloc(bytes);
  double *grad_1_gt4 = (double *) malloc(bytes);
  double *grad_2_gt4 = (double *) malloc(bytes);
  double *grad_0_gt5 = (double *) malloc(bytes);
  double *grad_1_gt5 = (double *) malloc(bytes);
  double *grad_2_gt5 = (double *) malloc(bytes);
  double *grad_0_At0 = (double *) malloc(bytes);
  double *grad_1_At0 = (double *) malloc(bytes);
  double *grad_2_At0 = (double *) malloc(bytes);
  double *grad_0_At1 = (double *) malloc(bytes);
  double *grad_1_At1 = (double *) malloc(bytes);
  double *grad_2_At1 = (double *) malloc(bytes);
  double *grad_0_At2 = (double *) malloc(bytes);
  double *grad_1_At2 = (double *) malloc(bytes);
  double *grad_2_At2 = (double *) malloc(bytes);
  double *grad_0_At3 = (double *) malloc(bytes);
  double *grad_1_At3 = (double *) malloc(bytes);
  double *grad_2_At3 = (double *) malloc(bytes);
  double *grad_0_At4 = (double *) malloc(bytes);
  double *grad_1_At4 = (double *) malloc(bytes);
  double *grad_2_At4 = (double *) malloc(bytes);
  double *grad_0_At5 = (double *) malloc(bytes);
  double *grad_1_At5 = (double *) malloc(bytes);
  double *grad_2_At5 = (double *) malloc(bytes);
  double *grad_0_Rsc = (double *) malloc(bytes);
  double *grad_1_Rsc = (double *) malloc(bytes);
  double *grad_2_Rsc = (double *) malloc(bytes);
  double *grad_0_Rsch = (double *) malloc(bytes);
  double *grad_1_Rsch = (double *) malloc(bytes);
  double *grad_2_Rsch = (double *) malloc(bytes);
  double *grad_0_Atr = (double *) malloc(bytes);
  double *grad_1_Atr = (double *) malloc(bytes);
  double *grad_2_Atr = (double *) malloc(bytes);
  double *grad_0_Aij0 = (double *) malloc(bytes);
  double *grad_1_Aij0 = (double *) malloc(bytes);
  double *grad_2_Aij0 = (double *) malloc(bytes);
  double *grad_0_Aij1 = (double *) malloc(bytes);
  double *grad_1_Aij1 = (double *) malloc(bytes);
  double *grad_2_Aij1 = (double *) malloc(bytes);
  double *grad_0_Aij2 = (double *) malloc(bytes);
  double *grad_1_Aij2 = (double *) malloc(bytes);
  double *grad_2_Aij2 = (double *) malloc(bytes);
  double *grad_0_Aij3 = (double *) malloc(bytes);
  double *grad_1_Aij3 = (double *) malloc(bytes);
  double *grad_2_Aij3 = (double *) malloc(bytes);
  double *grad_0_Aij4 = (double *) malloc(bytes);
  double *grad_1_Aij4 = (double *) malloc(bytes);
  double *grad_2_Aij4 = (double *) malloc(bytes);
  double *grad_0_Aij5 = (double *) malloc(bytes);
  double *grad_1_Aij5 = (double *) malloc(bytes);
  double *grad_2_Aij5 = (double *) malloc(bytes);
  double *grad_0_Btr = (double *) malloc(bytes);
  double *grad_1_Btr = (double *) malloc(bytes);
  double *grad_2_Btr = (double *) malloc(bytes);
  double *grad_0_Bij0 = (double *) malloc(bytes);
  double *grad_1_Bij0 = (double *) malloc(bytes);
  double *grad_2_Bij0 = (double *) malloc(bytes);
  double *grad_0_Bij1 = (double *) malloc(bytes);
  double *grad_1_Bij1 = (double *) malloc(bytes);
  double *grad_2_Bij1 = (double *) malloc(bytes);
  double *grad_0_Bij2 = (double *) malloc(bytes);
  double *grad_1_Bij2 = (double *) malloc(bytes);
  double *grad_2_Bij2 = (double *) malloc(bytes);
  double *grad_0_Bij3 = (double *) malloc(bytes);
  double *grad_1_Bij3 = (double *) malloc(bytes);
  double *grad_2_Bij3 = (double *) malloc(bytes);
  double *grad_0_Bij4 = (double *) malloc(bytes);
  double *grad_1_Bij4 = (double *) malloc(bytes);
  double *grad_2_Bij4 = (double *) malloc(bytes);
  double *grad_0_Bij5 = (double *) malloc(bytes);
  double *grad_1_Bij5 = (double *) malloc(bytes);
  double *grad_2_Bij5 = (double *) malloc(bytes);
  double *grad_0_Ci0 = (double *) malloc(bytes);
  double *grad_1_Ci0 = (double *) malloc(bytes);
  double *grad_2_Ci0 = (double *) malloc(bytes);
  double *grad_0_Ci1 = (double *) malloc(bytes);
  double *grad_1_Ci1 = (double *) malloc(bytes);
  double *grad_2_Ci1 = (double *) malloc(bytes);
  double *grad_0_Ci2 = (double *) malloc(bytes);
  double *grad_1_Ci2 = (double *) malloc(bytes);
  double *grad_2_Ci2 = (double *) malloc(bytes);
  double *grad2_0_0_gt0 = (double *) malloc(bytes);
  double *grad2_0_1_gt0 = (double *) malloc(bytes);
  double *grad2_0_2_gt0 = (double *) malloc(bytes);
  double *grad2_1_1_gt0 = (double *) malloc(bytes);
  double *grad2_1_2_gt0 = (double *) malloc(bytes);
  double *grad2_2_2_gt0 = (double *) malloc(bytes);
  double *grad2_0_0_gt1 = (double *) malloc(bytes);
  double *grad2_0_1_gt1 = (double *) malloc(bytes);
  double *grad2_0_2_gt1 = (double *) malloc(bytes);
  double *grad2_1_1_gt1 = (double *) malloc(bytes);
  double *grad2_1_2_gt1 = (double *) malloc(bytes);
  double *grad2_2_2_gt1 = (double *) malloc(bytes);
  double *grad2_0_0_gt2 = (double *) malloc(bytes);
  double *grad2_0_1_gt2 = (double *) malloc(bytes);
  double *grad2_0_2_gt2 = (double *) malloc(bytes);
  double *grad2_1_1_gt2 = (double *) malloc(bytes);
  double *grad2_1_2_gt2 = (double *) malloc(bytes);
  double *grad2_2_2_gt2 = (double *) malloc(bytes);
  double *grad2_0_0_gt3 = (double *) malloc(bytes);
  double *grad2_0_1_gt3 = (double *) malloc(bytes);
  double *grad2_0_2_gt3 = (double *) malloc(bytes);
  double *grad2_1_1_gt3 = (double *) malloc(bytes);
  double *grad2_1_2_gt3 = (double *) malloc(bytes);
  double *grad2_2_2_gt3 = (double *) malloc(bytes);
  double *grad2_0_0_gt4 = (double *) malloc(bytes);
  double *grad2_0_1_gt4 = (double *) malloc(bytes);
  double *grad2_0_2_gt4 = (double *) malloc(bytes);
  double *grad2_1_1_gt4 = (double *) malloc(bytes);
  double *grad2_1_2_gt4 = (double *) malloc(bytes);
  double *grad2_2_2_gt4 = (double *) malloc(bytes);
  double *grad2_0_0_gt5 = (double *) malloc(bytes);
  double *grad2_0_1_gt5 = (double *) malloc(bytes);
  double *grad2_0_2_gt5 = (double *) malloc(bytes);
  double *grad2_1_1_gt5 = (double *) malloc(bytes);
  double *grad2_1_2_gt5 = (double *) malloc(bytes);
  double *grad2_2_2_gt5 = (double *) malloc(bytes);
  double *grad2_0_0_chi = (double *) malloc(bytes);
  double *grad2_0_1_chi = (double *) malloc(bytes);
  double *grad2_0_2_chi = (double *) malloc(bytes);
  double *grad2_1_1_chi = (double *) malloc(bytes);
  double *grad2_1_2_chi = (double *) malloc(bytes);
  double *grad2_2_2_chi = (double *) malloc(bytes);
  double *grad2_0_0_alpha = (double *) malloc(bytes);
  double *grad2_0_1_alpha = (double *) malloc(bytes);
  double *grad2_0_2_alpha = (double *) malloc(bytes);
  double *grad2_1_1_alpha = (double *) malloc(bytes);
  double *grad2_1_2_alpha = (double *) malloc(bytes);
  double *grad2_2_2_alpha = (double *) malloc(bytes);
  double *grad2_0_0_beta0 = (double *) malloc(bytes);
  double *grad2_0_1_beta0 = (double *) malloc(bytes);
  double *grad2_0_2_beta0 = (double *) malloc(bytes);
  double *grad2_1_1_beta0 = (double *) malloc(bytes);
  double *grad2_1_2_beta0 = (double *) malloc(bytes);
  double *grad2_2_2_beta0 = (double *) malloc(bytes);
  double *grad2_0_0_beta1 = (double *) malloc(bytes);
  double *grad2_0_1_beta1 = (double *) malloc(bytes);
  double *grad2_0_2_beta1 = (double *) malloc(bytes);
  double *grad2_1_1_beta1 = (double *) malloc(bytes);
  double *grad2_1_2_beta1 = (double *) malloc(bytes);
  double *grad2_2_2_beta1 = (double *) malloc(bytes);
  double *grad2_0_0_beta2 = (double *) malloc(bytes);
  double *grad2_0_1_beta2 = (double *) malloc(bytes);
  double *grad2_0_2_beta2 = (double *) malloc(bytes);
  double *grad2_1_1_beta2 = (double *) malloc(bytes);
  double *grad2_1_2_beta2 = (double *) malloc(bytes);
  double *grad2_2_2_beta2 = (double *) malloc(bytes);
  double *grad2_0_0_Rsc = (double *) malloc(bytes);
  double *grad2_0_1_Rsc = (double *) malloc(bytes);
  double *grad2_0_2_Rsc = (double *) malloc(bytes);
  double *grad2_1_1_Rsc = (double *) malloc(bytes);
  double *grad2_1_2_Rsc = (double *) malloc(bytes);
  double *grad2_2_2_Rsc = (double *) malloc(bytes);
  double *grad2_0_0_Rsch = (double *) malloc(bytes);
  double *grad2_0_1_Rsch = (double *) malloc(bytes);
  double *grad2_0_2_Rsch = (double *) malloc(bytes);
  double *grad2_1_1_Rsch = (double *) malloc(bytes);
  double *grad2_1_2_Rsch = (double *) malloc(bytes);
  double *grad2_2_2_Rsch = (double *) malloc(bytes);
  double *grad2_0_0_Atr = (double *) malloc(bytes);
  double *grad2_0_1_Atr = (double *) malloc(bytes);
  double *grad2_0_2_Atr = (double *) malloc(bytes);
  double *grad2_1_1_Atr = (double *) malloc(bytes);
  double *grad2_1_2_Atr = (double *) malloc(bytes);
  double *grad2_2_2_Atr = (double *) malloc(bytes);
  double *grad2_0_0_Aij0 = (double *) malloc(bytes);
  double *grad2_0_1_Aij0 = (double *) malloc(bytes);
  double *grad2_0_2_Aij0 = (double *) malloc(bytes);
  double *grad2_1_1_Aij0 = (double *) malloc(bytes);
  double *grad2_1_2_Aij0 = (double *) malloc(bytes);
  double *grad2_2_2_Aij0 = (double *) malloc(bytes);
  double *grad2_0_0_Aij1 = (double *) malloc(bytes);
  double *grad2_0_1_Aij1 = (double *) malloc(bytes);
  double *grad2_0_2_Aij1 = (double *) malloc(bytes);
  double *grad2_1_1_Aij1 = (double *) malloc(bytes);
  double *grad2_1_2_Aij1 = (double *) malloc(bytes);
  double *grad2_2_2_Aij1 = (double *) malloc(bytes);
  double *grad2_0_0_Aij2 = (double *) malloc(bytes);
  double *grad2_0_1_Aij2 = (double *) malloc(bytes);
  double *grad2_0_2_Aij2 = (double *) malloc(bytes);
  double *grad2_1_1_Aij2 = (double *) malloc(bytes);
  double *grad2_1_2_Aij2 = (double *) malloc(bytes);
  double *grad2_2_2_Aij2 = (double *) malloc(bytes);
  double *grad2_0_0_Aij3 = (double *) malloc(bytes);
  double *grad2_0_1_Aij3 = (double *) malloc(bytes);
  double *grad2_0_2_Aij3 = (double *) malloc(bytes);
  double *grad2_1_1_Aij3 = (double *) malloc(bytes);
  double *grad2_1_2_Aij3 = (double *) malloc(bytes);
  double *grad2_2_2_Aij3 = (double *) malloc(bytes);
  double *grad2_0_0_Aij4 = (double *) malloc(bytes);
  double *grad2_0_1_Aij4 = (double *) malloc(bytes);
  double *grad2_0_2_Aij4 = (double *) malloc(bytes);
  double *grad2_1_1_Aij4 = (double *) malloc(bytes);
  double *grad2_1_2_Aij4 = (double *) malloc(bytes);
  double *grad2_2_2_Aij4 = (double *) malloc(bytes);
  double *grad2_0_0_Aij5 = (double *) malloc(bytes);
  double *grad2_0_1_Aij5 = (double *) malloc(bytes);
  double *grad2_0_2_Aij5 = (double *) malloc(bytes);
  double *grad2_1_1_Aij5 = (double *) malloc(bytes);
  double *grad2_1_2_Aij5 = (double *) malloc(bytes);
  double *grad2_2_2_Aij5 = (double *) malloc(bytes);
  double *grad2_0_0_Btr = (double *) malloc(bytes);
  double *grad2_0_1_Btr = (double *) malloc(bytes);
  double *grad2_0_2_Btr = (double *) malloc(bytes);
  double *grad2_1_1_Btr = (double *) malloc(bytes);
  double *grad2_1_2_Btr = (double *) malloc(bytes);
  double *grad2_2_2_Btr = (double *) malloc(bytes);
  double *grad2_0_0_Bij0 = (double *) malloc(bytes);
  double *grad2_0_1_Bij0 = (double *) malloc(bytes);
  double *grad2_0_2_Bij0 = (double *) malloc(bytes);
  double *grad2_1_1_Bij0 = (double *) malloc(bytes);
  double *grad2_1_2_Bij0 = (double *) malloc(bytes);
  double *grad2_2_2_Bij0 = (double *) malloc(bytes);
  double *grad2_0_0_Bij1 = (double *) malloc(bytes);
  double *grad2_0_1_Bij1 = (double *) malloc(bytes);
  double *grad2_0_2_Bij1 = (double *) malloc(bytes);
  double *grad2_1_1_Bij1 = (double *) malloc(bytes);
  double *grad2_1_2_Bij1 = (double *) malloc(bytes);
  double *grad2_2_2_Bij1 = (double *) malloc(bytes);
  double *grad2_0_0_Bij2 = (double *) malloc(bytes);
  double *grad2_0_1_Bij2 = (double *) malloc(bytes);
  double *grad2_0_2_Bij2 = (double *) malloc(bytes);
  double *grad2_1_1_Bij2 = (double *) malloc(bytes);
  double *grad2_1_2_Bij2 = (double *) malloc(bytes);
  double *grad2_2_2_Bij2 = (double *) malloc(bytes);
  double *grad2_0_0_Bij3 = (double *) malloc(bytes);
  double *grad2_0_1_Bij3 = (double *) malloc(bytes);
  double *grad2_0_2_Bij3 = (double *) malloc(bytes);
  double *grad2_1_1_Bij3 = (double *) malloc(bytes);
  double *grad2_1_2_Bij3 = (double *) malloc(bytes);
  double *grad2_2_2_Bij3 = (double *) malloc(bytes);
  double *grad2_0_0_Bij4 = (double *) malloc(bytes);
  double *grad2_0_1_Bij4 = (double *) malloc(bytes);
  double *grad2_0_2_Bij4 = (double *) malloc(bytes);
  double *grad2_1_1_Bij4 = (double *) malloc(bytes);
  double *grad2_1_2_Bij4 = (double *) malloc(bytes);
  double *grad2_2_2_Bij4 = (double *) malloc(bytes);
  double *grad2_0_0_Bij5 = (double *) malloc(bytes);
  double *grad2_0_1_Bij5 = (double *) malloc(bytes);
  double *grad2_0_2_Bij5 = (double *) malloc(bytes);
  double *grad2_1_1_Bij5 = (double *) malloc(bytes);
  double *grad2_1_2_Bij5 = (double *) malloc(bytes);
  double *grad2_2_2_Bij5 = (double *) malloc(bytes);
  double *grad2_0_0_Ci0 = (double *) malloc(bytes);
  double *grad2_0_1_Ci0 = (double *) malloc(bytes);
  double *grad2_0_2_Ci0 = (double *) malloc(bytes);
  double *grad2_1_1_Ci0 = (double *) malloc(bytes);
  double *grad2_1_2_Ci0 = (double *) malloc(bytes);
  double *grad2_2_2_Ci0 = (double *) malloc(bytes);
  double *grad2_0_0_Ci1 = (double *) malloc(bytes);
  double *grad2_0_1_Ci1 = (double *) malloc(bytes);
  double *grad2_0_2_Ci1 = (double *) malloc(bytes);
  double *grad2_1_1_Ci1 = (double *) malloc(bytes);
  double *grad2_1_2_Ci1 = (double *) malloc(bytes);
  double *grad2_2_2_Ci1 = (double *) malloc(bytes);
  double *grad2_0_0_Ci2 = (double *) malloc(bytes);
  double *grad2_0_1_Ci2 = (double *) malloc(bytes);
  double *grad2_0_2_Ci2 = (double *) malloc(bytes);
  double *grad2_1_1_Ci2 = (double *) malloc(bytes);
  double *grad2_1_2_Ci2 = (double *) malloc(bytes);
  double *grad2_2_2_Ci2 = (double *) malloc(bytes);
