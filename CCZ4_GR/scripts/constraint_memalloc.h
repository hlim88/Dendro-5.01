  const unsigned int bytes = n * sizeof(double);
  double *grad_0_psi = (double *) malloc(bytes);
  double *grad_1_psi = (double *) malloc(bytes);
  double *grad_2_psi = (double *) malloc(bytes);
  double *grad_0_Gh0 = (double *) malloc(bytes);
  double *grad_1_Gh0 = (double *) malloc(bytes);
  double *grad_2_Gh0 = (double *) malloc(bytes);
  double *grad_0_Gh1 = (double *) malloc(bytes);
  double *grad_1_Gh1 = (double *) malloc(bytes);
  double *grad_2_Gh1 = (double *) malloc(bytes);
  double *grad_0_Gh2 = (double *) malloc(bytes);
  double *grad_1_Gh2 = (double *) malloc(bytes);
  double *grad_2_Gh2 = (double *) malloc(bytes);
  double *grad_0_K = (double *) malloc(bytes);
  double *grad_1_K = (double *) malloc(bytes);
  double *grad_2_K = (double *) malloc(bytes);
  double *grad_0_theta_z4 = (double *) malloc(bytes);
  double *grad_1_theta_z4 = (double *) malloc(bytes);
  double *grad_2_theta_z4 = (double *) malloc(bytes);
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
  double *grad2_0_0_psi = (double *) malloc(bytes);
  double *grad2_0_1_psi = (double *) malloc(bytes);
  double *grad2_0_2_psi = (double *) malloc(bytes);
  double *grad2_1_1_psi = (double *) malloc(bytes);
  double *grad2_1_2_psi = (double *) malloc(bytes);
  double *grad2_2_2_psi = (double *) malloc(bytes);
