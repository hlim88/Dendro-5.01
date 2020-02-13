import dendro_ccz4
from sympy import *

###################################################################
# initialize
###################################################################

l1, l2, l3, l4, eta = symbols('lambda[0] lambda[1] lambda[2] lambda[3] eta')
lf0, lf1 = symbols('lambda_f[0] lambda_f[1]')
k1, k2 = symbols('kappa[0] kappa[1]')
p_expo = symbols('p_expo') 
#psi_to_neg_p_expo = symbols('psi_to_neg_p_expo') 

# declare variables
a   = dendro_ccz4.scalar("alpha", "[pp]")
psi = dendro_ccz4.scalar("psi", "[pp]")
K   = dendro_ccz4.scalar("K", "[pp]")

Gh  = dendro_ccz4.vec3("Gh", "[pp]")
b   = dendro_ccz4.vec3("beta", "[pp]")
B   = dendro_ccz4.vec3("B", "[pp]")

gt  = dendro_ccz4.sym_3x3("gt", "[pp]")
At  = dendro_ccz4.sym_3x3("At", "[pp]")

Gh_rhs  = dendro_ccz4.vec3("Gh_rhs", "[pp]")

# Additional variable for CCZ4
theta_z4 = dendro_ccz4.scalar("theta_z4","[pp]")

# Lie derivative weight (tensor density weights) 
weight    = - Rational(2,3)
weight_Gh =   Rational(2,3)

# specify the functions for computing first and second derivatives
d = dendro_ccz4.set_first_derivative('grad')    # first arg is direction
d2s = dendro_ccz4.set_second_derivative('grad2')  # first 2 args are directions
ad = dendro_ccz4.set_advective_derivative('agrad')  # first arg is direction
kod = dendro_ccz4.set_kreiss_oliger_dissipation('kograd')

d2 = dendro_ccz4.d2

# generate metric related quantities
dendro_ccz4.set_metric(gt)
igt = dendro_ccz4.get_inverse_metric()

C1 = dendro_ccz4.get_first_christoffel()
C2 = dendro_ccz4.get_second_christoffel()
C3 = dendro_ccz4.get_full_christoffel(psi,p_expo) 

R, R_DZ, R_DZ_for_TF_op, CalGt = dendro_ccz4.compute_ricci(Gh, psi, p_expo)


###################################################################
# Density and current terms from matter. 
# For a simple, vacuum BBH, all values are zero. 
###################################################################
#rho_ADM = 0.0
#Jtd_ADM = Matrix([[0,0,0]])
#S_ADM= Matrix([[0,0,0],[0,0,0],[0,0,0]])


###################################################################
# evolution equations
###################################################################

# The alpha (lapse) equation 
a_rhs = l1*dendro_ccz4.lie(b, a) - 2*a*(K - 2*theta_z4) 

# The beta (shift) equation 
b_rhs = [   S(3)/4 * (lf0 + lf1*a) * B[i] 
          + l2 * dendro_ccz4.vec_j_ad_j(b, b[i]) for i in dendro_ccz4.e_i ] 

# The metric evo eqn 
gt_rhs = dendro_ccz4.lie(b, gt, weight) - 2*a*At 

# The evo eqn for the conformal factor, called psi here 
psi_rhs = dendro_ccz4.lie(b, psi) + Rational(2,3)*psi*( sum([d(i,b[i]) for i in dendro_ccz4.e_i]) - a*K ) / p_expo  

AikAkj = Matrix([sum([At[i, k] * sum([dendro_ccz4.inv_metric[k, l]*At[l, j] for l in dendro_ccz4.e_i]) for k in dendro_ccz4.e_i]) for i, j in dendro_ccz4.e_ij])

#psi_to_neg_p_expo = psi**(-p_expo) 

# The (conformal, traceless) extrinsic curvature (A tilde) evo eqn 
At_rhs = dendro_ccz4.lie(b, At, weight) + (psi**(-p_expo))*dendro_ccz4.trace_free( a*R_DZ_for_TF_op - dendro_ccz4.DiDj(a)) + a*((K-2*theta_z4)*At - 2*AikAkj.reshape(3, 3)) 
#At_rhs = dendro_ccz4.lie(b, At, weight) + (psi_to_neg_p_expo)*dendro_ccz4.trace_free( a*R_DZ_for_TF_op - dendro_ccz4.DiDj(a)) + a*((K-2*theta_z4)*At - 2*AikAkj.reshape(3, 3)) 

# The trK evo eqn 
K_rhs = dendro_ccz4.lie(b, K) - dendro_ccz4.laplacian(a,psi,p_expo) + a*( psi**(-p_expo) * sum([ igt[i,j]*R_DZ[i,j] for i,j in dendro_ccz4.e_ij ]) + K*(K - 2.0*theta_z4) - 3*k1*(1+k2)*theta_z4) 
#K_rhs = dendro_ccz4.lie(b, K) - dendro_ccz4.laplacian(a,psi,p_expo) + a*( psi_to_neg_p_expo * sum([ igt[i,j]*R_DZ[i,j] for i,j in dendro_ccz4.e_ij ]) + K*(K - 2.0*theta_z4) - 3*k1*(1+k2)*theta_z4 ) 

# Intermediate variable used in the Gh eqn 
At_UU = dendro_ccz4.up_up(At)

# The evo eqn for Gh 
Gh_rhs = Matrix([sum(b[j]*ad(j,Gh[i]) for j in dendro_ccz4.e_i) for i in dendro_ccz4.e_i]) - \
         Matrix([sum(Gh[j]*d(j,b[i]) for j in dendro_ccz4.e_i) for i in dendro_ccz4.e_i]) + \
         Rational(2,3)*Matrix([ Gh[i] * sum(d(j,b[j]) for j in dendro_ccz4.e_i)  for i in dendro_ccz4.e_i ]) + \
         Matrix([sum([igt[j, k] * d2(j, k, b[i]) + igt[i, j] * d2(j, k, b[k])/3 for j, k in dendro_ccz4.e_ij]) for i in dendro_ccz4.e_i]) - \
         Matrix([sum([2*At_UU[i, j]*d(j, a) for j in dendro_ccz4.e_i]) for i in dendro_ccz4.e_i]) + \
         Matrix([sum([2*a*dendro_ccz4.C2[i, j, k]*At_UU[j, k] for j,k in dendro_ccz4.e_ij]) for i in dendro_ccz4.e_i]) + \
         Matrix([a*sum([( 3*p_expo/psi*At_UU[i,j]*d(j, psi) - Rational(4,3)*dendro_ccz4.inv_metric[i, j]*d(j, K)) for j in dendro_ccz4.e_i]) for i in dendro_ccz4.e_i]) + \
         Matrix([ 2*( sum([ igt[i,j]*(a*d(j,theta_z4) - theta_z4*d(j,a)) for j in dendro_ccz4.e_i]) - 0.5*a*(k1+Rational(2,3)*K)*(Gh[i]-CalGt[i]) ) for i in dendro_ccz4.e_i ])

Gh_rhs = [item for sublist in Gh_rhs.tolist() for item in sublist]

B_rhs = [Gh_rhs[i] - eta * B[i] +
         l3 * dendro_ccz4.vec_j_ad_j(b, B[i]) -
         l4 * dendro_ccz4.vec_j_ad_j(b, Gh[i]) 
         for i in dendro_ccz4.e_i]

theta_z4_rhs = dendro_ccz4.lie(b, theta_z4) + 0.5*a * ( psi**(-p_expo)*sum([ igt[i,j] * R_DZ[i,j] for i,j in dendro_ccz4.e_ij ]) + 2*K*K/3 - dendro_ccz4.sqr(At) - 2*theta_z4*K) - 0.5*psi**(-p_expo)*sum([ (Gh[i]-CalGt[i])*d(i,a) for i in dendro_ccz4.e_i ]) - a*( k1*(2+k2)*theta_z4 )
#theta_z4_rhs = dendro_ccz4.lie(b, theta_z4) + 0.5*a * (psi_to_neg_p_expo*sum([ igt[i,j] * R_DZ[i,j] for i,j in dendro_ccz4.e_ij ]) + 2*K*K/3 - dendro_ccz4.sqr(At) - 2*theta_z4*K) - 0.5*(psi_to_neg_p_expo)*sum([ (Gh[i]-CalGt[i])*d(i,a) for i in dendro_ccz4.e_i ]) - a*( k1*(2+k2)*theta_z4 )



###################################################################
# generate code
###################################################################

outs = [a_rhs, b_rhs, gt_rhs, psi_rhs, At_rhs, K_rhs, Gh_rhs, B_rhs, theta_z4_rhs]
vnames = ['a_rhs', 'b_rhs', 'gt_rhs', 'psi_rhs', 'At_rhs', 'K_rhs', 'Gh_rhs', 'B_rhs', 'theta_z4_rhs']
#dendro_ccz4.generate_debug(outs, vnames)
dendro_ccz4.generate(outs, vnames, '[pp]')
#numVars=len(outs)
#for i in range(0,numVars):
#    dendro_ccz4.generate_separate([outs[i]],[vnames[i]],'[pp]')
