####################################################################
# Mar.16.2018
# Python script for generating CCZ4 rhs
# Note that we use psi function instead of chi to get more options
####################################################################

import dendro_ccz4
from sympy import *
from sympy.physics.vector.vector import Vector
from sympy.printing.dot import dotprint
###################################################################
# initialize
###################################################################

l1, l2, l3, l4, eta = symbols('lambda[0] lambda[1] lambda[2] lambda[3] eta')
lf0, lf1 = symbols('lambda_f[0] lambda_f[1]')
k1, k2 = symbols('kappa[0] kappa[1]')
p_expo = symbols('p_expo') # Exponential for conformal factor

# declare variables
a   = dendro_ccz4.scalar("alpha", "[pp]")
psi = dendro_ccz4.scalar("psi", "[pp]") # This is a supplement for psi
K   = dendro_ccz4.scalar("K", "[pp]")

Gh  = dendro_ccz4.vec3("Gh", "[pp]")
b   = dendro_ccz4.vec3("beta", "[pp]")
B   = dendro_ccz4.vec3("B", "[pp]")

gt  = dendro_ccz4.sym_3x3("gt", "[pp]")
At  = dendro_ccz4.sym_3x3("At", "[pp]")


# note: these are just the symbolic vars that is being used to generate the
# Gh_rhs by stages

_Gh_rhs_s1  = dendro_ccz4.vec3("Gh_rhs_s1_", "[pp]")
_Gh_rhs_s2  = dendro_ccz4.vec3("Gh_rhs_s2_", "[pp]")
_Gh_rhs_s3  = dendro_ccz4.vec3("Gh_rhs_s3_", "[pp]")
_Gh_rhs_s4  = dendro_ccz4.vec3("Gh_rhs_s4_", "[pp]")
_Gh_rhs_s5  = dendro_ccz4.vec3("Gh_rhs_s5_", "[pp]")
_Gh_rhs_s6  = dendro_ccz4.vec3("Gh_rhs_s6_", "[pp]")
_Gh_rhs_s7  = dendro_ccz4.vec3("Gh_rhs_s7_", "[pp]")
_Gh_rhs_s8  = dendro_ccz4.vec3("Gh_rhs_s8_", "[pp]")
_CalGt  = dendro_ccz4.vec3("CalGt", "[pp]")
_Gh_rhs  = dendro_ccz4.vec3("Gh_rhs", "[pp]")


# Gh_rhs staged vars that is being used to generate the code.
At_UU  = dendro_ccz4.sym_3x3("At_UU", "[pp]")
CalGt  = dendro_ccz4.vec3("CalGt", "[pp]")
Gh_rhs_s1  = dendro_ccz4.vec3("Gh_rhs_s1_", "[pp]")
Gh_rhs_s2  = dendro_ccz4.vec3("Gh_rhs_s2_", "[pp]")
Gh_rhs_s3  = dendro_ccz4.vec3("Gh_rhs_s3_", "[pp]")
Gh_rhs_s4  = dendro_ccz4.vec3("Gh_rhs_s4_", "[pp]")
Gh_rhs_s5  = dendro_ccz4.vec3("Gh_rhs_s5_", "[pp]")
Gh_rhs_s6  = dendro_ccz4.vec3("Gh_rhs_s6_", "[pp]")
Gh_rhs_s7  = dendro_ccz4.vec3("Gh_rhs_s7_", "[pp]")
Gh_rhs_s8  = dendro_ccz4.vec3("Gh_rhs_s8_", "[pp]")

Gh_rhs  = dendro_ccz4.vec3("Gh_rhs", "[pp]")

# Additional variable for ccz4
# HL : I guess we do not need to make stage of theta_z4 variable

theta_z4 = dendro_ccz4.scalar("theta_z4","[pp]")
#pZ4 = dendro_ccz4.vec3("pZ4","[pp]")

# Lie derivative weight
weight = -Rational(2,3)
weight_Gh = Rational(2,3)

# specify the functions for computing first and second derivatives
d = dendro_ccz4.set_first_derivative('grad')    # first argument is direction
d2s = dendro_ccz4.set_second_derivative('grad2')  # first 2 arguments are directions
ad = dendro_ccz4.set_advective_derivative('agrad')  # first argument is direction
kod = dendro_ccz4.set_kreiss_oliger_dissipation('kograd')

d2 = dendro_ccz4.d2

#f = Function('f')

# generate metric related quantities
dendro_ccz4.set_metric(gt)
igt = dendro_ccz4.get_inverse_metric()

C1 = dendro_ccz4.get_first_christoffel()
C2 = dendro_ccz4.get_second_christoffel()
#what's this...tried to comment it out and python compilation fails
C3 = dendro_ccz4.get_full_christoffel(psi, p_expo)
R, R_DZ, R_DZ_for_TF_op, CalGt = dendro_ccz4.compute_ricci(Gh, psi, p_expo)

# Define value pi
# TODO : Just use math package?

###################################################################
# Density and currents terms from matter.
# For simple BBH, All values are zero
###################################################################

rho_ADM = 0.0
Jtd_ADM = Matrix([[0,0,0]])
S_ADM= Matrix([[0,0,0],[0,0,0],[0,0,0]])

###################################################################
# evolution equations
###################################################################
# Projected Z vector 
# HL : This should not be declared as dendro var because we are not
# using them. But required to calculate rhs equations. Define here
# and put zeros initially but need to check
pZ4 = Matrix([[0,0,0]])

# Lapse equation
a_rhs = l1*dendro_ccz4.lie(b, a) - 2*a*(K - 2*theta_z4) 

# Shift equation
b_rhs = [ S(3)/4 * (lf0 + lf1*a) * B[i] +
        l2 * dendro_ccz4.vec_j_ad_j(b, b[i])
         for i in dendro_ccz4.e_i ] 

# Metric evolution equation
gt_rhs = dendro_ccz4.lie(b, gt, weight) - 2*a*At

# Evolution equation for conformal factor, psi for CCZ4
psi_rhs = dendro_ccz4.lie(b, psi) + Rational(2,3)*psi*( sum([d(i,b[i]) for i in dendro_ccz4.e_i]) - a*K ) / p_expo


AikAkj = Matrix([sum([At[i, k] * sum([dendro_ccz4.inv_metric[k, l]*At[l, j] for l in dendro_ccz4.e_i]) for k in dendro_ccz4.e_i]) for i, j in dendro_ccz4.e_ij])

DiZjDjZi = Matrix([d(j,pZ4[j])+d(i,pZ4[j]) - 2*sum([C3[l,i,j]*pZ4[l] for l in dendro_ccz4.e_i]) for i,j in dendro_ccz4.e_ij])

# Extrinsic curvature (conformal, traceless) evolution equation
At_rhs = dendro_ccz4.lie(b, At, weight) + (psi**(-p_expo))*dendro_ccz4.trace_free( a*R_DZ_for_TF_op - dendro_ccz4.DiDj(a)) + a*((K-2*theta_z4)*At - 2*AikAkj.reshape(3, 3))

# Trace of K evolution equation
K_rhs = dendro_ccz4.lie(b, K) - dendro_ccz4.laplacian(a,psi,p_expo) + a*( psi**(-p_expo) * sum([ igt[i,j]*R_DZ[i,j] for i,j in dendro_ccz4.e_ij ]) + K*(K - 2.0*theta_z4) - 3*k1*(1+k2)*theta_z4)

# Intermediate value
At_UU = dendro_ccz4.up_up(At)

Gh_rhs_s1= ([sum(b[j]*ad(j,Gh[i]) for j in dendro_ccz4.e_i) for i in dendro_ccz4.e_i])
Gh_rhs_s2= ([sum(Gh[j]*d(j,b[i]) for j in dendro_ccz4.e_i) for i in dendro_ccz4.e_i])
Gh_rhs_s3= ([Gh[i] * sum(d(j,b[j]) for j in dendro_ccz4.e_i)  for i in dendro_ccz4.e_i ])
Gh_rhs_s4= ([sum([igt[j, k] * d2(j, k, b[i]) + igt[i, j] * d2(j, k, b[k])/3 for j, k in dendro_ccz4.e_ij]) for i in dendro_ccz4.e_i])
Gh_rhs_s5= ([sum([2*At_UU[i, j]*d(j, a) for j in dendro_ccz4.e_i]) for i in dendro_ccz4.e_i])
Gh_rhs_s6= ([sum([2*a*dendro_ccz4.C2[i, j, k]*At_UU[j, k] for j,k in dendro_ccz4.e_ij]) for i in dendro_ccz4.e_i])
Gh_rhs_s7= ([sum([a*(3*p_expo/psi*At_UU[i,j]*d(j, psi) - Rational(4,3)*dendro_ccz4.inv_metric[i, j]*d(j, K)) for j in dendro_ccz4.e_i]) for i in dendro_ccz4.e_i])
Gh_rhs_s8=([2*(sum([igt[i,j]*(a*d(j,theta_z4)-theta_z4*d(j,a)) for j in dendro_ccz4.e_i])-0.5*a*(k1+Rational(2,3)*K)*(Gh[i]-CalGt[i])) for i in dendro_ccz4.e_i])

'''Gh_rhs = Matrix([sum(b[j]*ad(j,Gh[i]) for j in dendro_ccz4.e_i) for i in dendro_ccz4.e_i]) - \
         Matrix([sum(CalGt[j]*d(j,b[i]) for j in dendro_ccz4.e_i) for i in dendro_ccz4.e_i]) + \
         Rational(2,3)*Matrix([ CalGt[i] * sum(d(j,b[j]) for j in dendro_ccz4.e_i)  for i in dendro_ccz4.e_i ]) + \
         Matrix([sum([igt[j, k] * d2(j, k, b[i]) + igt[i, j] * d2(j, k, b[k])/3 for j, k in dendro_ccz4.e_ij]) for i in dendro_ccz4.e_i]) - \
         Matrix([sum([2*At_UU[i, j]*d(j, a) for j in dendro_ccz4.e_i]) for i in dendro_ccz4.e_i]) + \
         Matrix([sum([2*a*dendro_ccz4.C2[i, j, k]*At_UU[j, k] for j,k in dendro_ccz4.e_ij]) for i in dendro_ccz4.e_i]) - \
         Matrix([sum([a*(3/chi*At_UU[i,j]*d(j, chi) + Rational(4,3)*dendro_ccz4.inv_metric[i, j]*d(j, K)) for j in dendro_ccz4.e_i]) for i in dendro_ccz4.e_i])
         # + kod(i,Gh[i])'''


Gh_rhs = Matrix(_Gh_rhs_s1) - \
         Matrix(_Gh_rhs_s2) + \
         Rational(2,3)*Matrix(_Gh_rhs_s3) + \
         Matrix(_Gh_rhs_s4) - \
         Matrix(_Gh_rhs_s5) + \
         Matrix(_Gh_rhs_s6) - \
         Matrix(_Gh_rhs_s7) + \
         Matrix(_Gh_rhs_s8)



# + kod(i,Gh[i])


Gh_rhs = [item for sublist in Gh_rhs.tolist() for item in sublist]

B_rhs = [_Gh_rhs[i] - eta * B[i] +
         l3 * dendro_ccz4.vec_j_ad_j(b, B[i]) -
         l4 * dendro_ccz4.vec_j_ad_j(b, Gh[i]) + 0*kod(i,B[i])
         for i in dendro_ccz4.e_i]

theta_z4_rhs = dendro_ccz4.lie(b, theta_z4) + 0.5*a * ( psi**(-p_expo)*sum([ igt[i,j] * R_DZ[i,j] for i,j in dendro_ccz4.e_ij ]) + 2*K*K/3 - dendro_ccz4.sqr(At) - 2*theta_z4*K) - 0.5*psi**(-p_expo)*sum([ (Gh[i]-CalGt[i])*d(i,a) for i in dendro_ccz4.e_i ]) - a*( k1*(2+k2)*theta_z4 )

#_I = gt*igt
#print(simplify(_I))

#_I = gt*dendro_ccz4.inv_metric
#print(simplify(_I))


###
# Substitute ...
#for expr in [a_rhs, b_rhs[0], b_rhs[1], b_rhs[2], B_rhs[0], B_rhs[1], B_rhs[2], K_rhs, chi_rhs, Gh_rhs[0], Gh_rhs[1], Gh_rhs[2], gt_rhs[0], gt_rhs[0,0], gt_rhs[1,1], gt_rhs[2,2], gt_rhs[0,1], gt_rhs[0,2], gt_rhs[1,2], At_rhs[0,0], At_rhs[0,1], At_rhs[0,2], At_rhs[1,1], At_rhs[1,2], At_rhs[2,2]]:
#    for var in [a, b[0], b[1], b[2], B[0], B[1], B[2], chi, K, gt[0,0], gt[0,1], gt[0,2], gt[1,1], gt[1,2], gt[2,2], Gh[0], Gh[1], Gh[2], At[0,0], At[0,1], At[0,2], At[1,1], At[1,2], At[2,2]]:
#        expr.subs(d2(1,0,var), d2(0,1,var))
#        expr.subs(d2(2,1,var), d2(1,2,var))
#        expr.subs(d2(2,0,var), d2(0,2,var))
#
#print (a_rhs)
#print (G_rhs)

###################################################################
# generate code
###################################################################

outs = [a_rhs, b_rhs, gt_rhs, psi_rhs, At_rhs, K_rhs, CalGt, Gh_rhs_s1, Gh_rhs_s2, Gh_rhs_s3, Gh_rhs_s4, Gh_rhs_s5, Gh_rhs_s6, Gh_rhs_s7, Gh_rhs_s8, Gh_rhs, B_rhs, theta_z4_rhs]
vnames = ['a_rhs', 'b_rhs', 'gt_rhs', 'psi_rhs', 'At_rhs', 'K_rhs', 'CalGt', 'Gh_rhs_s1_', 'Gh_rhs_s2_', 'Gh_rhs_s3_', 'Gh_rhs_s4_', 'Gh_rhs_s5_', 'Gh_rhs_s6_', 'Gh_rhs_s7_', 'Gh_rhs_s8_', 'Gh_rhs', 'B_rhs', 'theta_z4_rhs']
#dendro_ccz4.generate_debug(outs, vnames)
#dendro_ccz4.generate(outs, vnames, '[pp]')
numVars=len(outs)
for i in range(0,numVars):
    dendro_ccz4.generate_separate([outs[i]],[vnames[i]],'[pp]')

#dendro_ccz4.generate_separate([Gh_rhs],['Gh_rhs'],'[pp]')
#dendro_ccz4.generate([CalGt, Gh_rhs_s1, Gh_rhs_s2, Gh_rhs_s3, Gh_rhs_s4, Gh_rhs_s5, Gh_rhs_s6, Gh_rhs_s7, Gh_rhs,B_rhs],['CalGt', 'Gh_rhs_s1', 'Gh_rhs_s2', 'Gh_rhs_s3', 'Gh_rhs_s4', 'Gh_rhs_s5', 'Gh_rhs_s6', 'Gh_rhs_s7', 'Gh_rhs', 'B_rhs'],'[pp]')
