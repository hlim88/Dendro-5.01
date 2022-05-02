####################################################################
# May.2020
# Quad grav rhs generator new version
#####################################################################

import dendro
from sympy import *

###################################################################
# initialize
###################################################################

l1, l2, l3, l4, eta = symbols('lambda[0] lambda[1] lambda[2] lambda[3] eta')
lf0, lf1 = symbols('lambda_f[0] lambda_f[1]')

#QG related constants
a_const = symbols('a_const')
b_const = symbols('b_const')
qg_mass0_sq = symbols('qg_mass0_sq')
qg_mass2_sq = symbols('qg_mass2_sq')

PI = 3.14159265358979323846
kappa = 1/(16*PI)

# Additional parameters for damping term
R0 = symbols('QUADGRAV_ETA_R0')
ep1, ep2 = symbols('QUADGRAV_ETA_POWER[0] QUADGRAV_ETA_POWER[1]')

# declare variables (BSSN vars)
a   = dendro.scalar("alpha", "[pp]")
chi = dendro.scalar("chi", "[pp]")
K   = dendro.scalar("K", "[pp]")

Gt  = dendro.vec3("Gt", "[pp]")
b   = dendro.vec3("beta", "[pp]")
B   = dendro.vec3("B", "[pp]")

gt  = dendro.sym_3x3("gt", "[pp]")
At  = dendro.sym_3x3("At", "[pp]")

Gt_rhs  = dendro.vec3("Gt_rhs", "[pp]")

# Ricci scalar, R 
Rsc = dendro.scalar("Rsc", "[pp]")
# Aux Ricci scalar, R^
Rsch = dendro.scalar("Rsch", "[pp]")


#TODO : Clear up for documentation
# Ricci tensor, R_ab
#Rab = dendro.sym_3x3("Rab", "[pp]")
# Aux Ricci tensor, V_ab
#Vab = dendro.sym_3x3("Vab", "[pp]")

# Spatial projection of Ricci tensor related quantities
# From R_ab 
Atr = dendro.scalar("Atr", "[pp]")
Aij  = dendro.sym_3x3("Aij", "[pp]")
# From V_ab 
Btr = dendro.scalar("Btr", "[pp]")
Bij  = dendro.sym_3x3("Bij", "[pp]")

# Additional constraint as evolutions vars
Ci = dendro.vec3("Ci","[pp]")
#Ei = dendro.vec3("Ei","[pp]")

# Lie derivative weight
weight = -Rational(2,3)
weight_Gt = Rational(2,3)

# specify the functions for computing first and second derivatives
d = dendro.set_first_derivative('grad')    # first argument is direction
d2s = dendro.set_second_derivative('grad2')  # first 2 arguments are directions
ad = dendro.set_advective_derivative('agrad')  # first argument is direction
kod = dendro.set_kreiss_oliger_dissipation('kograd')
#covD1 = dendro.DiOd('covd1')
#covD2 = dendro.DiTd('covd2')

d2 = dendro.d2

#f = Function('f')

# generate metric related quantities
dendro.set_metric(gt)
igt = dendro.get_inverse_metric()

C1 = dendro.get_first_christoffel()
C2 = dendro.get_second_christoffel()
#what's this...tried to comment it out and python compilation fails
C2_spatial = dendro.get_complete_christoffel(chi)
R, Rt, Rphi, CalGt = dendro.compute_ricci(Gt, chi)
###################################################################
# evolution equations
###################################################################

# Gauge part
# For lapse, 1+log. For shift, modified gamma driver
a_rhs = l1*dendro.lie(b, a) - 2*a*K + 0*dendro.kodiss(a)

b_rhs = [ S(3)/4 * (lf0 + lf1*a) * B[i] +
        l2 * dendro.vec_j_ad_j(b, b[i])
         for i in dendro.e_i ] + 0*dendro.kodiss(b)

eta_func = R0*sqrt(sum([igt[i,j]*d(i,chi)*d(j,chi) for i,j in dendro.e_ij]))/((1-chi**ep1)**ep2)

B_rhs = [Gt_rhs[i] - eta_func * B[i] +
         l3 * dendro.vec_j_ad_j(b, B[i]) -
         l4 * dendro.vec_j_ad_j(b, Gt[i]) + 0*kod(i,B[i])
         for i in dendro.e_i]

# Metric and extrinsic curvature

gt_rhs = dendro.lie(b, gt, weight) - 2*a*At + 0*dendro.kodiss(gt)

chi_rhs = dendro.lie(b, chi, weight) + Rational(2,3) * (chi*a*K) + 0*dendro.kodiss(chi)

AikAkj = Matrix([sum([At[i, k] * sum([dendro.inv_metric[k, l]*At[l, j] for l in dendro.e_i]) for k in dendro.e_i]) for i, j in dendro.e_ij])

#NOTE : "CAUTION" THIS IS DIFFERENT THEN Atr for Ricci. 
#NOTE : The rho, S, and Sij in BSSN eqns are not "physical" matter.
#Introduce not conformally transformed metric again
gs = gt/chi
igs = igt/chi

# Define additional constraints
# NOTE : Ci is considered as evolution variable but Ei is treated algebraically
AiUjD = dendro.up_down(Aij)*chi
Ei = - Matrix([[K*Ci[0], K*Ci[1], K*Ci[2]]])
# Ei is determined by the spatial projection of RHS of Eqn.47

rho_qg = Rsc/4
Si_qg = Matrix([[ -Ci[0], -Ci[1], -Ci[2]]])
Sij_qg = Matrix([Aij[i,j] + gs[i,j]*Atr/3 + gs[i,j]*Rsc/4 for i,j in dendro.e_ij]).reshape(3,3)
S_qg = sum([sum([Sij_qg[i,j]*igs[i,j] for i in dendro.e_i]) for j in dendro.e_i])

At_rhs = dendro.lie(b, At, weight) + chi*dendro.trace_free( a*R - dendro.DiDj(a)-8*pi*Sij_qg) + a*(K*At - 2*AikAkj.reshape(3, 3)) + 0*dendro.kodiss(At)

K_rhs = dendro.lie(b, K) - dendro.laplacian(a,chi) + a*(K*K/3 + dendro.sqr(At)) + 4*pi*a*(rho_qg + S_qg) + 0*dendro.kodiss(K)

At_UU = dendro.up_up(At)

Gt_rhs = Matrix([sum(b[j]*ad(j,Gt[i]) for j in dendro.e_i) for i in dendro.e_i]) - \
         Matrix([sum(CalGt[j]*d(j,b[i]) for j in dendro.e_i) for i in dendro.e_i]) + \
         Rational(2,3)*Matrix([ CalGt[i] * sum(d(j,b[j]) for j in dendro.e_i)  for i in dendro.e_i ]) + \
         Matrix([sum([igt[j, k] * d2(j, k, b[i]) + igt[i, j] * d2(j, k, b[k])/3 for j, k in dendro.e_ij]) for i in dendro.e_i]) - \
         Matrix([sum([2*At_UU[i, j]*d(j, a) for j in dendro.e_i]) for i in dendro.e_i]) + \
         Matrix([sum([2*a*dendro.C2[i, j, k]*At_UU[j, k] for j,k in dendro.e_ij]) for i in dendro.e_i]) - \
         Matrix([sum([a*(3/chi*At_UU[i,j]*d(j, chi) + Rational(4,3)*dendro.inv_metric[i, j]*d(j, K)) for j in dendro.e_i]) for i in dendro.e_i])
         # + kod(i,Gt[i])

Gt_rhs = [item for sublist in Gt_rhs.tolist() for item in sublist]

# Ricci tensor and scalar

# normal vector n^a
n_vec = Matrix([[1/a, -b[0]/a, -b[1]/a, -b[2]/a]])

# Define acceleration (n^c \del_c n_a) HL : this is equal to 1/a*D_i a
a_acc = Matrix([[1,0,0]])

# Some precomputation/definition

Aij_UU = dendro.up_up(Aij)*(chi*chi)
Ci_U = Matrix([sum([Ci[j]*igs[i,j] for j in dendro.e_i]) for i in dendro.e_i])
Kij = At + 1/3*gt*K
Kki = dendro.up_down(Kij)*chi

# Ricci scalar, R
# Eqn.23
Rsc_rhs = dendro.lie(b, Rsc) - a*Rsch

# Aux Ricci scalar, R^
# Eqn.24
Rsch_rhs = dendro.lie(b, Rsch) -  K*Rsch - qg_mass0_sq*Rsc - 2*0*(rho_qg - S_qg)

# Ricci tensor
# From R_ab
# Eqn.38
Atr_rhs = dendro.lie(b, Atr) + a*(2*sum([a_acc[i]*Ci[i] for i in dendro.e_i]) - Btr) 
# Eqn.39
Aij_rhs = a*(2/3*gt*sum([a_acc[k]*Ci[k] for k in dendro.e_i]) - Bij) 

# From V_ab
# Eqn.42
#Btr_rhs = dendro.lie(b, Btr) +2*a*sum([a_acc[k]*Ei[k] for k in dendro.e_i]) 
Btr_rhs = dendro.lie(b, Btr) 

# Eqn.43
Bij_rhs = -Btr*Kij

#Eqn.41
Ci_rhs = Matrix([sum(b[j]*ad(j,Ci[i]) - Ei[i] for j in dendro.e_i) for i in dendro.e_i])
Ci_rhs = [item for sublist in Ci_rhs.tolist() for item in sublist]

#RHS of Eqn.43, same argument from Aij_rhs is applicable for this 

#TODO : Additional constraints, C_k, E_k go here if we want to evolve and monitor

#_I = gt*igt
#print(simplify(_I))

#_I = gt*dendro.inv_metric
#print(simplify(_I))


###
# Substitute ...
#for expr in [a_rhs, b_rhs[0], b_rhs[1], b_rhs[2], B_rhs[0], B_rhs[1], B_rhs[2], K_rhs, chi_rhs, Gt_rhs[0], Gt_rhs[1], Gt_rhs[2], gt_rhs[0], gt_rhs[0,0], gt_rhs[1,1], gt_rhs[2,2], gt_rhs[0,1], gt_rhs[0,2], gt_rhs[1,2], At_rhs[0,0], At_rhs[0,1], At_rhs[0,2], At_rhs[1,1], At_rhs[1,2], At_rhs[2,2]]:
#    for var in [a, b[0], b[1], b[2], B[0], B[1], B[2], chi, K, gt[0,0], gt[0,1], gt[0,2], gt[1,1], gt[1,2], gt[2,2], Gt[0], Gt[1], Gt[2], At[0,0], At[0,1], At[0,2], At[1,1], At[1,2], At[2,2]]:
#        expr.subs(d2(1,0,var), d2(0,1,var))
#        expr.subs(d2(2,1,var), d2(1,2,var))
#        expr.subs(d2(2,0,var), d2(0,2,var))
#
#print (a_rhs)
#print (G_rhs)


###################################################################
# generate code
###################################################################

outs = [a_rhs, b_rhs, gt_rhs, chi_rhs, At_rhs, K_rhs, Gt_rhs, B_rhs, Rsc_rhs, Rsch_rhs, Atr_rhs, Aij_rhs, Btr_rhs, Bij_rhs]
vnames = ['a_rhs', 'b_rhs', 'gt_rhs', 'chi_rhs', 'At_rhs', 'K_rhs', 'Gt_rhs', 'B_rhs', 'Rsc_rhs','Rsch_rhs','Atr_rhs', 'Aij_rhs', 'Btr_rhs', 'Bij_rhs']
#dendro.generate_debug(outs, vnames)
dendro.generate(outs, vnames, '[pp]')
#numVars=len(outs)
#for i in range(0,numVars):
#    dendro.generate_separate([outs[i]],[vnames[i]],'[pp]')
