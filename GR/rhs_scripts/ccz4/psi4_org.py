####################################################################
#
# Date : Dec.12.2017
# Python script that generates Psi4 for Gravitational waves
# 
####################################################################

#!/usr/bin/env/ python3

import dendro
from sympy import *

###################################################################
# initialize
###################################################################

# Declare variables
# This is BSSN variables that we need for Psi4 calculation
chi = dendro.scalar("chi","[pp]")
K   = dendro.scalar("K","[pp]")
Gt  = dendro.vec3("Gt","[pp]")
gt  = dendro.sym_3x3("gt","[pp]")
At  = dendro.sym_3x3("At","[pp]")

# Lie derivative weight
weight = -2/3
weight_Gt = 2/3

# specify the functions for computing first and second derivatives
d = dendro.set_first_derivative('grad')    # first argument is direction
d2 = dendro.set_second_derivative('grad2')  # first 2 arguments are directions
ad = dendro.set_advective_derivative('agrad')  # first argument is direction

# generate metric related quantities 
dendro.set_metric(gt)
igt = dendro.get_inverse_metric()

C1 = dendro.get_first_christoffel()
C2 = dendro.get_second_christoffel()
#what's this...tried to comment it out and python compilation fails 
C2_spatial = dendro.get_complete_christoffel(chi)
R, Rt, Rphi, CalGt = dendro.compute_ricci(Gt, chi)

###################################################################
# Calculate teatrad used in Psi4 calculation
####################################################################

# Define coordinates
x, y, z = symbols('x, y, z')

# Some other values
invsqrt2 = 0.7071067811865475244
inv_chi = 1/chi

# Define vectors for tetrad
r_vec = Matrix([[x,y,z]])
theta = Matrix([[x*z,y*z,-(x*x+y*y)]])
phi = Matrix([[-y,x,0.0]])

# Using Gram-Schmidt to make orthonormal basis. Here we use the non-conformally scaled metric
gd = gt*inv_chi


# For r_vec
inner_product = 0.0
inner_product = sum([sum([gd[i,j] * r_vec[i] * r_vec[j] for i in dendro.e_i]) for j in dendro.e_i])

r_vec /= sqrt(inner_product)

# For theta
inner_product_1 = 0.0
inner_product_2 = 0.0

inner_product_1 = sum([sum([gd[i,j] * theta[i] * theta[j] for i in dendro.e_i]) for j in dendro.e_i])
inner_product_2 = sum([sum([gd[i,j] * r_vec[i] * theta[j] for i in dendro.e_i]) for j in dendro.e_i])
theta -= inner_product_2 * r_vec
theta /= sqrt(inner_product_1 - inner_product_2 * inner_product_2)


# For phi
inner_product_1 = 0.0
inner_product_2 = 0.0
inner_product_3 = 0.0

inner_product_1 = sum([sum([gd[i,j] * phi[i] * phi[j] for i in dendro.e_i]) for j in dendro.e_i])
inner_product_2 = sum([sum([gd[i,j] * r_vec[i] * phi[j] for i in dendro.e_i]) for j in dendro.e_i])
inner_product_3 = sum([sum([gd[i,j] * theta[i] * phi[j] for i in dendro.e_i]) for j in dendro.e_i])

phi -= inner_product_2 * r_vec + inner_product_3 * theta
phi /= sqrt(inner_product_1 - inner_product_2 * inner_product_2 - inner_product_3 * inner_product_3)

# End of tetrad construction. The relevant pieces of the tetrad then beome (note that these
# could be viewed as 4D vectors with zero timelike components and each is contravariant form)

###################################################################
# Calculate the Weyl scalar, Psi4 for graviational wave extraction
###################################################################

# Tetrad variables for calculate psi-4
r_np = Matrix([[r_vec[0],r_vec[1],r_vec[2]]])
m_np_real = Matrix([[theta[0],theta[1],theta[2]]])*invsqrt2
m_np_img = Matrix([[phi[0],phi[1],phi[2]]])*invsqrt2

At_UD = dendro.up_down(At)

# Some auxilary variables
# MM and NN are 3x3 symmetric matrices and
# MR and NR are 3x3 anti-symmetric matrices

MM = Matrix([m_np_real[i]*m_np_real[j] - m_np_img[i]*m_np_img[j] for i,j in dendro.e_ij]) 
MM = MM.reshape(3,3)
NN = Matrix([m_np_real[i]*m_np_img[j] + m_np_real[j]*m_np_img[i] for i,j in dendro.e_ij])
NN = NN.reshape(3,3)  
MR = Matrix([m_np_real[i]*r_np[j] - m_np_real[j]*r_np[i] for i,j in dendro.e_ij])
MR = MR.reshape(3,3) 
NR = Matrix([m_np_img[i]*r_np[j] - m_np_img[j]*r_np[i] for i,j in dendro.e_ij])
NR = NR.reshape(3,3)  

# Additional intermediate variables

A_vec = Matrix([[sum([At[j,0]*r_np[j] for j in dendro.e_i]), sum([At[j,1]*r_np[j] for j in dendro.e_i]),sum([At[j,2]*r_np[j] for j in dendro.e_i])]])

Uu = Matrix([sum([m_np_real[k] * (d(j, At[k,i]) + sum([C1[m,k,i] * At_UD[m,j] for m in dendro.e_i])) for k in dendro.e_i]) for i,j in dendro.e_ij])
Uu = Uu.reshape(3,3)
Vv = Matrix([sum([m_np_img[k] * (d(j, At[k,i]) + sum([C1[m,k,i] * At_UD[m,j] for m in dendro.e_i])) for k in dendro.e_i]) for i,j in dendro.e_ij])
Vv = Vv.reshape(3,3)

r_d_chi = sum([r_np[i] * d(i, chi) for i in dendro.e_i]) 

A_temp = inv_chi * inv_chi * ( sum([A_vec[i] * r_np[i] for i in dendro.e_i]) +1.0/3.0 * K * chi + 0.5 * r_d_chi ) 

m_real_d_chi = sum([m_np_real[i] * d(i, chi) for i in dendro.e_i])  
m_img_d_chi = sum([m_np_img[i] * d(i, chi) for i in dendro.e_i]) 

m_real_A_vec = sum([m_np_real[i] * A_vec[i] for i in dendro.e_i]) 
m_img_A_vec = sum([m_np_img[i] * A_vec[i] for i in dendro.e_i])  


# Calculate Psi4

psi4_1_real = sum([(inv_chi * Rphi[i,i] + Rt[i,i]) * MM[i,i] for i in dendro.e_i]) + 2*sum([sum([(inv_chi * Rphi[i,j] + Rt[i,j]) * MM[i,j] for j in range(i+1,3)]) for i in range(0,2)]) 
psi4_1_img = sum([(inv_chi * Rphi[i,i] + Rt[i,i]) * NN[i,i] for i in dendro.e_i]) + 2*sum([sum([(inv_chi * Rphi[i,j] + Rt[i,j]) * NN[i,j] for j in range(i+1,3)]) for i in range(0,2)]) 

psi4_2_real = A_temp * (sum([At[i,i] * MM[i,i] for i in dendro.e_i]) + 2*sum([sum([At[i,j] * MM[i,j] for j in range(i+1,3)]) for i in range(0,2)]))
psi4_2_img = A_temp * (sum([At[i,i] * NN[i,i] for i in dendro.e_i]) + 2*sum([sum([At[i,j] * NN[i,j] for j in range(i+1,3)]) for i in range(0,2)]))  

psi4_3_real = inv_chi * sum([sum([MR[i,j]* Uu[i,j] - NR[i,j]*Vv[i,j] for i in dendro.e_i]) for j in dendro.e_i])  
psi4_3_img = inv_chi * sum([sum([NR[i,j]* Uu[i,j] + MR[i,j]*Vv[i,j] for i in dendro.e_i]) for j in dendro.e_i]) 

psi4_4_real = inv_chi * inv_chi * (m_real_A_vec * (m_real_A_vec + 0.5 * m_real_d_chi) - m_img_A_vec * (m_img_A_vec + 0.5 * m_img_d_chi))  
psi4_4_img = inv_chi * inv_chi * (m_real_A_vec * (m_img_A_vec - 0.5 * m_img_d_chi ) + m_img_A_vec * (m_real_A_vec - 0.5 * m_real_d_chi))  

# Adding previous auxilary Psi4 calculations

psi4_real = psi4_1_real + psi4_2_real - psi4_3_real - psi4_4_real
psi4_img = - (psi4_1_img + psi4_2_img - psi4_3_img - psi4_4_img)

# Adding constraints

# Output for this should be included psi4_real and psi4_img as double precision  

###################################################################
# generate code
###################################################################

outs = [psi4_real, psi4_img]
vnames = ['psi4_real', 'psi4_img']
dendro.generate(outs, vnames, '[pp]')
