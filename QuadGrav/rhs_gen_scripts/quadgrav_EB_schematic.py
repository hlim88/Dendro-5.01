####################################################################
# 19-Feb 2024
# form of electric and magnetic part (EE and BB) of the Weyl tensor
# in notation of other quadgrav scripts (see quadgrav_const.py)
#####################################################################

## declare variables to hold the electric and magnetic part of the Weyl tensor
## probably we do not need this, since we do not treat these as evolution variables
#EE = dendro.sym_3x3("EE", "[pp]")
#BB = dendro.sym_3x3("BB", "[pp]")

# TODO: @Hyun: let's both check the implementation, also against the equatios in the current notes


# import LeviCivita
from sympy import LeviCivita


# determine EE (electric Weyl) from other 3+1 variables
# TODO: @Hyun: please check how Rij is implemented

EEij_part1 = K*Kij

EEij_part2 = Matrix([
	sum([
		Kij[i,m]*Kki[m,j] 
		for m in dendro.e_i
	]) 
	for i,j in dendro.e_ij
]).reshape(3,3)

EEij_part3 = Rij # should be the *3D spatial* Ricci scalar

EEij_part4 = -1/2*(
	Aij 
	+ 4/3*Atr*gt
	+ 1/3*Rsc*gt
)

EEij = EEij_part1 + EEij_part2 + EEij_part3 + EEij_part4


# determine BB (magnetic Weyl) from other 3+1 variables
# TODO: @Hyun: not sure if this call to Levi-Civita works but I tested the scipy implementation and it seems to work for any set of adjacent real numbers

BBij = Matrix([
	sum([
		sum([
			1/2*(LeviCivita[k,l,i]*DiKkj[k,j,l] + LeviCivita[k,l,j]*DiKkj[k,i,l])
			for k in dendro.e_i
		]) 
		for l in dendro.e_i
	]) 
	for i,j in dendro.e_ij
]).reshape(3,3)


















