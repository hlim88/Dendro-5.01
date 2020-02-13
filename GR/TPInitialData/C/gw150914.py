####################################################################
#
#  This python script calculates the constants used by the tpid
#  code to generate initial data for the GW150914 event.
# 
#  This script prints C variables for the code. You have to change
#  [] to {} by hand.
#
####################################################################
from math import *

D = 10.0                    # Separation
q = 36.0/29.0               # Mass Ratio
M = 1.0                     # Total mass
chip = [0, 0,  0.31]        # Dimensionsless spin of  + BH (x0>0, more massive)
chim = [0, 0, -0.46]        # Dimensionsless spin of - BH (x0<0, less massive)
Pr = -0.00084541526517121   # Radial linear momentum
Pphi = 0.09530152296974252  # Azimuthal linear momentum

mp = M * q/(1.0+q)
mm = M * 1.0/(1.0+q)
xp = D * mm
xm = -D * mp
half_D = D / 2.0
center_offset = [xp - half_D, 0, 0]

Sp = [x * mp**2 for x in chip]
Sm = [x * mm**2 for x in chim]

Pp = [Pr, Pphi, 0.0]
Pm = [x * (-1.0) for x in Pp]

print("  const double target_M_plus = ", mp, ";")
print("  const double target_M_minus = ", mm, ";")
print("  double par_m_plus = ", mp, ";")
print("  double par_m_minus = ", mm, ";")
print("  const double par_b = ", half_D, ";")
print("  const double center_offset[3] = ", center_offset, ";")
print("  const double par_P_plus[3] = ", Pp, ";")
print("  const double par_P_minus[3] = ", Pm, ";")
print("  const double par_S_plus[3] = ", Sp, ";")
print("  const double par_S_minus[3] = ", Sm, ";")

