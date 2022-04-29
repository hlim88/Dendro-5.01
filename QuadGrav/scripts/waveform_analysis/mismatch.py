import numpy as np
import lalsimutils
import lal
import lalsimulation as lalsim

import geometric_to_MKS
import re

SEG_LENGTH=4.0  # length of segment in seconds

def nextPow2(length):        
    """                                                  
    Find next power of 2 <= length of data                                                                              
    """                                                    
    return int(2**np.ceil(np.log2(length)))

def str2z(s):
    s = s.decode("utf-8")
    ss = re.sub('[\(\)]','',s)
    sr, si = ss.split(',')
    return complex(float(sr), float(si))

def str2r(s):
    s = s.decode("utf-8")
    ss = re.sub('[\(\)]','',s)
    sr, si = ss.split(',')
    return float(sr)

def str2i(s):
    s = s.decode("utf-8")
    ss = re.sub('[\(\)]','',s)
    sr, si = ss.split(',')
    return float(si)

def overlap(wave1, wave2, mass):
    """
    Purpose: Calculate the overlap of two waveforms. Both waveforms are assumed
    to be psi_4 rather than the strain.
    wave1, wave2: Dictionaries with entries "time", and "psi4" corresponding to
    a uniform time grid and the complex psi_4, respectively.
    mass: the total mass in units of solar masses
    """

    # The mass normalization only affects the time array (any multiplicative
    # change to psi_4 is canceled and can therefore be ignored

    wave1_nt = geometric_to_MKS.GeometricTime_To_MKS_Time(wave1["time"],
            mass*geometric_to_MKS.msun)
    wave2_nt = geometric_to_MKS.GeometricTime_To_MKS_Time(wave2["time"],
            mass*geometric_to_MKS.msun)
  

    # both waveforms are interpolated onto a standard grid
    t_final = np.min((wave1_nt[-1], wave2_nt[-1]))
    dt = np.min((wave1_nt[1]-wave1_nt[0], wave2_nt[1]-wave2_nt[0]))
    npts = int(t_final / dt)
  
    tgrid = np.linspace(0.0, npts * dt, npts+1)
    wave1_psi4_grid = np.interp(tgrid, wave1_nt, wave1["psi4"])
    wave2_psi4_grid = np.interp(tgrid, wave2_nt, wave2["psi4"])
  
    Psi4T1 = lal.CreateCOMPLEX16TimeSeries("Complex overlap",lal.LIGOTimeGPS(0.), 0., dt, lal.DimensionlessUnit, npts+1)
    Psi4T2 = lal.CreateCOMPLEX16TimeSeries("Complex overlap",lal.LIGOTimeGPS(0.), 0., dt, lal.DimensionlessUnit, npts+1)
  
    Psi4T1.data.data = wave1_psi4_grid
    Psi4T2.data.data = wave2_psi4_grid
  
    # The final time is set by the desired segment length and the restriction that
    # the number of points must be a power of 2. We enlarge the grids by adding
    # zeros
    nsize = int(SEG_LENGTH / dt)
    pow2npts = nextPow2(nsize)
    Psi4T1 = lal.ResizeCOMPLEX16TimeSeries(Psi4T1, 0, pow2npts)
    Psi4T2 = lal.ResizeCOMPLEX16TimeSeries(Psi4T2, 0, pow2npts)


    Psi4F1 = lalsimutils.DataFourier(Psi4T1)
    Psi4F2 = lalsimutils.DataFourier(Psi4T2)
  
  
    Tlen = pow2npts * dt
    psd=lalsim.SimNoisePSDaLIGOZeroDetHighPower
    IP = lalsimutils.CreateCompatibleComplexOverlap(Psi4F1,deltaF = 1.0 / Tlen,
          psd=psd,fLow=20,fMax=2000,analyticPSD_Q=True,interpolate_max=True,
          waveform_is_psi4=True)
  
    rho_1 = IP.norm(Psi4F1)
    rho_2 = IP.norm(Psi4F2)
    inner_12 = IP.ip(Psi4F1,Psi4F2)/rho_1/rho_2
    return inner_12


OBS=5
L=2
M=0

# Load DendroGR data
# Can't use str2z directly due to incompatibility with Transpose
dendro_re  = np.genfromtxt("/tmp/bssn_prof_GW_l{0}_m{1}.dat".format(L, M),
        converters={2: str2r, 3: str2r, 4: str2r, 5: str2r, 6 : str2r, 7 :
            str2r}).T
dendro_im  = np.genfromtxt("/tmp/bssn_prof_GW_l{0}_m{1}.dat".format(L, M),
        converters={2: str2i, 3: str2i, 4: str2i, 5: str2i, 6 : str2i, 7 :
            str2i}).T

de_time = dendro_re[1]
de_psi4 = dendro_re[OBS + 2]  + 1j * dendro_im[OBS +2]

# Load LazEv data        
le_time, le_real, le_imag, _, _ =\
        np.loadtxt("/tmp/pk_rad_obs_{0}_psi4_{1}_{2}.tl".format(OBS, L, M)).T
le_psi4 = -le_real - 1j * le_imag

wave1={}
wave1["time"] = le_time
wave1["psi4"] = le_psi4

wave2={}
wave2["time"] = de_time
wave2["psi4"] = de_psi4

for mass in range(10, 205, 5):
    print (mass, overlap(wave1, wave2, mass))
