""" From the g & r light curves at very early times,
estimate the blackbody fit. """

import numpy as np
import matplotlib.pyplot as plt
from load_lc import get_lc

dt, filt, mag, emag = get_lc()
tgrid = 0.05 # time to interpolate onto

# Get the r-band color
choose = np.logical_and(filt=='r', emag < 99)
r = np.interp(tgrid, dt[choose], mag[choose])

# Get the g-band color
choose = np.logical_and(filt=='g', emag < 99)
g = np.interp(tgrid, dt[choose], mag[choose])

# g is brighter than r. So the lower limit is that
# the temperature corresponds to the peak wavelength
# being at the g-band filter.
# Using the filter profile service 
# (http://svo2.cab.inta-csic.es/svo/theory/fps/) I get

# lambda_eff
gwav = 4722.7
gwav_cm = gwav*1E-8
rwav = 6339.6 
rwav_cm = rwav*1E-8

# Wien's law: 
A = 0.002897755 # m K
Tlimit = A / (gwav/10**(10))
# so the minimum temperature is T = 6136 Kelvin

# use the mag ratio to get a flux ratio
flux_ratio = 10**((g-r)/(-2.5)) # g/r

# use the flux ratio to solve for the temperature
# B_lambda(T) = a/lambda^5 / (e^(b/(lambda*T))-1)
h = 6.63E-27
c = 3E10
k = 1.38E-16
b = h*c/k

const = flux_ratio * (gwav/rwav)**5 
temp = np.linspace(Tlimit,Tlimit*2,10000)
func = (np.exp(b/(rwav_cm*temp))-1)/(np.exp(b/(gwav_cm*temp))-1)-const

teff = np.interp(0,func,temp)
# 6928 Kelvin
