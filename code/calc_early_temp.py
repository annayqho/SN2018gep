""" From the g & r light curves at very early times,
estimate the blackbody fit. """

import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15
from load_lc import get_lc
from astropy.modeling.blackbody import blackbody_lambda, blackbody_nu

d = Planck15.luminosity_distance(z=0.03154).cgs.value

dt, filt, mag, emag = get_lc()
tgrid = 0.05 # time to interpolate onto

# Get the r-band color
choose = np.logical_and(filt=='r', emag < 99)
npts = len(dt[choose])
ndraw = 10000
# make 10000 new versions of the magnitude array,
# repopulated with the magnitude from the Gaussian dist
# with STD of that uncertainty
mag_new = np.zeros((npts, ndraw))
for ii, magval in enumerate(mag[choose]):
    mag_new[ii,:] = np.random.normal(
            loc=magval, scale=emag[choose][ii], size=ndraw)
rs = np.zeros(ndraw)
for ii in np.arange(ndraw):
    rs[ii] = np.interp(tgrid, dt[choose], mag_new[:,ii])
rmean = np.mean(rs)
ermean = np.std(rs)

# Get the g-band color
choose = np.logical_and(filt=='g', emag < 99)
npts = len(dt[choose])
ndraw = 10000
# make 10000 new versions of the magnitude array,
# repopulated with the magnitude from the Gaussian dist
# with STD of that uncertainty
mag_new = np.zeros((npts, ndraw))
for ii, magval in enumerate(mag[choose]):
    mag_new[ii,:] = np.random.normal(
            loc=magval, scale=emag[choose][ii], size=ndraw)
gs = np.zeros(ndraw)
for ii in np.arange(ndraw):
    gs[ii] = np.interp(tgrid, dt[choose], mag_new[:,ii])
gmean = np.mean(gs)
egmean = np.std(gs)

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

# Generate a flux density given a temperature and a radius
T = 7000
R = 1E15
c = 3E10
wl = np.array([gwav, rwav])
bb = blackbody_lambda(wl, 5000).value # erg/AA/cm2/s/sr
flam = (R/d)**2 * np.pi * bb
wl_cm = wl*1E-8
fnu = (wl_cm**2/c)*flam
modmag = -2.5*np.log10(fnu)-48.60

plt.scatter(wl, modmag)

# # use the mag ratio to get a flux ratio
# flux_ratio = 10**((g-r)/(-2.5)) # g/r
# 
# # use the flux ratio to solve for the temperature
# # B_lambda(T) = a/lambda^5 / (e^(b/(lambda*T))-1)
# h = 6.63E-27
# c = 3E10
# k = 1.38E-16
# b = h*c/k
# 
# const = flux_ratio * (gwav/rwav)**5 
# temp = np.linspace(Tlimit,Tlimit*2,10000)
# func = (np.exp(b/(rwav_cm*temp))-1)/(np.exp(b/(gwav_cm*temp))-1)-const
# 
# teff = np.interp(0,func,temp)
# # 6928 Kelvin

plt.show()
