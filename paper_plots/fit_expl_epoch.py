""" 
Fit for the explosion epoch by fitting a polynomial 
to the g-band light curve """

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.cosmology import Planck15

d = Planck15.luminosity_distance(z=0.033).cgs.value

DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/phot"
f = DATA_DIR + "/ZTF18abukavn_opt_phot.dat"
dat = np.loadtxt(f, dtype=str, delimiter=' ')
instr = dat[:,0]
mjd = dat[:,1].astype(float)
filt = dat[:,2]
mag = dat[:,3].astype(float)
emag = dat[:,4].astype(float)

det = np.logical_and(mag<99, ~np.isnan(mag))
tofit = np.logical_and(instr=='P48+ZTF', filt=='g')
choose = np.logical_and(det, tofit)
dt = mjd[choose]-mjd[choose][0]

# Initialize the figure
fig,ax = plt.subplots(1,1,figsize=(5,3))

# Convert from mag to flux
# AB flux zero point for g-band is 3631 Jy or 1E-23 erg/cm2/s/Hz
flux = 10**(-(mag[choose] + 48.6)/2.5)
lum = 4 * np.pi * d**2 * flux
# Uncertainty in magnitude is roughly the fractional uncertainty on the flux
eflux = emag[choose]*flux
elum = 4 * np.pi * d**2 * eflux

# Plot the light curve
ax.errorbar(dt, lum/1E28, yerr=elum/1E28, c='k', ms=5, fmt='.')

# Fit a polynomial (quadratic)
sec = dt < 3

out,cov = np.polyfit(dt[sec], lum[sec], deg=2, w=1/elum[sec], cov=True)
xlab = np.linspace(-1,2)
ylab = out[0]*xlab**2 + out[1]*xlab + out[2]
plt.plot(xlab, ylab/1E28, c='k', ls='--')

ax.set_ylabel(r"$L_\nu$ [$10^{28}$ erg/s/Hz]", fontsize=16)
ax.set_xlabel("Days since JD 2458370.6634", fontsize=16)
ax.yaxis.set_tick_params(labelsize=14)
ax.xaxis.set_tick_params(labelsize=14)
ax.set_xlim(-3, 8)
ax.set_ylim(-1,5)

# Print fitting parameters
a = out[0]
ea = np.sqrt(cov[0][0])
b = out[1]
eb = np.sqrt(cov[1][1])
c = out[2]
ec = np.sqrt(cov[2][2])

plt.tight_layout()
plt.savefig("quadfit.png")
#plt.show()
