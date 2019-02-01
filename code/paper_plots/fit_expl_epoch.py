""" 
Fit for the explosion epoch by fitting a polynomial 
to the g-band light curve """

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.cosmology import Planck15

d = Planck15.luminosity_distance(z=0.033).cgs.value

DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/phot"
f = DATA_DIR + "/ZTF18abukavn_opt_phot.dat"
dat = np.loadtxt(f, dtype=str, delimiter=' ')
instr = dat[:,0]
jd = dat[:,1].astype(float)
filt = dat[:,2]
mag = dat[:,3].astype(float)
emag = dat[:,4].astype(float)

det = np.logical_and(mag<99, ~np.isnan(mag))
tofit = np.logical_and(instr=='P48+ZTF', filt=='g')
choose = np.logical_and(det, tofit)
# JD of first (r-band) detection
t0 = 2458370.6634
dt = jd[choose]-t0

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
ax.errorbar(
        dt, lum/1E28, yerr=elum/1E28, c='k', ms=10, fmt='.', label="$g$-band")

# Plot the r-band non-detection
#ax.scatter(2458370.6408-t0, (10**(-(20.47+48.6)/2.5))/1E28, marker='.', c='r')
ax.arrow((2458370.6408-t0)*24, 0.2,
        0, -0.2, length_includes_head=True,
        head_width=0.1, head_length=0.03, fc='k', ec='k')

# Inset axis
axins = inset_axes(
        ax, 2, 1, loc=1,
        bbox_to_anchor=(0.98,0.98),
        bbox_transform=ax.transAxes)
axins.arrow((2458370.6408-t0)*24, 0.2,
        0, -0.2, length_includes_head=True,
        head_width=0.1, head_length=0.03, fc='k', ec='k')
# axins.scatter(
#         (2458370.6408-t0)*24, (10**(-(20.47+48.6)/2.5))/1E28,
#         marker='.', c='r')
axins.errorbar(dt*24, lum/1E28, yerr=elum/1E28, c='k', ms=10, fmt='.')
axins.axhline(y=0, c='k', lw=0.5)
axins.set_xlim(-0.05*24,0.05*24)
axins.set_ylim(-0.1,0.3)

# Fit a polynomial (quadratic)
sec = dt < 3

out,cov = np.polyfit(dt[sec], lum[sec], deg=2, w=1/elum[sec], cov=True)
xlab = np.linspace(-1,2)
ylab = out[0]*xlab**2 + out[1]*xlab + out[2]
ax.plot(xlab, ylab/1E28, c='k', ls='--')
axins.plot(xlab*24, ylab/1E28, c='k', ls='--')

ax.set_ylabel(r"$L_\nu$ [$10^{28}$ erg/s/Hz]", fontsize=16)
ax.set_xlabel("Days since first (r-band) detection", fontsize=16)
axins.set_xlabel("Hours since first $r$-band det.", fontsize=14)
ax.yaxis.set_tick_params(labelsize=14)
axins.yaxis.set_tick_params(labelsize=12)
ax.xaxis.set_tick_params(labelsize=14)
axins.xaxis.set_tick_params(labelsize=12)
ax.set_xlim(-0.5, 6)
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
