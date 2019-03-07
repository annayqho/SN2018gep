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
# JD of first (r-band) detection
t0 = 2458370.6634
dt = jd-t0

# Initialize the figure
fig,axarr = plt.subplots(2,1,figsize=(4,6))

# Convert from mag to flux
# AB flux zero point for g-band is 3631 Jy or 1E-23 erg/cm2/s/Hz
flux = 10**(-(mag + 48.6)/2.5)
lum = 4 * np.pi * d**2 * flux
# Uncertainty in magnitude is roughly the fractional uncertainty on the flux
eflux = emag*flux
elum = 4 * np.pi * d**2 * eflux

ax = axarr[0]
# Plot the zoom-in
gband = np.logical_and(instr=='P48+ZTF', filt=='g')
choose = np.logical_and(det, gband)
ax.errorbar(
        dt[choose]*24*60, lum[choose]/1E28, yerr=elum[choose]/1E28, 
        c='k', ms=5, fmt='s', label="P48 $g$", zorder=2)

# Plot the r-band detections
rband = np.logical_and(instr=='P48+ZTF', filt=='r')
choose = np.logical_and(det, rband)
ax.errorbar(
        dt[choose]*24*60, lum[choose]/1E28, yerr=elum[choose]/1E28, 
        ms=5, fmt='o', mfc='white', mec='grey', label="P48 $r$", c='grey',
        zorder=0)

# Show the last non-detection
#ax.axvline(x=-32.5, ls=':', c='grey')
# ax.text(-32.5, 0.3, 'ND', fontsize=14,
#         horizontalalignment='center',
#         verticalalignment='center')
#ax.axvline(x=-23, ls='--', c='k')
# ax.text(-18, 0.3, '$t_0$', fontsize=16, horizontalalignment='center',
#         verticalalignment='center')

# Format this box
ax.set_xlim(-25, 80)
ax.set_ylim(0, 0.35)
ax.set_ylabel(r"$L_\nu$ [$10^{28}$ erg/s/Hz]", fontsize=16)
ax.set_xlabel("Minutes since first detection", fontsize=16)
ax.yaxis.set_tick_params(labelsize=14)
ax.xaxis.set_tick_params(labelsize=14)
ax.legend(loc='lower right', fontsize=14)

ax = axarr[1]
# Plot the full light curve
gband = np.logical_and(instr=='P48+ZTF', filt=='g')
choose = np.logical_and(det, gband)
ax.errorbar(
        dt[choose], lum[choose]/1E28, yerr=elum[choose]/1E28, 
        c='k', ms=5, fmt='s', label="P48 $g$", zorder=2)

# Plot the r-band non-detection
#ax.scatter(2458370.6408-t0, (10**(-(20.47+48.6)/2.5))/1E28, marker='.', c='r')
# ax.arrow((2458370.6408-t0), 0.2,
#         0, -0.2, length_includes_head=True,
#         head_width=0.1, head_length=0.03, fc='k', ec='k')

# Fit a polynomial (quadratic)
sec = dt[choose] < 3
out,cov = np.polyfit(
        dt[choose][sec], lum[choose][sec], 
        deg=2, w=1/elum[choose][sec], cov=True)
xlab = np.linspace(-1,2)
ylab = out[0]*xlab**2 + out[1]*xlab + out[2]
ax.plot(xlab, ylab/1E28, c='k', ls='--')
axarr[0].plot(xlab*24*60, ylab/1E28, c='k', ls='--')

# Plot the r-band detections
rband = np.logical_and(instr=='P48+ZTF', filt=='r')
choose = np.logical_and(det, rband)
flux = 10**(-(mag[choose] + 48.6)/2.5)
lum = 4 * np.pi * d**2 * flux
dt = jd[choose]-t0
eflux = emag[choose]*flux
elum = 4 * np.pi * d**2 * eflux
ax.errorbar(
        dt, lum/1E28, yerr=elum/1E28, 
        ms=5, fmt='o', mfc='white', mec='grey', label="P48 $r$", c='grey',
        zorder=0)

# Vertical line for the first UVOT epoch
textor = 'vertical' # textorientation
ax.axvline(x=0.48, lw=2, c='lightblue') # UVOT
ax.text(
        0.48, 4.2, 'Swift', fontsize=14, horizontalalignment='right',
        rotation=textor)
ax.axvline(x=0.7, lw=2, c='lightblue') # LT
ax.text(
        0.7, 4.2, 'LT', fontsize=14, horizontalalignment='right',
        rotation=textor)
ax.axvline(x=1.0, lw=2, c='lightblue') # P200/P60
ax.text(
        1.0, 0.8, 'P200', fontsize=14, horizontalalignment='right',
        rotation=textor)
ax.text(1.0, 4.2, 'P60', fontsize=14, horizontalalignment='right',
        rotation=textor)
ax.axvline(x=1.7, lw=2, c='lightblue') # LT
ax.text(1.7, 4.2, 'LT', fontsize=14, horizontalalignment='right',
        rotation=textor)
ax.axvline(x=2.0, lw=2, c='lightblue') # P200
ax.text(2.0, 2.0, 'P200', fontsize=14, horizontalalignment='right',
        rotation=textor)
ax.axvline(x=2.7, lw=2, c='lightblue') # LT
ax.text(2.7, 4.2, 'LT', fontsize=14, horizontalalignment='right',
        rotation=textor)

ax.legend(loc='lower right', fontsize=14)
ax.set_ylabel(r"$L_\nu$ [$10^{28}$ erg/s/Hz]", fontsize=16)
ax.set_xlabel("Days since first detection", fontsize=16)
ax.yaxis.set_tick_params(labelsize=14)
ax.xaxis.set_tick_params(labelsize=14)
ax.set_xlim(-0.2, 2.2)
ax.set_ylim(-0.5,4.7)

# Print fitting parameters
a = out[0]
ea = np.sqrt(cov[0][0])
b = out[1]
eb = np.sqrt(cov[1][1])
c = out[2]
ec = np.sqrt(cov[2][2])

plt.tight_layout()
plt.savefig("early_data.png")
#plt.show()
