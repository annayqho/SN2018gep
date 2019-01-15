""" Plot all data from < 3 days into one figure """

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15
import glob

DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/phot"

def plot_band(x, y, ey, filt, f, col, fcol, ext, shape):
    choose = np.logical_and(filt==f, x < 1)
    ax.errorbar(
            x[choose], y[choose]-ext, ey[choose], fmt=shape, ms=5,
            mec=col, mfc=fcol, c=col, label=f, zorder=8)
    order = np.argsort(x[choose])
    ax.plot(x[choose][order], y[choose][order]-ext, c=col)


fig,ax = plt.subplots(1,1,figsize=(8,5))

f = DATA_DIR + "/ZTF18abukavn_opt_phot.dat"
dat = np.loadtxt(f, dtype=str, delimiter=' ')
instr = dat[:,0]
jd = dat[:,1].astype(float)
filts = dat[:,2]
mag = dat[:,3].astype(float)
emag = dat[:,4].astype(float)

det = np.logical_and(mag<99, ~np.isnan(mag))
nondet = np.logical_or(mag==99, np.isnan(mag))
zp = 2458370.6473
dt = jd-zp

plot_band(
        dt[det], mag[det], emag[det], filts[det], 
        'u', '#f98e09', '#f98e09', 0.045, 'o')
plot_band(
        dt[det], mag[det], emag[det], filts[det], 
        'g', '#57106e', '#57106e', 0.035, 'D')
plot_band(
        dt[det], mag[det], emag[det], filts[det], 
        'r', 'k', 'k', 0.024, 's')
plot_band(
        dt[det], mag[det], emag[det], filts[det], 
        'i', 'grey', 'grey', 0.018, 'P')
# value of limiting mag: 20.47
ax.arrow(
       2458370.6408-zp, 19.97, 0, 0.5, length_includes_head=True,
       head_width=0.01, head_length=0.1, fc='k', ec='k')
# There's no z-band in the first day

# UVOT light curves
DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/uv"
d = Planck15.luminosity_distance(z=0.033).cgs.value

f = DATA_DIR + "/UVOT_lightcurve_maghist_full_hostsub.ascii"
dat = np.loadtxt(f, dtype=str, delimiter=' ')
jd = dat[:,0].astype(float)
filts = dat[:,1]
fnu_mjy = dat[:,2].astype(float)
# convert mJy to AB magnitudes
mab = -2.5 * np.log10(fnu_mjy * 1E-3) + 8.90 
efnu = dat[:,3].astype(float) 
emab = 2.5 * (efnu / fnu_mjy)
zp = 2458370.6473
dt = jd-zp

plot_band(
        dt, mab, emab, filts, 
        'UVM2', 'k', 'white', 0, 's')
plot_band(
        dt, mab, emab, filts, 
        'UVW2', '#57106e', 'white', 0, 'o')
plot_band(
        dt, mab, emab, filts, 
        'U', '#f98e09', 'white', 0, 'D')
plot_band(
        dt, mab, emab, filts, 
        'UVW1', 'magenta', 'white', 0, '^')
plot_band(
        dt, mab, emab, filts, 
        'B', 'grey', 'white', 0, 'P')
plot_band(
        dt, mab, emab, filts, 
        'V', 'black', 'white', 0, 'x')

ax2 = ax.twinx()
ax2.set_ylabel(
        "Absolute Magnitude",
        fontsize=14, rotation=270, labelpad=15.0)
y_f = lambda y_i: y_i-Planck15.distmod(z=0.033).value
ymin, ymax = ax.get_ylim()
ax2.set_ylim((y_f(ymin), y_f(ymax)))
ax2.plot([],[])
ax2.tick_params(axis='both', labelsize=14)

ax.set_ylabel("Apparent Magnitude", fontsize=16)
ax.set_xlabel(
    r"Days since $t_0=$JD 2458370.6473 (UT 2018 Sept 09.15)", fontsize=16)
ax.yaxis.set_tick_params(labelsize=14)
ax.xaxis.set_tick_params(labelsize=14)
ax.legend(loc='lower right', fontsize=12, ncol=5)
#ax.set_xscale('log')
ax.set_xlim(-0.1, 1)
ax.set_ylim(15, 21)
ax.invert_yaxis()
ax2.invert_yaxis()

plt.tight_layout()
plt.savefig("early_lc.png")
#plt.show()
