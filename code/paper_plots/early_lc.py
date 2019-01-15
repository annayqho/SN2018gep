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

fig,ax = plt.subplots(1,1,figsize=(8,5))

f = DATA_DIR + "/ZTF18abukavn_opt_phot.dat"
dat = np.loadtxt(f, dtype=str, delimiter=' ')
instr = dat[:,0]
jd = dat[:,1].astype(float)
filt = dat[:,2]
mag = dat[:,3].astype(float)
emag = dat[:,4].astype(float)

det = np.logical_and(mag<99, ~np.isnan(mag))
nondet = np.logical_or(mag==99, np.isnan(mag))
zp = 2458370.6473
dt = jd-zp

rcol = 'k'
ucol = '#f98e09'

band = filt=='u'
ext = 0.045
choose = np.logical_and(det, band)
ax.errorbar(
        dt[choose], mag[choose]-ext, emag[choose], fmt='o', ms=5,
        mec=ucol, mfc=ucol, c=ucol, label='u', zorder=8)

band = filt=='g'
ext = 0.035
choose = np.logical_and(det, band)
ax.errorbar(
        dt[choose], mag[choose]-ext, emag[choose], fmt='o', ms=5,
        mec='#57106e', mfc='white', c='#57106e', label='g')

band = filt=='r'
ext = 0.024
choose = np.logical_and(det, band)
ax.errorbar(
        dt[choose], mag[choose]-ext, emag[choose], fmt='s', ms=6,
        mec=rcol, mfc=rcol, c=rcol, label='r', zorder=9)
choose = np.logical_and(nondet, band)
# value of limiting mag: 20.47
ax.arrow(
        2458370.6408-zp, 19.97, 0, 0.5, length_includes_head=True,
        head_width=0.01, head_length=0.1, fc='k', ec='k')

band = filt=='i'
ext = 0.018
choose = np.logical_and(det, band)
ax.errorbar(
        dt[choose], mag[choose]-ext, emag[choose], fmt='v', ms=5,
        c='grey', label='i')

band = filt=='z'
ext = 0.013
choose = np.logical_and(det, band)
ax.errorbar(
        dt[choose], mag[choose]-ext, emag[choose], fmt='s', ms=6,
        mec='k', mfc='white', label='z', c='k', zorder=1)

# UVOT light curves

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
ax.legend(loc='upper left', fontsize=12)
#ax.set_xscale('log')
ax.set_xlim(-0.1, 3)
ax.invert_yaxis()
ax2.invert_yaxis()

plt.tight_layout()
#plt.savefig("lc.png")
plt.show()
