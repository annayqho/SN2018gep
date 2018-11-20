""" Christoffer's host-subtracted photometry
Plot the light curve! """

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
mjd = dat[:,1].astype(float)
filt = dat[:,2]
mag = dat[:,3].astype(float)
emag = dat[:,4].astype(float)

det = np.logical_and(mag<99, ~np.isnan(mag))
nondet = np.logical_or(mag==99, np.isnan(mag))
zp = 2458370.6473
dt = mjd-zp

rcol = 'k'
ucol = '#f98e09'

band = filt=='u'
choose = np.logical_and(det, band)
ax.errorbar(
        dt[choose], mag[choose], emag[choose], fmt='o', ms=5,
        mec=ucol, mfc=ucol, c=ucol, label='u', zorder=8)

band = filt=='g'
choose = np.logical_and(det, band)
ax.errorbar(
        dt[choose], mag[choose], emag[choose], fmt='o', ms=5,
        mec='#57106e', mfc='white', c='#57106e', label='g')

band = filt=='r'
choose = np.logical_and(det, band)
ax.errorbar(
        dt[choose], mag[choose], emag[choose], fmt='s', ms=6,
        mec=rcol, mfc=rcol, c=rcol, label='r', zorder=9)
choose = np.logical_and(nondet, band)
# value of limiting mag: 20.47
ax.arrow(
        2458370.6408-zp, 19.97, 0, 0.5, length_includes_head=True,
        head_width=1, head_length=0.1, fc='k', ec='k')

# zoomed-in window showing the earliest non-detection and detection
axins = inset_axes(
        ax, 2, 1, loc=1,
        bbox_to_anchor=(0.87,0.98),
        bbox_transform=ax.transAxes)
choose = np.logical_and(det, band)
axins.errorbar(
    dt[choose]*24, mag[choose], emag[choose], fmt='s', ms=6,
    mec=rcol, mfc=rcol, c=rcol, label='r', zorder=9)
choose = np.logical_and(nondet, band)
axins.arrow(
        2458370.6408-zp, 19.97, 0, 0.5, length_includes_head=True,
        head_width=0.2, head_length=0.3, fc='k', ec='k')
band = filt=='g'
choose = np.logical_and(det, band)
axins.errorbar(
        dt[choose]*24, mag[choose], emag[choose], fmt='o', ms=5,
        mec='#57106e', mfc='white', c='#57106e', label='g')
axins.set_xlim(-0.3,3)
axins.set_ylim(18,21)
axins.tick_params(axis='both', labelsize=12)
axins.set_xlabel(r"Hours since $t_0$", fontsize=12)
axins.invert_yaxis()
ax.plot([-1, -1], [21, 18], c='k', lw=0.5)
ax.plot([1, 1], [21, 18], c='k', lw=0.5)
ax.plot([-1, 1], [18, 18], c='k', lw=0.5)
ax.plot([-1, 1], [21, 21], c='k', lw=0.5)
#mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
 
band = filt=='i'
choose = np.logical_and(det, band)
ax.errorbar(
        dt[choose], mag[choose], emag[choose], fmt='v', ms=5,
        c='grey', label='i')
 
band = filt=='z'
choose = np.logical_and(det, band)
ax.errorbar(
        dt[choose], mag[choose], emag[choose], fmt='s', ms=6,
        mec='k', mfc='white', label='z', c='k', zorder=1)

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
ax.legend(loc='upper right', fontsize=12)
#ax.set_xscale('log')
ax.invert_yaxis()
ax2.invert_yaxis()

#plt.tight_layout()
plt.savefig("lc.png")
#plt.show()
