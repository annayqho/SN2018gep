""" Christoffer's host-subtracted photometry
Plot the light curve! """

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15
import glob


DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/uv"
d = Planck15.luminosity_distance(z=0.033).cgs.value

fig,ax = plt.subplots(1,2,figsize=(8,5), sharey=True, sharex=True)

f = DATA_DIR + "/UVOT_lightcurve_maghist_full_hostsub.ascii"
dat = np.loadtxt(f, dtype=str, delimiter=' ')
jd = dat[:,0].astype(float)
filt = dat[:,1]
fnu_mjy = dat[:,2].astype(float)
fact = 1E-3 * 1E-23 * 4 * np.pi * d**2 / 1E28
lnu = fnu_mjy * fact 
elnu = dat[:,3].astype(float) * fact 

zp = 2458370.6473
dt = jd-zp

choose = filt=='UVM2'
ax[0].errorbar(
        dt[choose], lnu[choose], elnu[choose], fmt='s', ms=6,
        mec='k', mfc='k', c='k', label='UVM2', zorder=9)

choose = filt=='UVW2'
ax[1].errorbar(
        dt[choose], lnu[choose], elnu[choose], fmt='s', ms=6,
        mec='k', mfc='k', c='k', label='UVW2', zorder=9)

choose = filt=='U'
ax[0].errorbar(
        dt[choose], lnu[choose], elnu[choose], fmt='o', ms=5,
        mec='#f98e09', mfc='#f98e09', c='#f98e09', label='U', zorder=8)

choose = filt=='UVW1'
ax[1].errorbar(
        dt[choose], lnu[choose], elnu[choose], fmt='o', ms=5,
        mec='#f98e09', mfc='#f98e09', c='#f98e09', label='UVW1')
 
choose = filt=='B'
ax[1].errorbar(
        dt[choose], lnu[choose], elnu[choose], fmt='v', ms=7,
        mec='#57106e', mfc='white', c='#57106e', label='B')

choose = filt=='V'
ax[0].errorbar(
        dt[choose], lnu[choose], elnu[choose], fmt='v', ms=7,
        mec='#57106e', mfc='white', c='#57106e', label='V')

ax[0].set_ylabel(r"UVOT $L_\nu$ ($10^{28}$ erg/s/Hz)", fontsize=16)
for ii in np.arange(2):
    ax[ii].yaxis.set_tick_params(labelsize=14)
    ax[ii].xaxis.set_tick_params(labelsize=14)
    ax[ii].legend(loc='upper right', fontsize=12)
ax[0].set_xlabel("Days since JD 2458370.6473", fontsize=16)
ax[1].set_xlabel("(UT 2018 Sept 09.15)", fontsize=16)

plt.subplots_adjust(wspace=0)
#plt.tight_layout()
plt.savefig("uv_lc.png")
#plt.show()
