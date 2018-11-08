""" A luminosity/risetime plot """

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import numpy as np
from astropy.cosmology import Planck15


fig,ax = plt.subplots(1,1,figsize=(6,5))

# AT2018gep: NEED TO UPDATE!
# This source: rise time to peak bolometric luminosity
# let's say 3 days for now, but this is a placeholder until we get the
# bolometric light curve
# the luminosity is a placeholder as well!
trise = 3
plum = 3E44
ax.scatter(
        trise, plum, marker='*', s=300, 
        facecolors='black', edgecolors='black')
ax.text(
        trise, plum, "AT2018gep", fontsize=12, 
        verticalalignment='bottom', 
        horizontalalignment='left')

# AT2018cow
# peak bol from Margutti 2018
# rise time from Perley 2018: 2-3 days
trise = 2.5
plum = 4E44
ax.errorbar(
        trise, plum, xerr=0.5, fmt='o', ms=10, 
        mfc='black', mec='black', c='k')
ax.text(
        trise, plum, "AT2018cow", fontsize=12, 
        verticalalignment='bottom', 
        horizontalalignment='left')


# iPTF16asu
# peak bol from Whitesides 2017
# this is a strict lower limit
# since they only account for the obs. flux
plum = 3.4E43
eplum = 0.3E43
trise = 3.97
etrise = 0.19
ax.errorbar(
        trise, plum, xerr=etrise, yerr=eplum,
        fmt='o', ms=10, mfc='black', mec='black', c='k')
ax.text(
        trise, plum, "iPTF16asu", fontsize=12, 
        verticalalignment='bottom', 
        horizontalalignment='left')


ax.set_ylim(1E41, 1E45)
ax.set_xlim(0.4, 100)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylabel("Peak Bolometric Luminosity $L_\mathrm{bol}$", fontsize=16)
ax.set_xlabel(
        r"Days from Explosion to Peak $L_\mathrm{bol}$", fontsize=16)
ax.yaxis.set_tick_params(labelsize=14)
ax.xaxis.set_tick_params(labelsize=14)
ax.legend(loc='upper right', fontsize=12)

plt.tight_layout()
plt.savefig("lum_rise.png")
#plt.show()
