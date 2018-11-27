""" A luminosity/risetime plot """

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import numpy as np
from astropy.cosmology import Planck15


def at2018gep():
    trise = 1391.0399966526045 / 86400 # convert from seconds to days
    plum = 3E44
    ax.scatter(
            trise, plum, marker='*', s=300, 
            facecolors='black', edgecolors='black')
    ax.text(
            trise*1.1, plum, "AT2018gep", fontsize=14, 
            verticalalignment='bottom', 
            horizontalalignment='left')
    ax.arrow(
            trise, plum, -0.006, 0, length_includes_head=True,
            head_width=plum/5, head_length=0.001, fc='k')


def at2018cow():
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

fig,ax = plt.subplots(1,1,figsize=(6,5))

at2018gep()
at2018cow()

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

# KN (AT2017gfo)
plum = 1E42
trise = 0.5
ax.scatter(
        trise, plum, 
        marker='o', s=100, facecolor='black', edgecolor='black')
ax.text(
        trise, plum, "AT2017gfo", fontsize=12, 
        verticalalignment='bottom', 
        horizontalalignment='left')



ax.set_ylim(1E41, 1E45)
ax.set_xlim(0.005, 100)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylabel("Peak Bolometric Luminosity $L_\mathrm{bol}$", fontsize=16)
ax.set_xlabel(
        r"Days from Explosion to Peak $L_\mathrm{bol}$", fontsize=16)
ax.yaxis.set_tick_params(labelsize=14)
ax.xaxis.set_tick_params(labelsize=14)
ax.legend(loc='upper right', fontsize=12)

plt.tight_layout()
#plt.savefig("lum_rise.png")
plt.show()
