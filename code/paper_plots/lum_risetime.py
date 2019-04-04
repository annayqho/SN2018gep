""" A luminosity/risetime plot """

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import numpy as np
from astropy.cosmology import Planck15


def sn2018gep(axarr):
    trise = [3, 0.2]
    plum = [-20, 3E44]

    # left panel: rise to max g-band 
    for ii,ax in enumerate(axarr):
        ax.scatter(
                trise[ii], plum[ii], marker='*', s=300, 
                facecolors='black', edgecolors='black')
        ax.text(
                trise[ii]*1.1, plum[ii], "SN2018gep", fontsize=14, 
                verticalalignment='bottom', 
                horizontalalignment='left')


def at2018cow(ax):
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


fig,axarr = plt.subplots(1,2,figsize=(10,5), sharex=True)
sn2018gep(axarr)


axarr[0].set_ylabel("$M_g$", fontsize=16)
axarr[0].set_xlabel(
    r"$t_2$ [days]", fontsize=16)

# Right panel: peak bolometric luminosity
axarr[0].set_xlabel(
    r"$t_{1/2}$ [days]", fontsize=16)


#at2018cow(ax)

ax = axarr[1]
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
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylabel("$L_\mathrm{bol}$", fontsize=16)
ax.legend(loc='upper right', fontsize=12)



for ax in axarr:
    ax.set_xlim(0.05, 100)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.yaxis.set_tick_params(labelsize=14)

fig.tight_layout()
#plt.savefig("lum_rise.png")
plt.show()
