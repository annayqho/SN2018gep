""" A luminosity/risetime plot """

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import numpy as np
from astropy.cosmology import Planck15


def sn2018gep(axarr):
    """ Rise time and peak mag, peak lbol, for SN2018gep """
    # left panel is g-band mag
    # right panel is bolometric luminosity
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


def at2018cow(axarr):
    """ Rise time and peak mag, peak lbol, for AT2018cow """
    # rise time will be upper limit in both cases 
    trise = np.array([2, 2.5])
    plum = np.array([-20.5, 3.5E44])
    eplum = np.array([0.03, 3.5E44*(6E9/8.5E10)])
    apos = trise/3
    hw = np.array([0.02, 1E44])

    for ii,ax in enumerate(axarr):
        ax.errorbar(
                trise[ii], plum[ii], yerr=eplum[ii], fmt='o', ms=10, c='k')
        ax.text(
                trise[ii], plum[ii], "AT2018cow", fontsize=12, 
                verticalalignment='bottom', 
                horizontalalignment='left')
        ax.arrow(
                trise[ii], plum[ii], -apos[ii], 0, length_includes_head=True,
                head_width=hw[ii], head_length=apos[ii]/3, fc='k')


def iptf16asu(axarr):
    """ 
    rise time, peak Mg, peak Lbol for iPTF16asu
    """
    # rise time to Mg is a lower limit because we only observe it
    # rising by 1 mag, not by 2 mag
    # we do resolve the bolometric rise
    trise = [1.71, 1.4]

    # we do resolve Mg, but Lbol is a strict lower limit
    plum = [-20.4, 3.4E43]
    eplum = [0.09, 0.3E43]

    for ii, ax in enumerate(axarr):
        ax.errorbar(
                trise[ii], plum[ii], yerr=eplum[ii],
                fmt='o', ms=10, mfc='black', mec='black', c='k')
        ax.text(
                trise[ii], plum[ii], "iPTF16asu", fontsize=12, 
                verticalalignment='bottom', 
                horizontalalignment='left')
    axarr[0].arrow(
            trise[0], plum[0], -0.8, 0, length_includes_head=True,
            head_width=0.02, head_length=0.2, fc='k')
    axarr[1].arrow(
            trise[1], plum[1], 0, 5E43, length_includes_head=True,
            head_width=0.4, head_length=1E44/5, fc='k')


fig,axarr = plt.subplots(1,2,figsize=(10,5), sharex=True)
sn2018gep(axarr)
at2018cow(axarr)
iptf16asu(axarr)

axarr[0].set_ylabel("Observed $M_g$", fontsize=16)
axarr[0].invert_yaxis()
axarr[0].set_xlabel(
    r"Observed $t_2$ [days]", fontsize=16)
axarr[1].set_xlabel(
    r"Observed $t_{1/2}$ [days]", fontsize=16)

ax = axarr[1]
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
