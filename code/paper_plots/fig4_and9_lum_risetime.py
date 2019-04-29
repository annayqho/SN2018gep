""" A luminosity/risetime plot 

Take the values from Rest (2018) and Margutti (2019) and add a few events
One panel is Fig 4, the other is Fig 9 in the paper
"""

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import numpy as np
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/code")
from load_lc import get_uv_lc, get_lc
from load_lum import load_lc
from astropy.cosmology import Planck15
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

datadir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data"


def sn2018gep(axarr):
    """ Rise time and peak mag, peak lbol, for SN2018gep """
    # left panel is g-band mag
    # right panel is bolometric luminosity
    trise = [3, 2]
    plum = [-20, 3E44]

    # left panel: rise to max g-band 
    for ii,ax in enumerate(axarr):
        ax.scatter(
                trise[ii], plum[ii], marker='*', s=300, 
                facecolors='black', edgecolors='black')


def at2018cow(axarr):
    """ Rise time and peak mag, peak lbol, for AT2018cow """
    axarr[0].errorbar(
            2.5, -20.15, yerr=0.25, c='k', marker='o')
    

def iptf16asu(axarr):
    """ 
    rise time, peak Mg, peak Lbol for iPTF16asu
    """
    # from the paper: 3.97 \pm 0.19 days in g-band
    # peak is -20.4
    trise = [3.97, 1.4]
    etrise = [0.19, 0.19] # ignore the second number

    # we do resolve Mg, but Lbol is a strict lower limit
    plum = [-20.4, 3.4E43]
    eplum = [0.09, 0.3E43]

    axarr[0].errorbar(
            trise[0], plum[0], marker='o', c='k',
            yerr=eplum[0], xerr=etrise[0])


def rest2018(axarr):
    """
    plot the data from Rest (2018)
    """
    dat = np.loadtxt(datadir + "/rest2018.csv", delimiter=',')
    axarr[0].scatter(dat[:,0], dat[:,1], c='k')


def margutti2018(axarr):
    """
    plot the rise time vs. bol lum from Margutti (2018)
    """
    dat = np.loadtxt(datadir + "/margutti_fbots.txt", delimiter=',')
    axarr[1].scatter(dat[:,0], dat[:,1], c='k', label="FBOT")
    dat = np.loadtxt(datadir + "/margutti_slsn.txt", delimiter=',')
    axarr[1].scatter(dat[:,0], dat[:,1], edgecolor='k', facecolor='white', marker='s', label="SLSN")
    dat = np.loadtxt(datadir + "/margutti_IIn.txt", delimiter=',')
    axarr[1].scatter(dat[:,0], dat[:,1], marker='x', c='k', label="IIn")
    dat = np.loadtxt(datadir + "/margutti_IIL.txt", delimiter=',')
    axarr[1].scatter(dat[:,0], dat[:,1], marker='x', c='#e55c30', label="IIL")
    dat = np.loadtxt(datadir + "/margutti_IIP.txt", delimiter=',')
    axarr[1].scatter(dat[:,0], dat[:,1], marker='x', c='#84206b', label="IIP")
    dat = np.loadtxt(datadir + "/margutti_Ibc.txt", delimiter=',')
    axarr[1].scatter(
            dat[:,0], dat[:,1], marker='D', edgecolor='#84206b', facecolor='white', label="Ibc")


fig,axarr = plt.subplots(1,2,figsize=(10,5))
#sn2018gep(axarr)
at2018cow(axarr)
iptf16asu(axarr)
rest2018(axarr)
margutti2018(axarr)

ax = axarr[0]
ax.set_ylabel("Peak Abs Mag (Optical Filter)", fontsize=16)
ax.set_xlim(1,20)
ax.set_ylim(-21, -16)
ax.set_xscale('log')
ax.invert_yaxis()
ax.set_xlabel(
    r"Rise Time (rest frame days)", fontsize=16)

ax = axarr[1]
ax.set_xlim(0.48,100)
ax.set_ylim(1E41, 1E45)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylabel("Peak Bolometric Luminosity (erg/s)", fontsize=16)
#ax.legend(loc='lower left', fontsize=12)
ax.set_xlabel(
    r"Rise Time (rest frame days)", fontsize=16)

for ax in axarr:
    ax.xaxis.set_tick_params(labelsize=14)
    ax.yaxis.set_tick_params(labelsize=14)

fig.tight_layout()
#plt.savefig("lum_rise_wogep.eps", format='eps', dpi=1000)

plt.show()
