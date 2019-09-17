""" A plot of rise time (Rest frame days)
vs observer-frame peak abs mag in g-band 
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
textsize=14


def sn2018gep(ax):
    """ Rise time and peak g-band mag for SN2018gep """
    trise = 3
    plum = -19.9

    ax.scatter(
            trise, plum, marker='*', s=500, 
            facecolors='black', edgecolors='black')
    ax.text(trise*1.1, plum/1.005, "18gep (Ic-BL)", fontsize=textsize)


def at2018cow(ax):
    """ Rise time and peak mag """
    x = 1.5
    y = -20.4
    ax.errorbar(
            x, y, xerr=0.5, yerr=0.03, c='k', marker='o')
    ax.text(x, y*1.005, "18cow", fontsize=textsize)
    

def iptf16asu(ax):
    """ 
    rise time, peak Mg

    I calculated this by taking the quadratic fit to the early
    light curve in Whitesides+2017, and their estimate for the peak time.
    I calculated the t_1/2, rise using that quadratic function.
    I took their measurement of the g-band max magnitude.
    """
    # from the paper: 3.97 \pm 0.19 days in g-band
    # peak is -20.4
    trise = 2.7

    # we do resolve Mg, but Lbol is a strict lower limit
    plum = -20.4

    ax.errorbar(
            trise, plum, yerr=0.1, xerr=0.19, marker='o', c='k')
    ax.text(
            trise, plum*1.005, "16asu (Ic-BL)", fontsize=textsize)


def rest2018(ax):
    """
    plot the data from Rest (2018)
    """
    trise = 1.4
    plum = -18.8
    ax.scatter(
            trise, plum, marker='o', c='k')
    ax.text(
            trise, plum*1.005, "15K", fontsize=textsize)


def drout(ax):
    """ using rest-frame g-band """
    x = 1.0
    y = -17.5
    ax.errorbar(
           x, y, xerr=0.1, yerr=0.11, c='k', marker='o')
    ax.text(x, y*1.01, "10ah", fontsize=textsize)

    x = 1.0
    y = -18.2
    ax.errorbar(
           x, y, xerr=0.1, yerr=0.11, c='k', marker='o')
    ax.text(x, y/1.01, "10bjp", fontsize=textsize)
    
    x = 2.9
    y = -19.5
    ax.errorbar(
           x, y, xerr=0.1, yerr=0.08, c='k', marker='o')
    ax.text(x, y*1.005, "11qr", fontsize=textsize)
    
    x = 2.2
    y = -19.1
    ax.errorbar(
           x, y, yerr=0.1, c='k', marker='o')
    ax.arrow(
           x, y, -0.3, 0, color='k', head_width=0.1, head_length=0.1)
    ax.text(x, y*1.005, "12bv", fontsize=textsize)
    
    x = 1.0
    y = -18.3
    ax.errorbar(
           x, y, yerr=0.1, c='k', marker='o')
    ax.arrow(
           x, y, -0.3, 0, color='k', head_width=0.1, head_length=0.1)
    ax.text(x, y*1.005, "12brf", fontsize=textsize)
    

def arcavi(ax):
    """ using rest-frame g-band """
    x = 1.12
    y = -20.7
    ax.errorbar(
           x, y, yerr=0.1, c='k', marker='o')
    ax.arrow(
           x, y, -0.3, 0, color='k', head_width=0.1, head_length=0.1)
    ax.text(x, y*1.005, "04D4ec", fontsize=textsize)

    x = 3
    y = -20.08
    ax.errorbar(
           x, y, yerr=0.1, c='k', marker='o')
    ax.text(x, y*1.005, "05D2bk", fontsize=textsize)

    x = 4.6
    y = -20.7
    ax.errorbar(
           x, y, yerr=0.1, c='k', marker='o')
    ax.text(x, y*1.005, "06D1hc", fontsize=textsize)


fig,ax = plt.subplots(1,1,figsize=(5,5))
sn2018gep(ax)
at2018cow(ax)
iptf16asu(ax)
drout(ax)
arcavi(ax)
rest2018(ax)

ax.set_ylabel("Peak Mag (Rest-frame $g$-band)", fontsize=16)
ax.set_xlim(0,6)
ax.set_ylim(-21, -17)
#ax.set_xscale('log')
ax.invert_yaxis()
ax.set_xlabel(
    r"$t_{\mathrm{1/2,rise}}$ (rest frame days)", fontsize=16)

ax.xaxis.set_tick_params(labelsize=14)
ax.yaxis.set_tick_params(labelsize=14)

fig.tight_layout()
plt.savefig("lum_rise.eps", format='eps', dpi=1000)

#plt.show()
