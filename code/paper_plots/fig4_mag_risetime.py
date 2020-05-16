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
            trise, plum, marker='o', c='red')
    ax.text(trise*1.02, plum, "18gep", fontsize=textsize, color='red')


def at2018cow(ax):
    """ Rise time and peak mag """
    x = 1.43
    y = -20.87
    ax.errorbar(
            x, y, xerr=0.5, yerr=0.03, c='red', marker='o')
    ax.text(x, y*1.005, "18cow", fontsize=textsize, color='red')
    

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
    plum = -20.2

    ax.errorbar(
            trise, plum, yerr=0.1, xerr=0.19, marker='o', c='red')
    ax.text(
            trise, plum/1.005, "16asu", fontsize=textsize, color='red',
            horizontalalignment='right', verticalalignment='top')


def rest2018(ax):
    """
    plot the data from Rest (2018)
    """
    trise = 1.4
    plum = -18.8
    ax.scatter(
            trise, plum, marker='o', c='k')
    ax.text(
            trise, plum*1.005, "15K", fontsize=textsize, color='k')


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
    x = 3.81
    y = -20.26
    ax.errorbar(
           x, y, yerr=0.1, c='k', marker='o')
    ax.arrow(
           x, y, -0.3, 0, color='k', head_width=0.1, head_length=0.1)
    ax.text(x, y/1.005, "04D4ec", fontsize=textsize,
            verticalalignment='top', horizontalalignment='center')

    x = 2.9
    y = -20.39
    ax.errorbar(
           x, y, yerr=0.1, c='k', marker='o')
    ax.text(x, y*1.005, "05D2bk", fontsize=textsize)

    x = 4.6
    y = -20.28
    ax.errorbar(
           x, y, yerr=0.1, c='k', marker='o')
    ax.text(
            x, y*1.005, "06D1hc", fontsize=textsize, 
            horizontalalignment='center')


def ibn(ax):
    kwargs = dict(mec='#b2df8a', marker='s', mfc='white', c='#b2df8a')
    akwargs = dict(length_includes_head=True, head_width=0.05, color='#b2df8a')
    dx = -0.5
    ax.scatter(
            0,0,facecolor='white',
            edgecolor='#b2df8a', label="Ibn", marker='s')

    # SN 1999cq
    ax.errorbar(
        3.9, -19.73, yerr=0.10, **kwargs)

    # PTF11rfh
    # Shri says error bar too big
    # ax.errorbar(
    #         7, -20.49, xerr=4.7, yerr=0.99, color=ibn_c, marker=ibn_m, label="Ibn")

    # LSQ12btw
    ax.errorbar(3.8, -19.44, yerr=0.04, **kwargs)
    ax.arrow(3.3, -19.73, dx, 0, **akwargs)

    # PTF12ldy
    ax.errorbar(6.2, -19.20, xerr=2.0, yerr=0.02, **kwargs)

    # iPTF13beo
    ax.errorbar(1.6, -18.57, xerr=0.9, yerr=0.05, **kwargs)

    # LSQ13ccw
    ax.errorbar(4.7, -18.46, xerr=2.1, yerr=0.06, **kwargs)

    # iPTF14aki
    ax.errorbar(7.0, -19.30, xerr=1.9, yerr=0.03, **kwargs)

    # SN 2014bk
    ax.errorbar(9.9, -19.99, yerr=0.65, **kwargs)
    ax.arrow(9.9, -19.99, dx, 0, **akwargs)

    # SN 2015U
    ax.errorbar(8.8, -19.41, xerr=0.9, yerr=0.27, **kwargs)


    # iPTF15akq
    ax.errorbar(8.3, -18.62, xerr=2.7, yerr=0.31, **kwargs)


def ptf09uj(ax):
    """ Plot PTF09uj """
    ax.scatter(3.5, -19, marker='D', c='lightblue')
    ax.text(3.5+losx, -19+losy, "09uj", fontsize=12)

losx = 0.1
losy = -0.06
fig,ax = plt.subplots(1,1,figsize=(8,6))
sn2018gep(ax)
at2018cow(ax)
iptf16asu(ax)
drout(ax)
arcavi(ax)
rest2018(ax)
ibn(ax)
ptf09uj(ax)

# Add SN2006aj
ax.scatter(0.4, -18.5, c='red', marker='o')
ax.text(
        0.4, -18.6, "06aj", color='red', fontsize=14, 
        horizontalalignment='center', verticalalignment='bottom')
ax.scatter(0.4, -18.4, c='red', marker='o')
ax.text(
        0.4, -18.3, "20bvc", color='red', fontsize=14, 
        horizontalalignment='center', verticalalignment='top')

# Add Koala
ax.scatter(1.83, -20.6, c='red', marker='o')
ax.text(
        1.9, -20.6, "Koala", color='red', fontsize=14, 
        horizontalalignment='center', verticalalignment='bottom')

# Add SN2011kl
ax.scatter(4.97, -20.31, c='red', marker='o')
ax.text(
        5.05, -20.31, "SN2011kl", color='red', fontsize=14, 
        horizontalalignment='left', verticalalignment='top')

# Add DES16X1eho
x = 1.9
y = -20.39
ax.errorbar(
       x, y, xerr=0.6, yerr=0.09, c='k', marker='o')
ax.text(x/1.01, y/1.001, "16X1eho", fontsize=textsize,
    horizontalalignment='right', verticalalignment='top')

ax.set_ylabel("Peak Mag (Rest-frame $g$-band)", fontsize=16)
ax.set_xlim(0,10)
ax.set_ylim(-21.5, -17)
#ax.set_xscale('log')
ax.invert_yaxis()
ax.set_xlabel(
    r"$t_{\mathrm{1/2,rise}}$ (rest frame days)", fontsize=16)

ax.xaxis.set_tick_params(labelsize=14)
ax.yaxis.set_tick_params(labelsize=14)

fig.tight_layout()
plt.savefig("rise_mag.eps", format='eps', dpi=300)

#plt.show()
