""" A plot showing host galaxy properties of FBOTs
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


def sn2018gep(ax):
    """ Rise time and peak mag, peak lbol, for SN2018gep """
    ax.scatter(np.log10(8.11), np.log10(0.12), marker='*', s=300, c='k')


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
                trise[ii], plum[ii], yerr=eplum[ii], fmt='o', c='k')
    

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
                fmt='o', mfc='black', mec='black', c='k')


def bersten2018(ax):
    """
    plot the data from Bersten (2018)
    """
    dat = np.loadtxt(datadir + "/bersten.txt", delimiter=',')
    ax.scatter(dat[:,0], dat[:,1], c='k')


def tanaka2018(ax):
    dat = np.loadtxt(datadir + "/tanaka.txt", delimiter=',')
    ax.scatter(dat[:,0], dat[:,1], c='k')


fig,ax = plt.subplots(1,1,figsize=(5,5.4))
sn2018gep(ax)
bersten2018(ax)
tanaka2018(ax)
# sn2018gep(axarr)
# at2018cow(axarr)
# iptf16asu(axarr)
# rest2018(axarr)
# margutti2018(axarr)

ax.set_ylabel("Absolute magnitude", fontsize=16)
# ax.set_xlim(1,20)
# ax.set_ylim(-21, -16)
ax.set_xscale('log')
ax.invert_yaxis()
ax.set_xlabel(
    r"Rise timescale [day/mag]", fontsize=16)
ax2 = ax.twiny()
ax2.set_xlabel(
    r"Rise rate [mag/day]", fontsize=16)
x_f = lambda x_i: 1/x_i
xmin, xmax = ax.get_xlim()
ax2.set_xlim((x_f(xmin), x_f(xmax)))
ax2.set_xscale('log')
ax2.plot([],[])
ax2.tick_params(axis='both', labelsize=14)
ax.tick_params(axis='both', labelsize=14)




# 
# ax = axarr[1]
# ax.set_ylim(1E41, 1E45)
# ax.set_yscale('log')
# ax.set_xscale('log')
# ax.set_ylabel("Peak Bolometric Luminosity", fontsize=16)
# ax.legend(loc='lower left', fontsize=12)
# ax.set_xlabel(
#     r"Rise Time [rest frame days]", fontsize=16)
# 
# for ax in axarr:
#     ax.xaxis.set_tick_params(labelsize=14)
#     ax.yaxis.set_tick_params(labelsize=14)
# 
fig.tight_layout()
plt.savefig("lum_riserate.png")
#plt.show()
