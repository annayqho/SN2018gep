""" A luminosity/risetime plot 

In the left panel, show the t_{1/2,rise} in a single band
that is close to 3000 AA in the rest frame.
You can either do this in flux space,
or in magnitude space, in which case it's the time
to increase by 0.75 mag.

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



def tp_mag(dt, mag, z):
    """ get t1/2,peak and peak absolute mag from a LC 
    dt: times in days
    M: apparent mag
    """
    # get the peak absolute mag
    ipeak = np.argmin(mag)
    mpeak = mag[ipeak]
    Mpeak = mpeak - Planck15.distmod(z).value

    # get the time of peak mag
    tend = dt[ipeak]

    # now only keep the rise
    xrise = dt[0:ipeak]
    yrise = mag[0:ipeak]

    # interpolate to get the time at which m=mmax+0.75
    order = np.argsort(yrise)
    tstart = np.interp(mpeak+1, yrise[order], xrise[order])

    # find the time it takes to go from half-peak to peak
    trise = (tend-tstart)/(1+z)

    return trise, Mpeak


def tp_lum(dt, lum, z):
    """ get t12 peak and peak lum from a LC """
    ipeak = np.argmax(lum)
    lpeak = lum[ipeak]

    # get the time of peak lum
    tend = dt[ipeak]

    # now only keep the rise
    xrise = dt[0:ipeak]
    yrise = lum[0:ipeak]

    # interpolate to get the time at which lum = lpeak/2
    order = np.argsort(yrise)
    tstart = np.interp(lpeak/2, yrise[order], xrise[order])

    # find t1/2
    trise = (tend-tstart)/(1+z)

    return trise, lpeak


def sn2018gep(axarr):
    """ Rise time and peak mag, peak lbol, for SN2018gep """
    # for the left panel, load the U-band light curve from UVOT
    z = 0.03154
    dt, filt, mag, emag = get_lc() 
    choose = filt == 'U'
    trise, Mpeak = tp_mag(dt[choose], mag[choose], z=0.03154)
    
    axarr[0].scatter(trise, Mpeak, marker='*', s=300,
            facecolors='black', edgecolors='black')

    axarr[0].text(
            trise*1.1, Mpeak, "SN2018gep", fontsize=14, 
            verticalalignment='bottom', 
            horizontalalignment='left')
    
    # next, load in the bolometric light curve
    dt, lum, llum, ulum = load_lc()
    trise, lpeak = tp_lum(dt, lum, z)

    # plot the right panel
    axarr[1].scatter(trise, lpeak, marker='*', s=300,
            facecolors='black', edgecolors='black')

    axarr[1].text(
            trise, lpeak*1.1, "SN2018gep", fontsize=14, 
            verticalalignment='bottom', 
            horizontalalignment='center')


def at2018cow(axarr):
    """ Rise time and peak mag, peak lbol, for AT2018cow 
    
    using g-band filter
    """
    trise, Mpeak = tp_mag(
            np.array([-0.87, 2.15]), np.array([18.4, 12.9]), z=0.014)
    axarr[0].scatter(
            trise, Mpeak, marker='o', c='k')
    axarr[0].text(
            trise, Mpeak, "AT2018cow", fontsize=12, 
            verticalalignment='bottom', 
            horizontalalignment='left')

    trise = 58288.44-58284.13
    lpeak = 3.5E44
    axarr[1].scatter(
            trise, lpeak, marker='o', c='k')
    axarr[1].text(
            trise, lpeak, "AT2018cow", fontsize=12, 
            verticalalignment='center', 
            horizontalalignment='left')
    

def iptf16asu(axarr):
    """ 
    rise time, peak Mg, peak Lbol for iPTF16asu
    Using the observed g-band filter, which is 3909 rest-frame

    and pseudo-bolometric luminosity, which is also certainly
    a lower limit
    """
    dt = np.array([-2.46, -2.43, -1.61, -1.58, -1.56])
    mag = np.array([19.8, 19.69, 19.34, 19.25, 19.28])
    trise, Mpeak = tp_mag(dt, mag, z=0.187)
    axarr[0].scatter(trise, Mpeak, marker='o', c='k')
    axarr[0].text(trise, Mpeak, "iPTF16asu", fontsize=12, 
                verticalalignment='bottom', 
                horizontalalignment='left')

    # get 16asu LC
    dat = np.loadtxt(datadir + "/bol_lc/iPTF16asu_bolometric_luminosity.txt",
        delimiter=" ")  
    dt = dat[:,0]
    lum = dat[:,1]
    elum = dat[:,2]
    trise, lpeak = tp_lum(dt, lum, z=0.187)
    axarr[1].scatter(trise, lpeak, marker='o', c='k')
    axarr[1].text(trise, lpeak, "iPTF16asu", fontsize=12, 
                verticalalignment='bottom', 
                horizontalalignment='left')
    



def ps1(axarr):
    """ 
    Taken from Table 4 of Drout et al. (2014)
    use the g-band ones, since that's the closest to rest-frame 3000 AA
    only use the gold sample
    """
    trise = []
    Mpeak = []
    lpeak = []

    # PS1-10ah
    z = 0.074
    trise.append(1.0 / (1+z))
    Mpeak.append(-17.59)
    lpeak.append(4E42)

    # PS1-10bjp
    z = 0.113
    trise.append(1.0 / (1+z))
    Mpeak.append(-18.34)
    lpeak.append(8E42)

    # PS1-11qr
    z = 0.324
    trise.append(2.8/(1+z))
    Mpeak.append(-19.84)
    lpeak.append(3E43)

    # PS1-12bb
    z = 0.101
    trise.append(1.8/(1+z))
    Mpeak.append(-16.97)
    lpeak.append(2E42)

    # PS1-12bv
    z = 0.405
    trise.append(3.4/(1+z))
    Mpeak.append(-19.68)
    lpeak.append(2.5E43)

    # PS1-12brf
    z = 0.275
    trise.append(0.9/(1+z))
    Mpeak.append(-18.84)
    lpeak.append(1E43)

    axarr[0].scatter(trise, Mpeak, marker='x', c='k', label="PS1 Gold")
    axarr[1].scatter(trise, lpeak, marker='x', c='k', label="PS1 Gold")


def des(axarr):
    """ Dark Energy Survey transients, gold sample
    
    rise times are strongly upper limits
    """
    trise = np.array([2.6, 12.1, 6.4, 8.6, 9, 5.2, 5.4, 9.2, 3.3, 6.5, 4.9, 5.8, 7.7, 9.7, 3.8, 4.3, 3.5, 6.2, 4.5, 5.0])
    Mpeak = np.array([-19.42, -19.66, -19.69, -16.19, -15.76, -18.98, -16.91, -19.84, -17.30, -19.74, -19.76, -18.23, -19.64, -19.46, -18.41, -22.24, -19.84, -16.04, -20.97, -19.62])
    lpeak = 1E43 * np.array([2.16, 3.05, 1.01, 0.11, 0.08, 2.69, 0.22, 6.71, 0.35, 5.10, 3.56, 2.06, 2.85, 2.37, 1.09, 41.41, 3.87, 0.12, 9.24, 3.29])
    axarr[0].scatter(trise, Mpeak, marker='s', facecolor='white', edgecolor='k', label="DES Gold")
    axarr[1].scatter(trise, lpeak, marker='s', facecolor='white', edgecolor='k', label="DES Gold")


def subaru(axarr):
    """ Subaru rapidly rising transients
    """
    dat = np.loadtxt(datadir + "/subaru.txt", delimiter=',')
    Mpeak = dat[:,1] 
    daymag = dat[:,0] # day per mag
    t12 = 0.75*daymag
    axarr[0].scatter(t12, Mpeak, marker='D', edgecolor='k', facecolor='white', label="Subaru/HSC")


def snia(axarr):
    """ stolen from Fig 7 of Drout+ 14 """
    dat = np.loadtxt(
            "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/snia.txt", 
            delimiter=',')
    t12 = dat[:,0]
    lum = dat[:,1]
    xcenter = 10**(np.average([max(np.log10(t12)), min(np.log10(t12))]))
    ycenter = 10**(np.average([max(np.log10(lum)), min(np.log10(lum))]))
    axarr[1].scatter(t12, lum, alpha=0.5, s=50, c='#e55c30', lw=0)
    axarr[1].text(
        xcenter, ycenter, "SN Ia", fontsize=12,
        horizontalalignment='center', verticalalignment='center')


def slsn(axarr):
    """ stolen from Fig 7 of Drout+ 14 """
    dat = np.loadtxt(
            "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/slsn.csv", 
            delimiter=',')
    t12 = dat[:,0]
    lum = dat[:,1]
    xcenter = 10**(np.average([max(np.log10(t12)), min(np.log10(t12))]))
    ycenter = 10**(np.average([max(np.log10(lum)), min(np.log10(lum))]))
    axarr[1].scatter(t12, lum, alpha=0.5, s=50, c='#f6d746', lw=0)
    axarr[1].text(
        xcenter, ycenter, "SLSN-I", fontsize=12,
        horizontalalignment='center', verticalalignment='center')


def Ibc(axarr):
    """ stolen from Fig 7 of Drout+ 14 """
    dat = np.loadtxt(
            "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/Ibc.csv", 
            delimiter=',')
    t12 = dat[:,0]
    lum = dat[:,1]
    xcenter = 10**(np.average([max(np.log10(t12)), min(np.log10(t12))]))
    ycenter = 10**(np.average([max(np.log10(lum)), min(np.log10(lum))]))
    axarr[1].scatter(t12*2, lum, marker='s', alpha=0.5, s=50, c='#84206b', lw=0, label="Ibc")


def IIn(axarr):
    """ stolen from Fig 7 of Drout+ 14 """
    dat = np.loadtxt(
            "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/IIn.csv", 
            delimiter=',')
    t12 = dat[:,0]
    lum = dat[:,1]
    xcenter = 10**(np.average([max(np.log10(t12)), min(np.log10(t12))]))
    ycenter = 10**(np.average([max(np.log10(lum)), min(np.log10(lum))]))
    axarr[1].scatter(t12*2, lum, marker='D', alpha=0.5, s=50, c='grey', lw=0, label="IIn")

def IIPL(axarr):
    """ stolen from Fig 7 of Drout+ 14 """
    dat = np.loadtxt(
            "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/IIPn.csv", 
            delimiter=',')
    t12 = dat[:,0]
    lum = dat[:,1]
    xcenter = 10**(np.average([max(np.log10(t12)), min(np.log10(t12))]))
    ycenter = 10**(np.average([max(np.log10(lum)), min(np.log10(lum))]))
    axarr[1].scatter(t12*2, lum, marker='x', alpha=1, s=50, c='k', lw=1, label="II(P/L)")


fig,axarr = plt.subplots(1,2,figsize=(10,5), sharex=True)
sn2018gep(axarr)
at2018cow(axarr)
iptf16asu(axarr)
ps1(axarr)
des(axarr)
subaru(axarr)

ax = axarr[0]
ax.set_ylabel("Absolute Magnitude", fontsize=16)
ax.set_ylim(-21.5, -17)
ax.invert_yaxis()
ax.set_xlabel(
    r"Rest-frame $t_\mathrm{1/2, rise}$ [days]", fontsize=16)
ax.legend(loc='lower right', fontsize=12)

ax = axarr[1]
ax.set_ylim(1E41, 1E45)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylabel("$L_\mathrm{bol}$", fontsize=16)
ax.legend(loc='lower left', fontsize=12)
ax.set_xlabel(
    r"Rest-frame $t_\mathrm{1/2, rise}$ [days]", fontsize=16)

for ax in axarr:
    ax.set_xlim(0.05, 100)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.yaxis.set_tick_params(labelsize=14)

fig.tight_layout()
plt.savefig("lum_rise.png")
#plt.show()
