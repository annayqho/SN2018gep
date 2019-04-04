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


def sn2018gep(axarr):
    """ Rise time and peak mag, peak lbol, for SN2018gep """
    z = 0.03154

    # for the left panel, load the U-band light curve from UVOT
    dt, filt, mag, emag = get_lc() 
    choose = filt == 'U'

    # get the peak absolute mag
    ipeak = np.argmin(mag[choose])
    mpeak = mag[choose][ipeak]
    Mpeak = mpeak - Planck15.distmod(z=z).value

    # get the time of peak mag
    tend = dt[choose][ipeak]

    # now only keep the rise
    xrise = dt[choose][0:ipeak]
    yrise = mag[choose][0:ipeak]

    # interpolate to get the time at which m=mmax+0.75
    order = np.argsort(yrise)
    tstart = np.interp(mpeak+1, yrise[order], xrise[order])

    # find the time it takes to go from half-peak to peak
    trise = (tend-tstart)/(1+z)

    axarr[0].scatter(trise, Mpeak, marker='*', s=300,
            facecolors='black', edgecolors='black')

    axarr[0].text(
            trise*1.1, Mpeak, "SN2018gep", fontsize=14, 
            verticalalignment='bottom', 
            horizontalalignment='left')

    # next, load in the bolometric light curve
    dt, lum, llum, ulum = load_lc()

    # get the peak lum
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

    # plot the right panel
    axarr[1].scatter(trise, lpeak, marker='*', s=300,
            facecolors='black', edgecolors='black')

    axarr[1].text(
            trise*1.1, lpeak, "SN2018gep", fontsize=14, 
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


def ps1(axarr):
    """ equivalent rise times, bol lum for the PS1 sample

    extracted from Fig 7 of Drout+14
    they say that if they used a pseubolometric correction,
    the luminosities would be roughly a factor of 2 higher
    I think the more accurate thing is to show the real points,
    and then put an arrow
    """
    t12 = np.array([6.706972119044536, 10.84536995642165, 12.04667856053548,
        7.718257367165403, 5.769891242129365, 7.939867719009442, 
        10.923585374719543, 5.716243487271193, 8.477348029364208, 
        8.711492903371369])
    lum = np.array([2.4992169158787695e+43, 2.4586509528158067e+43, 
        1.261169269861789e+43, 1.0022730413557269e+43, 4.622323136168413e+42,
        4.2361218340535585e+42, 1.5844083921328312e+42, 3.0314287872197434e+43,
        2.77862997924052e+43, 7.567893474912679e+42])
    # show a box around these
    xwidth = max(t12)-min(t12)
    xcenter = 10**(np.average([max(np.log10(t12)), min(np.log10(t12))]))
    ywidth = max(lum)-min(lum)
    ycenter = 10**(np.average([max(np.log10(lum)), min(np.log10(lum))]))
    rect = Rectangle(
            xy=(min(t12), min(lum)),
            width=xwidth, height=ywidth)
    pc = PatchCollection([rect], facecolor='grey', alpha=0.5, edgecolor='k')
    axarr[1].add_collection(pc)
    axarr[1].arrow(
            xcenter, max(lum), 0, 3E43, length_includes_head=True,
            head_width=2, head_length=1E44/7, fc='k')
    axarr[1].text(
        xcenter, ycenter, "PS1", fontsize=12, 
        horizontalalignment='center', verticalalignment='center')


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
#at2018cow(axarr)
#iptf16asu(axarr)
#ps1(axarr)
#snia(axarr)
#slsn(axarr)
#Ibc(axarr)
#IIn(axarr)
#IIPL(axarr)

ax = axarr[0]
ax.set_ylabel("Abs. Mag. (Rest-frame 3000\,\AA)", fontsize=16)
ax.invert_yaxis()
ax.set_xlabel(
    r"Rest-frame $t_\mathrm{1/2, rise}$ [days]", fontsize=16)

ax = axarr[1]
ax.set_ylim(1E41, 1E45)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylabel("$L_\mathrm{bol}$", fontsize=16)
ax.legend(loc='lower left', fontsize=12)
ax.set_xlabel(
    r"Rest-frame $t_\mathrm{1/2, rise}$ [days]", fontsize=16)

for ax in axarr:
    ax.set_xlim(0.05, 1000)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.yaxis.set_tick_params(labelsize=14)

fig.tight_layout()
#plt.savefig("lum_rise.png")
plt.show()
