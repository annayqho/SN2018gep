""" 
Plot of low-frequency radio luminosity over time
Ideally as close to 8 GHz as possible, but also as well-sampled as possible
"""


import matplotlib
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/AT2018cow/code")
sys.path.append("/Users/annaho/Dropbox/Projects/Research/AT2018cow/data/radio_compilations/Zauderer2011")
from astropy.cosmology import Planck15
from get_radio import *
from scale_fluxes import sma_lc
from read_table import *


def plot_limits(ax, x, y, ratiox, ratioy, col):
    """ Plot two arrows from the point """
    ax.annotate('', xy=(x*ratiox, y), xytext=(x, y),
            arrowprops=dict(
                facecolor=col, headwidth=10, width=1, headlength=7))
    ax.annotate('', xy=(x, y*ratioy), xytext=(x, y),
            arrowprops=dict(
                facecolor=col, headwidth=10, width=1, headlength=7))


def plot_line(ax, d, t, nufnu, name, label, col, legend=False, zorder=1):
    """ Plot a line
    Parameters
    ----------
    nufnu: Hz * mJy 
    name: name of the source
    label: label to use as the legend (which also determines the col)
    """
    lum = nufnu * 1e-23 * 1e-3 * 4 * np.pi * d**2
    fs = 11
    nsize = 10 # normal size for points
    if name=='AT2018gep':
        marker='v'
        fcol = col
        s=70
    else:
        if label=='SN':
            marker='o'
            s=nsize
            fcol = col # fill color
        elif label=='GRB':
            marker='o'
            fcol = 'white' # unfilled
            s=nsize
        elif label=='SN Ic-BL':
            marker='s'
            fcol = col 
            s=nsize
        elif label=='TDE':
            marker='s'
            fcol = 'white' #unfilled
            s=nsize
        if legend:
            ax.plot(t, lum, c=col, ls='-', label=label, zorder=zorder)
        else:
            ax.plot(t, lum, c=col, ls='-', label=None, zorder=zorder)
    ax.scatter(
            t, lum, facecolor=fcol, edgecolor=col, 
            marker=marker, s=s, zorder=zorder)
    return lum


def plot_points(ax, d, nu, t, f, marker, name=None):
    """ Plot set of two points """
    lums = []
    for ii,nuval in enumerate(nu):
        if nuval > 90E9:
            lum = plot_point(ax, d, nuval, t[ii], f[ii], marker, name=name)
            lums.append(lum)
        else:
            lum = plot_point(ax, d, nuval, t[ii], f[ii], marker)
            lums.append(lum)
    ax.plot(
        t, lums, ls='--', c='k')
    return lums


def sn2009bb(ax, col, legend):
    """ expl date Mar 19 """
    nu = 8.46E9
    d = 1.237517263280789e+26
    t_apr = 11 + np.array([5.2, 8.2, 13.2, 15.1, 23.2, 29.1])
    t_may = 11 + 30 + np.array([3.1, 10.1, 13, 20.1, 27])
    t_jun = 11 + 30 + 31 + np.array([6, 17, 26])
    t_jul = 11 + 30 + 31 + 30 + np.array([18.9])
    t_aug = 11 + 30 + 31 + 30 + 31 + np.array([11.8])
    t = np.hstack((t_apr, t_may, t_jun, t_jul, t_aug))
    flux = np.array([24.681, 17.568, 16.349, 13.812, 8.881,
        7.714, 8.482, 6.824, 6.327, 3.294, 4.204, 3.203, 2.392,
        1.903, 1.032, 1.084])
    lum = nu * flux * 1E-3 * 1e-23 * 4 * np.pi * d**2 
    eflux = np.array([0.066, 0.088, 0.107, 0.114, 0.121, 0.095,
        0.098, 0.102, 0.151, 0.118, 0.060, 0.074, 0.082, 0.548, 0.104, 0.091])
    elum = nu * eflux * 1E-3 * 1e-23 * 4 * np.pi * d**2 
    ax.fill_between(
            t, lum-elum, lum+elum, color='orange', alpha=0.5)
    ax.text(t[0], lum[0], 'SN2009bb', fontsize=12, horizontalalignment='right')
    


def sn2012ap(ax, col, legend):
    """ Chakraborti et al. 2014 """
    # reading off their Fig 1
    # expl date is Feb 5 +/- 2 days
    # just giving it 10% uncertainties
    d = 40 * 3.086E24
    nu = np.array([7.5E9, 7.5E9, 7.5E9, 8.5E9])
    t = np.array([15, 18, 24, 34])
    f = np.array([3, 5.5, 4.8, 2.9])
    ef = 0.1*f
    
    lum = nu * f * 1E-3 * 1e-23 * 4 * np.pi * d**2 
    elum = nu * ef * 1E-3 * 1e-23 * 4 * np.pi * d**2
    ax.fill_between(
            t, lum-elum, lum+elum, color='orange', alpha=0.5)
    ax.text(
            t[-1], lum[-1], 'SN2012ap', 
            fontsize=12, horizontalalignment='left',
            verticalalignment='center')



def sn1998bw(ax, col, legend):
    """ SN 1998bw
    
    This is a bit complicated because there are two peaks in the light curve,
    but I am using a frequency that has a main peak rather than a frequecy
    with two clear distinct peaks...
    """
    d = 1.17E26 # cm
    nu = 2.3E9
    t = np.array(
            [11.7, 14.6, 15.7, 16.5, 17.8, 19.7, 21.6, 23.6, 25.9, 26.8, 
             28.8, 30.0, 32.9, 34.7, 36.8, 38.8, 40.0, 45.7, 51.7, 57.7, 
             64.7, 67.7, 80.5])
    f = np.array(
            [19.7, 22.3, 23.5, 23.9, 25.1, 25.3, 20.9, 22.9, 28.0, 28.7, 
            31.1, 31.3, 27.3, 33.5, 31.8, 31, 31.3, 26.8, 23.1, 18.5, 
            15.6, 15.6, 9.6])
    # From Shri's paper: noise is a quadrature sum of 0.5 mJy and 2% of flux
    ef = np.sqrt((0.02*f)**2 + 0.5**2)
    lum = nu * f * 1E-3 * 1e-23 * 4 * np.pi * d**2 
    elum = nu * ef * 1E-3 * 1e-23 * 4 * np.pi * d**2 
    ax.fill_between(
            t, lum-elum, lum+elum, color='grey', 
            label="GRB-selected", alpha=0.5)
    ax.text(t[0], lum[0], 'SN1998bw', fontsize=12, horizontalalignment='right')


def sn2010bh(ax, col, legend):
    """ GRB 100316D """
    d = Planck15.luminosity_distance(z=0.0593).cgs.value
    nu = 5.4E9

    # upper limits
    t = np.array([1.81, 10.57])
    f = np.array([78, 81])*1E-3
    lum = nu * f * 1E-3 * 1e-23 * 4 * np.pi * d**2 
    for ii,tval in enumerate(t):
        ax.arrow(
                tval, lum[ii]*1.5, 0, -0.5*lum[ii], color='grey', 
                length_includes_head=True,
                head_width=1, head_length=5E36) 


    tlims = np.hstack((t, 18.93))
    llims = np.hstack((f, 90E-3)) * nu * 1E-3 * 1e-23 * 4 * np.pi * d**2
    # connect the non-detections to the detections
    ax.plot(tlims, llims, c='grey', ls='--')
    ax.text(tlims[0], llims[0]/1.1, 'SN2010bh', fontsize=12,
            verticalalignment='top')

    t = np.array([18.93, 29.87, 36.97, 69.87])
    f = np.array([90, 128, 68, 50])*1E-3
    ef = np.array([20, 19, 21, 16])*1E-3
    lum = nu * f * 1E-3 * 1e-23 * 4 * np.pi * d**2 
    elum = nu * ef * 1E-3 * 1e-23 * 4 * np.pi * d**2 
    ax.fill_between(
            t, lum-elum, lum+elum, color='grey', alpha=0.5)


def sn2006aj(ax, col, legend):
    """ Soderberg 2006 """
    d = 145*3.086E24
    nu = 8.46E9
    t = np.array([1.87, 3, 3.83, 4.85, 6.97, 9.95, 
        16.74, 19.86, 21.96, 30.71])
    f = np.array([453, 381, 269, 280, 164, 39, 75, 48, 87, 32])*1E-3
    ef = np.array(
            [77, 60, 40, 47, 39, 25, 13, 14, 39, 20]) * 1E-3
    lum = nu * f * 1E-3 * 1e-23 * 4 * np.pi * d**2
    elum = nu * ef * 1E-3 * 1e-23 * 4 * np.pi * d**2
    ax.fill_between(
            t, lum-elum, lum+elum, color='grey',
            alpha=0.5)
    ax.text(1.02*t[-1], lum[-1], 'SN2006aj', fontsize=12, horizontalalignment='left')


def sn2003lw(ax, col, legend):
    """ Soderberg et al. 2004. I'm excluding the point on 4.43 day because
    it just seems off... it's basically a non-detection """
    d = Planck15.luminosity_distance(z=0.105).cgs.value
    t = np.array([1.60, 3.60, 8.46, 11.45, 13.46, 17.43, 19.45,
    22.48, 27.41, 31.41, 35.34, 39.37, 42.43, 52.32, 65.32])
    nu = 8.46E9
    flux = np.array([0.540, 0.249, 0.280, 0.304, 0.448, 0.457, 0.811,
    0.467, 0.675, 0.459, 0.308, 0.647, 0.664, 0.450, 0.533])
    eflux = np.array([0.062, 0.043, 0.049, 0.042, 0.039, 0.041, 0.040,
    0.046, 0.045, 0.047, 0.043, 0.045, 0.061, 0.044, 0.028])
    lum = nu * flux * 1E-3 * 1e-23 * 4 * np.pi * d**2 
    elum = nu * eflux * 1E-3 * 1e-23 * 4 * np.pi * d**2 
    ax.fill_between(
            t, lum-elum, lum+elum, color='grey', alpha=0.5)
    ax.text(t[0], lum[0], 'SN2003lw', fontsize=12, horizontalalignment='left')


#def sn2012bz(ax, col, legend):
#    """ Schulze et al. 2014 """
#    d = Planck15.luminosity_distance(z=0.283).cgs.value


def at2018gep(ax, col, legend):
    d = Planck15.luminosity_distance(z=0.033).cgs.value

    # Two VLA detections
    t = np.array([5, 16])
    nu = np.array([10, 9]) *1E9
    fnu = np.array([34, 24.4])
    efnu = np.array([4, 6.7])
    nufnu = nu*fnu
    elum = nu*efnu
    lum = nufnu * 1e-23 * 1e-6 * 4 * np.pi * d**2
    elum = elum * 1e-23 * 1e-6 * 4 * np.pi * d**2

    fs = 11
    nsize = 10 # normal size for points
    ax.errorbar(
            t, lum, yerr=elum, mfc='k', mec='k',
            marker='s', ms=8, zorder=10, c='k', label="AT2018gep")

    # One VLA non-detection (3-sigma)
    t = 59
    nu = 9 * 1E9
    fnu = 16
    lum = nu * fnu * 1e-23 * 1e-6 * 4 * np.pi * d**2
    ax.arrow(
            t, lum*1.5, 0, -0.5*lum, color='k', length_includes_head=True,
            head_width=1, head_length=5E35) 


def iptf16asu(ax, col, legend):
    nu = 10E9
    t = 30
    lum = 1E38
    ax.arrow(
            t, lum*1.5, 0, -0.5*lum, color='orange', length_includes_head=True,
            head_width=1, head_length=1E37) 
    ax.text(30.6, 1.3E38, "iPTF16asu", horizontalalignment='left',
            verticalalignment='center', fontsize=12)


def iptf17cw(ax, col, legend):
    d = Planck15.luminosity_distance(z=0.093).cgs.value
    nu = np.array([6.2E9, 6.2E9, 6.2E9, 6.2E9, 2.8E9, 14E9])
    t = np.array([12.6, 15.7, 21.6, 30.7, 41.6, 65.6])
    f = np.array([38.1, 30.4, 19.9, 22.4, 44, 20.2]) * 1E-3
    ef = np.array([7.3, 5.9, 6.2, 5.5, 14, 6.8]) * 1E-3
    lum = nu * f * 1e-3 * 1e-23 * 4 * np.pi * d**2 
    elum = nu * ef * 1E-3 * 1e-23 * 4 * np.pi * d**2 
    ax.fill_between(
            t, lum-elum, lum+elum, color='orange', 
            label="Opt.-selected", alpha=0.5)
    ax.text(t[0], lum[0], 'iPTF17cw', fontsize=12, horizontalalignment='right',
            verticalalignment='center')


if __name__=="__main__":
    fig, ax = plt.subplots(1, 1, figsize=(8,6), sharex=True, sharey=True)
    props = dict(boxstyle='round', facecolor='white')

    at2018gep(ax, 'black', None)

    sncol = 'orange'
    iptf16asu(ax, 'orange', None)
    sn2009bb(ax, 'orange', legend=True)
    sn1998bw(ax, sncol, None)
    iptf17cw(ax, sncol, None)
    sn2010bh(ax, sncol, None)
    sn2006aj(ax, sncol, None)
    sn2012ap(ax, sncol, None)
    #sn2003lw(ax, sncol, None)

    ax.set_ylabel(
            r"Luminosity $\nu L_{\nu}$ [erg\,s$^{-1}$]", 
            fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    ax.set_xlim(0, 65) 
    ax.set_ylim(3E36, 5E38)
    #ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r"Time [days; rest frame]", fontsize=16)
    ax.legend(fontsize=12, loc='upper right')

    plt.subplots_adjust(wspace=0.05)
    #plt.show()
    plt.savefig("radio_lc.png")
