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


def sn2003L(ax, col, legend):
    """ Soderberg et al
    Values at 8.5 GHz """
    d = 2.8432575937224894e+26
    nu_plt = 8.5E9
    nu, dt, f, ef = read_2003L()
    choose = nu == nu_plt
    lum = plot_line(
            ax[1], d, dt[choose], 1E-3*f[choose]*nu_plt, 
            'SN2003L', 'SN', col, legend)
    ax[1].text(dt[choose][-1]/1.05, lum[-1], 'SN2003L', fontsize=11,
            verticalalignment='center',
            horizontalalignment='left')
    

def sn1979c(ax, col, legend):
    """ Weiler 1986 and Weiler 1991
    This is a IIL 
    
    Too old to have the raw table on arXiv
    """
    d = 5.341805643483106e+25
    nu = 1.4E9 # 20cm
    t = np.array(
            [437,594,631,663,679,684,727,747,786,822,839,876,882,
             914,937,973,995,
             1026,1071,1091,1127,1156,1168,1212,1243,1277,1314,1358,1390,
             1415,1435,1466,1513,1565,1600,1634,1659,1698,1714,1750,1771,
             1931,2027])
    flux = np.array(
            [0.2,2.1,2.5,2.7,2.8,2.8,4.4,4.8,6.0,7.1,7.1,7.6,8.6,
             9.8,6.5,8.6,9.5,
             10.2,10.8,10.3,10.4,12.2,10.1,10.2,11.5,11.2,13.0,11.3,10.2,
             9.6,11.2,13.2,11.1,9.1,8.5,9.1,8.8,10.1,9.7,9.1,8.9,
             7.0,7.7])
    lum = plot_line(ax[1], d, t, nu*flux, 'SN1979c', 'SN', col, legend)
    ax[1].text(t[0]/1.05, lum[0], 'SN1979C', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')
    

def sn1993J(ax, col, legend):
    """ SN 1993J from Weiler et al. 
    This is the peak of the 99.4 GHz light curve
    There is also a 110 GHz point, but only one,
    so I can't get the peak.
    1.4 GHz / 20 cm, 4.9 GHz / 6 cm
    values come from the best-fit model,
    but by eye they are clearly pretty close
    """
    d = 1.1e25
    freq = 5E9
    nu, dt, f, ef, islim = read_1993J_low_freq()
    choose = np.logical_and(~islim, nu==freq)
    lum = plot_line(
            ax[1], d, dt[choose], freq*f[choose], 
            'SN1993J', 'SN', col, legend)
    ax[1].text(dt[choose][0]/1.05, lum[0], 'SN1993J', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')

    freq = 99.4E9
    nu, dt, f, ef, islim = read_1993J_high_freq()
    choose = np.logical_and(~islim, nu==freq)
    lum = plot_line(
            ax[0], d, dt[choose], freq*f[choose], 
            'SN1993J', 'SN', col, legend)
    ax[0].text(dt[choose][-1]*1.1, lum[-1], 'SN1993J', fontsize=11,
            verticalalignment='center',
            horizontalalignment='left')


def sn2011dh(ax, col, legend):
    """ SN 2011dh
    Horesh et al. 2013
    Krauss et al. 2012
    M51: d = 8.03 Mpc; expl date May 31.58
    """
    d = 2.5E25

    # HIGH FREQUENCY
    # use two freq: 107E9 and 93E9
    dt, nu, f, ef, islim = read_2011dh()
    choose = np.logical_and(~islim, np.logical_or(nu==107E9, nu==93E9))
    lum = plot_line(
            ax[0], d, dt[choose], nu[choose]*f[choose], 
            'SN2011dh', 'SN', col, legend)
    ax[0].text(dt[choose][0]/1.05, lum[0], 'SN2011dh', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')

    # LOW FREQUENCY
    # from Horesh 2013 and Krauss 2012
    dt_all = []
    f_all = []
    nu_all = []

    freq = 8.5E9
    choose = np.logical_and(~islim, nu==freq)
    dt_all.extend(dt[choose])
    f_all.extend(f[choose])
    nu_all.extend([freq]*sum(choose))

    freq = 6.7E9
    dt_all.extend([16.4, 20.4, 25.4, 35.3, 45.3, 58.2, 92.9])
    f_all.extend([4.09, 4.8, 5.98, 7.222, 6.987, 6.11, 3.941])
    nu_all.extend([freq]*7)

    dt_all = np.array(dt_all)
    f_all = np.array(f_all)
    nu_all = np.array(nu_all)

    lum = plot_line(
            ax[1], d, dt_all, nu_all*f_all, 
            'SN2011dh', 'SN', col, legend)
    ax[1].text(dt[choose][0]/1.05, lum[0], 'SN2011dh', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')


def grb030329(ax, col, legend):
    """ Sheth et al. 2003 
    Berger 2003
    Van der Horst et al. 2008
    
    Explosion day was obviously 03/29
    """
    z = 0.1686
    d = Planck15.luminosity_distance(z=z).cgs.value

    # LOW FREQUENCY

    # Berger: this is the best frequency to pick from this paper
    t = np.array(
            [0.58, 1.05, 2.65, 3.57, 4.76, 6.89, 7.68, 9.49, 11.90, 
                12.69, 14.87, 16.66, 18.72, 20.58, 25.70, 28.44, 31.51, 
                33.58, 36.52, 42.55, 44.55, 59.55, 66.53]) / (1+z)
    f = np.array(
            [3.50, 1.98, 8.50, 6.11, 9.68, 15.56, 12.55, 13.58, 17.70, 
                17.28, 19.15, 17.77, 15.92, 16.08, 15.34, 12.67, 13.55, 
                13.10, 10.64, 8.04, 8.68, 4.48, 4.92])
    nu = np.array([8.5E9]*len(f))

    # Van der Horst: best frequency is 2.3 GHz
    t = np.append(t, np.array([268.577, 306.753, 365.524, 420.168, 462.078, 
        583.683, 743.892, 984.163]) / (1+z))
    f = np.append(
            f, np.array([1613, 1389, 871, 933, 707, 543, 504, 318]) * 1E-3)
    nu = np.append(nu, np.array([2.3E9]*8))
    lum = plot_line(ax, d, t, nu*f, 'GRB030329', 'GRB', col, legend)
    ax.text(t[-8]*1.05, lum[-8]/1.1, 'GRB030329', fontsize=11,
            verticalalignment='bottom',
            horizontalalignment='left')
    

def grb130427A(ax, col, legend):
    """ Perley et al
    """
    z = 0.340
    d = Planck15.luminosity_distance(z=z).cgs.value

    freq = 5.10E9
    t = np.array([0.677, 2.04, 4.75, 9.71, 17.95, 63.78, 128.34]) / (1+z)
    f = np.array([1290, 1760, 648, 454, 263, 151, 86]) * 1E-3
    lum = plot_line(ax, d, t, freq*f, 'GRB130427A', 'GRB', col, legend)
    ax.text(t[0], lum[0]*1.1, 'GRB130427A', fontsize=11,
            verticalalignment='bottom',
            horizontalalignment='left')


def sn2007bg(ax, col, legend):
    """ Salas et al. 2013
    Peak is resolved for 4.86, 8.46 GHz
    """
    nu = 8.46E9
    d = Planck15.luminosity_distance(z=0.0346).cgs.value
    t = np.array(
            [13.8, 19.2, 26.1, 30.9, 41.3, 55.9, 66.8, 81.8, 98.8, 124, 
                144, 159.8, 189.9, 214.9, 250.9, 286.8, 314.8, 368.8, 
                386.8, 419.9, 566.9, 623.8, 720.8, 775.8, 863.8])
    f = np.array(
            [480, 753, 804, 728, 1257, 1490, 1390, 1325, 1131, 957, 
                621, 316, 379, 404, 783, 1669, 2097, 2200, 
                2852, 3344, 3897, 3891, 3842, 3641, 3408]) * 1E-3
    lum = plot_line(ax[1], d, t, nu*f, 'SN2007bg', 'SN', col, legend)
    ax[1].text(t[0]/1.05, lum[0], 'SN2007bg', fontsize=11,
            verticalalignment='bottom',
            horizontalalignment='right')


def sn2003bg(ax, col, legend):
    """ Soderberg et al.
    Peak is resolved for 22.5, 15, 8.46, 4.86, 1.43
    Again, there are two peaks...
    Let's choose the first peak, 8.46
    """
    nu = 8.46E9
    d = 6.056450393620008e+25

    t = np.array(
                [10, 12, 23, 35, 48, 58, 63, 73, 85, 91, 115, 129,
                132, 142, 157, 161, 181, 201, 214, 227, 242, 255,
                266, 285, 300, 326, 337, 351, 368, 405, 410, 424,
                434, 435, 493, 533, 632, 702, 756, 820, 902, 978])
    f = np.array(
                [2.51, 3.86, 12.19, 24.72, 40.34, 51.72, 49.64, 46.20,
                38.638, 33.85, 45.74, 53.94, 54.27, 54.83, 48.43,
                47.43, 35.76, 31.35, 28.67, 27.38, 24.57, 22.30,
                21.67, 21.31, 20.88, 20.33, 19.85, 18.84, 17.14,
                14.61, 14.49, 14.16, 13.25, 13.08, 10.04, 8.92,
                6.23, 6.18, 4.62, 3.93, 4.69, 4.48])
    lum = plot_line(ax[1], d, t, nu*f, 'SN2003bg', 'SN', col, legend)
    ax[1].text(t[0]/1.05, lum[0], 'SN2003bg', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')


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
    lum = plot_line(ax, d, t, nu*flux, 'SN2009bb', 'SN Ic-BL', col, legend)
    ax.text(t[0]/1.05, lum[0], '2009bb', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')


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
    lum = plot_line(ax, d, t, nu*f, 'SN1998bw', 'SN Ic-BL', col, legend)
    ax.text(t[0]/1.05, lum[0], '1998bw', fontsize=11,
            verticalalignment='bottom',
            horizontalalignment='right')


def sn2010bh(ax, col, legend):
    d = Planck15.luminosity_distance(z=0.0593).cgs.value
    nu = 5.4E9
    t = np.array([18.93, 29.87, 36.97, 69.87])
    f = np.array([90, 128, 68, 50])*1E-3
    lum = plot_line(ax, d, t, nu*f, 'SN2010bh', 'SN Ic-BL', col, legend)
    ax.text(t[-1], lum[-1]/1.3, 'SN2010bh', fontsize=11,
            verticalalignment='top',
            horizontalalignment='center')


def at2018gep(ax, col, legend):
    d = Planck15.luminosity_distance(z=0.033).cgs.value
    t = np.array([4, 5, 15, 16])
    nu = np.array([15, 10, 15, 9]) *1E9
    fnu = np.array([35, 34, 43, 24.4])
    nufnu = nu*fnu
    islim = np.array([True, False, True, False])
    lum = nufnu * 1e-23 * 1e-6 * 4 * np.pi * d**2
    fs = 11
    nsize = 10 # normal size for points
    ax.scatter(
            t[~islim], lum[~islim], facecolor='k', edgecolor='k',
            marker='*', s=100, zorder=10)
    ax.plot(t[~islim], lum[~islim], c='k')
    ax.text(t[-1], lum[-1]/2, 'AT2018gep', fontsize=11,
            horizontalalignment='right',
            verticalalignment='top')


def iptf16asu(ax, col, legend):
    nu = 10E9
    t = np.array([30, 200])
    f = np.array([1E38, 1E38])
    fs = 11
    nsize = 10 # normal size for points
    ax.scatter(
            t, f, facecolor=col, edgecolor='k',
            marker='v', s=70, zorder=10)
    ax.plot(
            t, f, c=col, ls='--')
    ax.text(
            t[-1]*1.4, f[-1], 'iPTF16asu', fontsize=11,
            horizontalalignment='left', verticalalignment='center')


def iptf17cw(ax, col, legend):
    d = Planck15.luminosity_distance(z=0.093).cgs.value
    nu = 6.2E9
    t = np.array([12.6, 15.7, 21.6, 30.7, 41.6])
    f = np.array([38.1, 30.4, 19.9, 22.4, 19])*1E-3
    lum = plot_line(ax, d, t, nu*f, 'iPTF17cw', 'SN Ic-BL', col, legend)
    ax.text(t[0]/1.05, lum[0], 'iPTF17cw', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')


def sn2006aj(ax, col, legend):
    d = 145*3.086E24
    nu = 8.46E9
    t = np.array([1.87, 3, 3.83, 4.85, 6.97, 7.94, 9.95, 12.88,
        16.74, 19.86, 21.96, 24.91, 30.71, 34.81, 41.74, 50.70, 104.52])
    f = np.array([453, 381, 269, 280, 164, 30, 39, 15, 75, 48, 87, 20, 32, 15,
        22, 25, 17])*1E-3
    lum = plot_line(ax, d, t, nu*f, 'SN2006aj', 'SN Ic-BL', col, legend)
    ax.text(t[0]/1.05, lum[0], 'SN2006aj', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')



if __name__=="__main__":
    fig, ax = plt.subplots(1, 1, figsize=(8,6), sharex=True, sharey=True)
    props = dict(boxstyle='round', facecolor='white')

    sncol = '#f98e09' # yellow
    grbcol = '#57106e' # purple

    at2018gep(ax, 'black', None)

    grb030329(ax, grbcol, legend=True)
    grb130427A(ax, grbcol, None)

    iptf16asu(ax, sncol, None)
    sn2009bb(ax, sncol, legend=True)
    sn1998bw(ax, sncol, None)
    iptf17cw(ax, sncol, None)
    sn2010bh(ax, sncol, None)
    sn2006aj(ax, sncol, None)

    ax.set_ylabel(
            r"Luminosity $\nu L_{\nu}$ [erg\,s$^{-1}$]", 
            fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    ax.set_xlim(0.3, 2000) 
    ax.set_ylim(1E34, 1E42)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r"Time [days; rest frame]", fontsize=16)
    ax.legend(fontsize=12, loc='upper right')

    plt.subplots_adjust(wspace=0.05)
    #plt.show()
    plt.savefig("radio_lc.png")
