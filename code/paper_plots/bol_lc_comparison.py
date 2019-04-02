""" Comparison of bolometric light curves """

import numpy as np
import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import sys
from astropy.table import Table
from astropy.cosmology import Planck15
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/code")
from load_lum import load_lc
from calc_nickel import get_qdep


# Where the bolometric light curve compilation lives
ddir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/bol_lc"


def at2018gep(ax):
    """ Bolometric LC of AT2018gep """
    dt, lum, llum, ulum = load_lc()
    ax.errorbar(dt, lum, yerr=[llum,ulum], fmt='s', c='k', zorder=10,
            label="SN2018gep")
    # ax.text(dt[4]/1.05, lum[4], 'AT2018gep', fontsize=14,
    #         horizontalalignment='right', verticalalignment='center')


def llgrb(ax):
    """ Comparison to LLGRB-SNe """
    col = 'grey'

    # SN1998bw
    dat = Table.read(ddir + "/sn1998bw.dat", format='ascii.fast_no_header')
    dt = dat['col1']
    lum = dat['col2']
    #ax.scatter(dt, lum, marker='.', c=col)
    ax.plot(dt, lum, c=col, ls='-', lw=2, alpha=0.5, label="Ic-BL")
    # ax.text(dt[-1]*1.01, lum[-1], 'SN1998bw', fontsize=14, 
    #         horizontalalignment='left',
    #         verticalalignment='center')

    # SN2010bh
    dat = Table.read(ddir + "/sn2010bh.dat", format='ascii.fast_no_header')
    dt = dat['col1']
    lum = dat['col2']
    #ax.scatter(dt, lum, marker='.', c=col)
    ax.plot(dt, lum, c=col, ls='-', lw=2, alpha=0.5, label='_nolegend_')
    # ax.text(dt[8], lum[8]/1.15, 'SN2010bh', fontsize=14, 
    #         horizontalalignment='center',
    #         verticalalignment='top')

    # SN 2006aj
    dat = Table.read(ddir + "/sn2006aj.dat", format='ascii.fast_no_header')
    dt = dat['col1']
    order = np.argsort(dt)
    dt = dt[order]
    lum = dat['col2'][order]
    ax.plot(dt, lum, c=col, ls='-', lw=2, alpha=0.5, label='_nolegend_')
    # ax.text(dt[6], lum[6]*1.02, 'SN2006aj', fontsize=14, 
    #         horizontalalignment='center',
    #         verticalalignment='bottom')

    # ax.text(0.9, 0.9, 'LLGRB-SNe', fontsize=14,
    #         horizontalalignment='right',
    #         verticalalignment='top', transform=ax.transAxes)


def fbot(ax):
    """ AT2018cow """
    lsun = 3.839E33
    dat = Table.read(
            ddir + "/at2018cow.dat", delimiter='&', format='ascii.fast_no_header')
    mjd = dat['col1']
    jd0 = 58285.441 # time of optical discovery
    dt = mjd-jd0
    lum_raw = dat['col2']
    lum = lsun * np.array(
        [val.split('^')[0] for val in lum_raw]).astype(float)
    ulum = lsun*(np.array(
        [val.split('^')[1].split('_')[0] for val in lum_raw]).astype(float))
    llum = lsun*(np.array(
        [val.split('^')[1].split('_')[1] for val in lum_raw]).astype(float))
    ax.plot(dt, lum, c='grey', ls='-', lw=2, alpha=0.5, label="AT2018cow")
    # ax.text(dt[4]*1.05, lum[4], 'AT2018cow', fontsize=14,
    #         horizontalalignment='left',
    #         verticalalignment='center')

    # ax.text(0.9, 0.9, 'FBOTs', fontsize=14,
    #         horizontalalignment='right',
    #         verticalalignment='top', transform=ax.transAxes)


def sn2008d():
    dat = Table.read(
            ddir + "/sn2008d.dat", delimiter='&', format='ascii.fast_no_header')
    dt = dat['col1']
    TBB_raw = dat['col2']
    TBB = np.array([val.split('+')[0] for val in TBB_raw])
    uTBB = np.array([val.split('+')[1].split('_')[0] for val in TBB_raw])
    lTBB = np.array([val.split('+')[1].split('_')[1] for val in TBB_raw])
    lum_raw = dat['col4']
    lum = 10**(np.array(
        [val.split('+')[0] for val in lum_raw]).astype(float))
    ulum = 10**(np.array(
        [val.split('+')[1].split('_')[0] for val in lum_raw]).astype(float))
    llum = 10**(np.array(
        [val.split('+')[1].split('_')[1] for val in lum_raw]).astype(float))
    ax.errorbar(dt, lum, yerr=[llum,ulum], fmt='.', c='grey')
    ax.plot(dt, lum, c='grey')
    #ax.fill_between(dt, y1=lum-ulum, y2=lum+llum, color='grey')
    ax.text(
            dt[8], lum[8], 'SN2008D', fontsize=14, 
            horizontalalignment='left',
            verticalalignment='top')


if __name__=="__main__":
    # Initialize figure
    fig,axarr = plt.subplots(1,2, figsize=(9,5), sharex=True, sharey=True)
    ax = axarr[0]
    llgrb(ax)
    at2018gep(ax)
    ax.legend(fontsize=12)

    # Show the nickel decay line
    x = np.linspace(15,40)
    t0 = 30 
    mni = 0.28
    y = get_qdep(x, t0, mni)
    ax.plot(x, y, ls='--', lw=0.5, c='k')
    ax.text(26, 3E42, "$M_\mathrm{Ni} = 0.28\,M_\odot$", fontsize=11)

    ax = axarr[1]
    at2018gep(ax)
    fbot(ax)
    ax.legend(fontsize=12)

    # Formatting
    for ax in axarr:
        ax.tick_params(axis='both', labelsize=14)
        ax.set_yscale('log')
        ax.set_xlim(-5, 40)
        ax.set_ylim(5E41, 1E45)
        ax.set_xlabel(r'Days since first light', fontsize=16)
    axarr[0].set_ylabel(r'$L_\mathrm{bol}$ (erg/s)', fontsize=16)

    plt.subplots_adjust(wspace=0)

    #plt.show()
    plt.savefig("bol_lc_comparison.png")
