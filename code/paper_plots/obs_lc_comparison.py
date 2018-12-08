""" Comparison of observed light curves """

import numpy as np
import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import sys
from astropy.table import Table
from astropy.cosmology import Planck15
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/code")
from load_lum import load_lc


# Where the bolometric light curve compilation lives
ddir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/lc"

def at2018gep(ax):
    """ LC of AT2018gep """
    DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/phot"
    f = DATA_DIR + "/ZTF18abukavn_opt_phot.dat"
    dat = np.loadtxt(f, dtype=str, delimiter=' ')
    instr = dat[:,0]
    jd = dat[:,1].astype(float)
    zp = 2458370.6473
    dt = jd-zp
    filt = dat[:,2]
    mag = dat[:,3].astype(float)
    emag = dat[:,4].astype(float)
    distmod = Planck15.distmod(z=0.033).value
    choose = np.logical_and(mag<99, filt=='u')
    order = np.argsort(dt[choose])
    ax.plot(
            dt[choose][order], mag[choose][order]-distmod, 
            ls='--', c='cyan')
    choose = np.logical_and(mag<99, filt=='g')
    order = np.argsort(dt[choose])
    ax.plot(
            dt[choose][order], mag[choose][order]-distmod, 
            ls='--', c='g')
    choose = np.logical_and(mag<99, filt=='r')
    order = np.argsort(dt[choose])
    ax.plot(
            dt[choose][order], mag[choose][order]-distmod, 
            ls='--', c='r')
    choose = np.logical_and(mag<99, filt=='i')
    order = np.argsort(dt[choose])
    ax.plot(
            dt[choose][order], mag[choose][order]-distmod, 
            ls='--', c='pink')
    choose = np.logical_and(mag<99, filt=='z')
    order = np.argsort(dt[choose])
    ax.plot(
            dt[choose][order], mag[choose][order]-distmod, 
            ls='--', c='grey')



def ksn2015k(ax):
    """ Comparison to KSN2015K """
    col = 'black'
    dat = Table.read(ddir + "/ksn2015k.txt", format='ascii.fast_no_header')
    dt = dat['col1'] + 3
    mag = dat['col2']
    ax.plot(dt, mag, c=col, ls='-', lw=3, alpha=1)

    ax.text(0.9, 0.9, 'KSN2015K', fontsize=14,
            horizontalalignment='right',
            verticalalignment='top', transform=ax.transAxes)


def iptf16asu(ax):
    """ Comparison to iPTF16asu """
    dat = np.loadtxt(ddir + "/lc_16asu.txt", delimiter='&', dtype='str')
    t = dat[:,0].astype(float)
    dt = t-t[0]
    band = np.array([val.strip() for val in dat[:,2]])
    mag_raw = dat[:,3]
    mag = np.zeros(len(mag_raw))
    emag = np.zeros(len(mag))
    for ii,val in enumerate(mag_raw):
        if '>' not in val:
            mag[ii] = float(val.split('$pm$')[0])
            emag[ii] = float(val.split('$pm$')[1])

    choose = np.logical_and(mag > 0, band == 'g')
    Mag = mag[choose]-39.862 # distance modulus

    ax.plot(dt[choose]-11, Mag, c='black', label='iPTF16asu', lw=3)
    ax.text(0.9, 0.9, 'iPTF16asu', fontsize=14,
            horizontalalignment='right',
            verticalalignment='top', transform=ax.transAxes)


if __name__=="__main__":
    # Initialize figure
    fig,axarr = plt.subplots(2,2, figsize=(9,6), sharex=True, sharey=True)
    ax = axarr[0,0]
    at2018gep(ax)
    ksn2015k(ax)
    ax.invert_yaxis()
    ax.set_ylabel(r'$M$', fontsize=16)
    ax.tick_params(axis='y', labelsize=14)

    ax = axarr[1,0]
    at2018gep(ax)
    iptf16asu(ax)
    ax.set_ylabel(r'$M$', fontsize=16)

    # Formatting
    for ax in axarr[1,:]:
        ax.tick_params(axis='both', labelsize=14)
        ax.set_xlim(-5, 40)
        ax.set_xlabel(r'Days since first light', fontsize=16)

    plt.subplots_adjust(hspace=0, wspace=0)

    #plt.show()
    plt.savefig("obs_lc_comparison.png")
