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


# Where the observed light curve compilation lives
ddir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/lc"

def at2018gep(ax, shift=3):
    """ LC of AT2018gep
    
    shifted to the left by some number of days. 
    for example, to align with g-band max,
    you would shift the whole thing left by 3 days
    """
    DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/phot"
    f = DATA_DIR + "/ZTF18abukavn_opt_phot.dat"
    dat = np.loadtxt(f, dtype=str, delimiter=' ')
    instr = dat[:,0]
    jd = dat[:,1].astype(float)
    zp = 2458370.6473
    dt = jd-zp-shift
    filt = dat[:,2]
    mag = dat[:,3].astype(float)
    emag = dat[:,4].astype(float)
    distmod = Planck15.distmod(z=0.033).value
    choose = np.logical_and(mag<99, filt=='u')
    order = np.argsort(dt[choose])

    ls=':'
    lw=1.0
    
    ax.plot(
            dt[choose][order], mag[choose][order]-distmod, 
            ls=ls, c='cyan', lw=lw)
    choose = np.logical_and(mag<99, filt=='g')
    order = np.argsort(dt[choose])
    ax.plot(
            dt[choose][order], mag[choose][order]-distmod, 
            ls=ls, c='g', lw=lw)
    choose = np.logical_and(mag<99, filt=='r')
    order = np.argsort(dt[choose])
    ax.plot(
            dt[choose][order], mag[choose][order]-distmod, 
            ls=ls, c='r', lw=lw)
    choose = np.logical_and(mag<99, filt=='i')
    order = np.argsort(dt[choose])
    ax.plot(
            dt[choose][order], mag[choose][order]-distmod, 
            ls=ls, c='pink', lw=lw)
    choose = np.logical_and(mag<99, filt=='z')
    order = np.argsort(dt[choose])
    ax.plot(
            dt[choose][order], mag[choose][order]-distmod, 
            ls=ls, c='grey', lw=lw)


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


def drout(axarr):
    """ Comparison to the transients from Drout+14 """
    # redshift key
    z = {'10ah':0.074, '10bjp':0.113, '11qr':0.324, '12bb':0.101, '12bv':0.405,
         '12brf':0.275, '11bbq': 0.646, 
         '13duy':0.270, '13dwm':0.245, '13ess':0.296}
    inputf = ddir + "/drout14.txt"
    dat = np.loadtxt(inputf,delimiter=';',dtype=str)
    names = np.array([val.strip() for val in dat[:,0]])
    unames = np.unique(names)
    for ii,ax in enumerate(axarr.reshape(-1)):
        name = unames[ii]
        choose = names == name
        redshift = z[name]
        # only do this if there is a known redshift
        filt = np.array(dat[:,1][choose]).astype(str)
        dt = dat[:,2][choose].astype(float)
        islim = dat[:,3][choose]
        mag = dat[:,4][choose].astype(float)
        distmod = Planck15.distmod(z=redshift).value
        isdet = np.array([val == " " for val in islim])
        isfilt = np.array(['g' in val for val in filt])
        c = np.logical_and(isdet, isfilt)
        ax.plot(
                dt[c], mag[c]-distmod, 
                ls='-', c='g')
        isfilt = np.array(['r' in val for val in filt])
        c = np.logical_and(isdet, isfilt)
        ax.plot(
                dt[c], mag[c]-distmod, 
                ls='-', c='r')
        isfilt = np.array(['i' in val for val in filt])
        c = np.logical_and(isdet, isfilt)
        ax.plot(
                dt[c], mag[c]-distmod, 
                ls='-', c='pink')
        isfilt = np.array(['z' in val for val in filt])
        c = np.logical_and(isdet, isfilt)
        ax.plot(
                dt[c], mag[c]-distmod, 
                ls='-', c='grey')
        ax.text(
                0.9, 0.9, 'PS1-%s' %name, fontsize=12, 
                transform=ax.transAxes, 
                horizontalalignment='right', verticalalignment='top')


def gen_lc_comparison():
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


if __name__=="__main__":
    # Drout LC comparison
    fig,axarr = plt.subplots(3,3, figsize=(9,6), sharex=True, sharey=True)
    drout(axarr)
    axarr[0,0].set_xlim(-5,20)
    axarr[0,0].set_ylim(-21,-15)
    axarr[0,0].invert_yaxis()
    axarr[2,1].set_xlabel("Days from g-band max", fontsize=14)
    axarr[1,0].set_ylabel("Absolute Magnitude", fontsize=14)
    for ax in axarr.reshape(-1):
        at2018gep(ax, shift=3)
        ax.tick_params(axis='both', labelsize=14)
    plt.subplots_adjust(hspace=0, wspace=0)
    plt.savefig("drout_obs_lc_comparison.png")
    #plt.show()
