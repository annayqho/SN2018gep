""" Compare light curves to the Arcavi sample """

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import numpy as np
import sys
from astropy.cosmology import Planck15
from plot_lc import get_lc


datadir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/lc"

# use wavelengths to get extinctions
bands = ['UVW2', 'UVM2', 'UVW1', 'U', 'u', 'B', 'g', 'V', 'r', 'i', 'z']
wl = {}
wl['UVW2'] = 1928
wl['UVM2'] = 2246
wl['UVW1'] = 2600
wl['U'] = 3465
wl['u'] = 3543
wl['B'] = 4392
wl['g'] = 4770
wl['V'] = 5468
wl['r'] = 6231
wl['i'] = 7625
wl['z'] = 9134
col = {}
col['UVW2'] = 'black'
col['UVM2'] = 'black'
col['UVW1'] = 'black'
col['U'] = 'blue'
col['u'] = 'darkblue'
col['B'] = 'lightblue'
col['g'] = 'green'
col['V'] = 'yellow'
col['r'] = 'orange'
col['i'] = 'red'
col['z'] = 'pink'


def plot_18gep(ax, band, c='k'):
    """ Plot 18gep lc for a given band """
    z = 0.03154
    dm = Planck15.distmod(z=z).value
    dt_gep, filt_gep, mag_gep, emag_gep = get_lc()
    isdet_gep = np.logical_and(mag_gep<99, ~np.isnan(mag_gep))
    choose = np.logical_and(filt_gep == band, isdet_gep)
    x = dt_gep[choose]/(1+z)
    y = mag_gep[choose]-dm
    # ax.errorbar(x, y,
    #         yerr=emag_gep[choose],
    #         c='k', fmt='s')
    order = np.argsort(x)
    ax.plot(
            x[order], y[order], c=c, 
            lw=0.5, label="18gep, $%s$" %band, ls='--')


def plot_arcavi(name, band):
    dat = np.loadtxt(datadir + "/arcavi_lc.txt", delimiter=';', dtype=str)
    names = np.array([val.strip() for val in dat[:,0]])
    filt = dat[:,2]
    jd = dat[:,3].astype(float)
    mag = dat[:,5].astype(float)
    emag = dat[:,6] # blanks for upper limits
    islim = emag == ''


def d4ec(ax):
    ax = axarr[0,1]
    name = 'SNLS04D4ec'
    z = 0.593
    offset = 3
    dm = Planck15.distmod(z=z).value
    choose = np.logical_and.reduce((names == name, ~islim, filt=='i'))
    ax.errorbar(
            (jd[choose]-jd[choose][0])/(1+z)+offset, mag[choose]-dm,
            yerr=emag[choose].astype(float), mec='k', mfc='white', fmt='s',
            label="D4ec, $i$", c='k')
    choose = np.logical_and.reduce((names == name, ~islim, filt=='g'))
    ax.errorbar(
            (jd[choose]-jd[choose][0])/(1+z)+offset, mag[choose]-dm,
            yerr=emag[choose].astype(float), mec='k', mfc='white', fmt='D',
            label="D4ec, $i$", c='k')
    plot_18gep(ax, 'g')
    plot_18gep(ax, 'UVW1', c='grey')


def d2bk(ax):
    ax = axarr[1,0]
    name = 'SNLS05D2bk'
    z = 0.699
    dm = Planck15.distmod(z=z).value
    choose = np.logical_and.reduce((names == name, ~islim, filt=='i'))
    # R-band light curve
    ax.errorbar(
            (jd[choose]-jd[choose][0])/(1+z), mag[choose]-dm,
            yerr=emag[choose].astype(float), c='k', fmt='s')
    plot_18gep(ax, 'g')


def d1hc(ax):
    ax = axarr[1,1]
    name = 'SNLS06D1hc'
    z = 0.555
    dm = Planck15.distmod(z=z).value
    choose = np.logical_and.reduce((names == name, ~islim, filt=='i'))
    # R-band light curve
    ax.errorbar(
            (jd[choose]-jd[choose][0])/(1+z), mag[choose]-dm,
            yerr=emag[choose].astype(float), c='k', fmt='s')
    plot_18gep(ax, 'g')


def iptf16asu(ax, shift):
    """ plot the bands of iPTF16asu compared to the closest filters in 18gep

    shift: scoot the transient along the x-axis
    """
    z = 0.187
    dat = np.loadtxt(datadir + "/lc_16asu.txt", delimiter='&', dtype='str') 
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

    # for each available band, plot a solid line in the appropriate color
    uband = np.unique(band)

    for b in uband:
        print("Plotting %s for iPTF16asu %s" %(b,wl[b]))
        choose = np.logical_and(mag > 0, band == b)
        if sum(dt[choose]<50) > 1:
            Mag = mag[choose]-39.862 # distance modulus
            ax.plot((dt[choose])/(1+z)-shift, Mag, c=col[b])

            # find the closest rest-frame band
            wl_rest = wl[b]/(1+z)
            print("In the rest frame, that's %s" %wl_rest)
            wl_all = np.array(list(wl.values()))
            wl_keys = np.array(list(wl.keys()))
            closest_band = wl_keys[np.argmin(np.abs(wl_all-wl_rest))]
            print("The closest band is %s, %s" %(closest_band, wl[closest_band]))
            plot_18gep(ax, closest_band, c=col[b])


if __name__=="__main__":
    fig,axarr = plt.subplots(3,2, figsize=(8,8), sharex=True, sharey=True)

    ax = axarr[0,0]
    iptf16asu(ax)

    ax.set_xlim(-5,50)
    ax.invert_yaxis()

    for ax in axarr[:,0]:
        ax.set_ylabel("Abs Mag", fontsize=14)

    for ax in axarr[2,:]:
        ax.set_xlabel("Rest-frame days", fontsize=14)

    for ax in axarr.reshape(-1):
        ax.yaxis.set_tick_params(labelsize=14)
        ax.xaxis.set_tick_params(labelsize=14)
        #ax.legend(fontsize=14)

    plt.show()
