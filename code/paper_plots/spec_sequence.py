""" Plot the spectral sequence """

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
from scipy.signal import savgol_filter
from astropy.table import Table
from astropy.cosmology import Planck15
from astropy.time import Time
import glob
from plot_lc import get_lc


def get_files():
    files = np.array(glob.glob(
    "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/spec/ZTF18abukavn/*.ascii"))
    dt = np.zeros(len(files))
    tels = []
    cols = np.array([""]*len(dt), dtype='U10')

    # Read in all of the files, pull out the corresponding dates, and sort by date
    t0 = 2458370.6634 # in JD
    for ii,f in enumerate(files):
        tel = f.split("_")[2]
        tels.append(tel)
        alldat = open(f).readlines()
        if tel == 'LT':
            for line in alldat:
                if 'DATE-OBS' in line:
                    obsdate = line[13:36]
                    t = Time(obsdate, format='isot').jd
                    dt[ii] = t-t0
            cols[ii] = 'magenta'
        elif tel == 'P200':
            for line in alldat:
                if 'UT shutter open' in line:
                    obsdate = line[12:35]
                    print(obsdate)
                    t = Time(obsdate, format='isot').jd
                    dt[ii] = t-t0
            cols[ii] = 'lightblue'
        elif tel == 'Keck1':
            for line in alldat:
                if 'DATE_BEG' in line:
                    obsdate = line[13:32]
                    t = Time(obsdate, format='isot').jd
                    dt[ii] = t-t0
            cols[ii] = 'red'
        elif tel == 'DCT':
            obsdate = '2018-09-14T00:00:00' # temporary
            t = Time(obsdate, format='isot').jd
            dt[ii] = t-t0
            cols[ii] = 'yellow'
        elif tel == 'NOT':
            obsdate = '2018-09-17T00:00:00' # temporary
            t = Time(obsdate, format='isot').jd
            dt[ii] = t-t0
            cols[ii] = 'green'
        elif tel == 'P60':
            for line in alldat:
                if 'MJD_OBS' in line:
                    obsdate = float(line[11:25])
                    t = Time(obsdate, format='mjd').jd
                    dt[ii] = t-t0
            cols[ii] = 'black'
        else:
            print("couldn't find telescope")
            print(tel)
    order = np.argsort(dt)
    files_sorted = files[order]
    dt_sorted = dt[order]
    tel_sorted = np.array(tels)[order]
    cols = cols[order]
    return files_sorted, dt_sorted, tel_sorted


def load_spec(f):
    """ load data from spec file """
    dat = np.loadtxt(f)
    wl = dat[:,0]
    flux = dat[:,1]
    return wl, flux


def plot_spec(ax, x, y, tel, epoch):
    """ plot the spectrum """
    choose_y = y >= 0
    choose_x = np.logical_and(x >= 3200, x<= 9300)
    choose = np.logical_and(choose_x, choose_y)
    ax.plot(
            x[choose], y[choose], c='grey', 
            drawstyle='steps-mid', lw=0.5, alpha=0.6)
    return ax


def plot_smoothed_spec(ax, x, y, tel, epoch):
    """ plot the smoothed spectrum """
    choose_y = y >= 0
    choose_x = np.logical_and(x >= 3200, x<= 9300)
    choose = np.logical_and(choose_x, choose_y)

    smoothed = savgol_filter(y[choose], 45, polyorder=2)

    ax.plot(
            x[choose], smoothed, c='black', 
            drawstyle='steps-mid', lw=1.0, alpha=1.0)
    dt_str = r"$\Delta t$=%s d" %str(np.round(epoch, 1))
    ax.text(
            0.98, 0.9, s=dt_str, 
            horizontalalignment='right', verticalalignment='center', 
            fontsize=14, transform=ax.transAxes)
    return ax


def choose_lines(z, dt):
    """ choose galaxy emission lines given the epoch """
    balmer = np.array([6564.61, 4862.68, 4341.68, 4102.89, 3970.072])
    oiii = np.array([4363, 4932.6, 4960.295, 5008.24]) # O III
    oii = np.array([3727.092, 3729.875])
    nii = np.array([6549.86])
    oi = np.array([6302.046, 6365.536])
    gal_wl = np.hstack((balmer, oiii, oii)) * (z+1)
    return gal_wl


def plot_lines(z, tel, dt):
    """ Plot galaxy emission lines for a particular redshift """
    if tel == 'LT':
        res = 18 # Angstrom, res at central wavelength
    else:
        res = 30 # temp 
    #elif tel == 'DBSP':
    gal_wl = choose_lines(z, dt)
    for val in gal_wl:
        plt.axvspan(
                val-res/2, val+res/2, ls='--', color='grey', lw=0.5, alpha=0.5)


def clip_lines(wl, flux, z, tel, dt):
    if tel == 'LT':
        res = 18 # Angstrom, res at central wavelength
        res = 30 # add a couple of Ang?
    else:
        res = 30 # placeholder, I don't actually know what it is
    gal_wl = choose_lines(z, dt)
    for line in gal_wl:
        choose = np.logical_and(wl >= line-res/2, wl <= line+res/2)
        flux = np.interp(wl, wl[~choose], flux[~choose]) # interp over features
    return wl, flux


def get_tellurics():
    start = np.array([7594, 6853])
    end = np.array([7678, 6950])
    return start, end


def clip_tellurics(wl, flux):
    start, end = get_tellurics()
    for ii,beg in enumerate(start):
        choose = np.logical_and(wl >= beg, wl <= end[ii])
        flux = np.interp(wl, wl[~choose], flux[~choose])
    return wl, flux


def plot_tellurics():
    col = 'pink'
    plt.axvspan(7594, 7678, ls='--', color=col, lw=0.5, alpha=0.5)
    plt.axvspan(6853, 6950, ls='--', color=col, lw=0.5, alpha=0.5)


def fluxcal(wl, flux, dt_spec):
    """ Flux-calibrate to R-band light curve """
    # get r-band LC
    dt, filt, det, mag, emag = get_lc()
    choose = np.logical_and(det, filt=='r')
    # interpolate to this epoch
    rval = np.interp(dt_spec, dt[choose], mag[choose])
    # TEMP: r is roughly 658nm +/- 138nm
    # TEMP: assume AB mag
    lam = 6580 # in angstroms
    c = 3E18 # angstrom/s
    fnu = 1E-23 * 3631 * 10**(rval/(-2.5)) # erg/s/cm2/Hz
    flam = fnu * (c/lam**2) # should be erg/s/cm2/AA
    # scale factor
    flam_meas = np.interp(lam, wl, flux)
    scale = (flam/flam_meas)/1E-15
    return wl, flux*scale


if __name__=="__main__":
    z = 0.03154

    files, epochs, tels = get_files()
    files = files[0:1]

    fig,ax = plt.subplots(
            1, 1, figsize=(8,10), sharex=True)

    # Giant y-axis label
    plt.ylabel(
            r"Scaled $F_{\lambda}$ + constant",
            fontsize=16)

    for ii,f in enumerate(files):
        tel = tels[ii]
        dt = epochs[ii]
        wl, flux = load_spec(f)
        wl, flux = clip_lines(wl, flux, z, tel, dt)
        wl, flux = clip_tellurics(wl, flux)
        wl, flux = fluxcal(wl, flux, dt)
        plot_spec(ax, wl, flux, tel, dt)
        plot_smoothed_spec(ax, wl, flux, tel, dt)
    plt.tick_params(axis='both', labelsize=14)
    plt.xlabel(r"Observed Wavelength (\AA)", fontsize=16)
    plt.xlim(3000, 10000)
    #ax.set_ylim(0,4)

    #plt.tight_layout()
    #plt.savefig("spec_first_third.png")
    plt.show()
    #plt.close()
