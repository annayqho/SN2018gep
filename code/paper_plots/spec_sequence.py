""" Plot the spectral sequence """

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15
from astropy.time import Time
import glob
from plot_lc import get_lc
import sys
sys.path.append("/Users/annaho/Github/Spectra")
from normalize import smooth_spec
from measure_snr import get_snr


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
        elif tel == 'omr.ascii':
            # first Xinglong spectrum
            t = Time('2018-09-21T11:15:10.0').jd
            dt[ii] = t-t0
        elif tel == 'Bfosc.ascii':
            # second Xinglong spectrum
            t = Time('2018-09-25T11:16:43.0').jd
            dt[ii] = t-t0
        else:
            print("couldn't find telescope")
            print(tel)
    order = np.argsort(dt)
    files_sorted = files[order]
    dt_sorted = dt[order]
    tel_sorted = np.array(tels)[order]
    cols = cols[order]
    return files_sorted, dt_sorted, tel_sorted


def get_res(tel):
    """ Here, this means the width of a line in Angstroms """
    if tel == 'LT':
        res = 18 # Angstrom, res at central wavelength
        res = 30 # add a couple of Ang?
    elif tel == 'P200':
        res = 10 # determined by eye from the spectrum
        # basically, width of a galaxy emission line is 10 AA
        # and each pixel is 1 AA
    elif tel == 'Keck1':
        res = 7*2 # determined by eye from spectrum
        # width of a line is around 7 pixels
        # and each pixel is 2 Angstroms
    elif tel == 'NOT':
        # width of a line is around 8 pixels
        # and each pixel is around 2.63 Ang
        res = 8*2.63
    elif tel == 'DCT':
        # width of a line is around 7 pixels
        # and each pixel is 2.2 Ang
        res = 7*2.2
    elif tel == 'P60':
        res = 20
    elif 'ascii' in tel:
        # Xinglong spectrum
        res = 26
    else:
        res = 1
    return res


def load_spec(f, tel):
    """ load data from spec file """
    dat = np.loadtxt(f)
    wl = dat[:,0]
    flux = dat[:,1]
    if tel == 'Keck':
        eflux = dat[:,3]
    else:
        # need to estimate uncertainty from scatter
        eflux = np.array([get_snr(wl, flux, 6000, 6200)]*len(wl))
    ivar = 1/eflux**2
    return wl, flux, ivar


def plot_spec(ax, x, y, tel, epoch):
    """ plot the spectrum """
    choose_x = np.logical_and(x >= 3200, x<= 9300)
    choose = choose_x
    ax.plot(
            x[choose], y[choose], c='grey', 
            drawstyle='steps-mid', lw=0.5, alpha=0.4)
    return ax


def plot_smoothed_spec(ax, x, y, ivar, tel, epoch):
    """ plot the smoothed spectrum """
    res = get_res(tel)
    choose_x = np.logical_and(x >= 3200, x<= 9300)
    choose = choose_x 
    smoothed = smooth_spec(x, y, ivar, res*3)
    ax.plot(
            x[choose], smoothed[choose], c='black', 
            drawstyle='steps-mid', lw=0.5, alpha=1.0)
    dt_str = r"+%s\,d" %str(np.round(epoch, 1))
    ax.text(
            x[choose][-1]+100, smoothed[choose][-1],  s=dt_str, 
            horizontalalignment='left', verticalalignment='center', 
            fontsize=14)
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


def plot_lines(ax, z, tel, dt):
    """ Plot galaxy emission lines for a particular redshift """
    res = get_res(tel)
    gal_wl = choose_lines(z, dt)
    for val in gal_wl:
        ax.axvspan(
                val-res/2, val+res/2, ls='--', color='grey', lw=0.5, alpha=0.5)


def clip_lines(wl, flux, z, tel, dt):
    res = get_res(tel)
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
    dt, filt, mag, emag = get_lc()
    det = np.logical_and(mag<99, ~np.isnan(mag))
    nondet = np.logical_or(mag==99, np.isnan(mag))
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
    #scale = (flam/flam_meas)/1E-15
    return wl, flux


if __name__=="__main__":
    z = 0.03154

    files, epochs, tels = get_files()
    start = 0
    end = 19
    files = files[start:end]
    epochs = epochs[start:end]
    tels = tels[start:end]
    nfiles = len(files)

    fig,axarr = plt.subplots(
            1, 2, figsize=(10,10), sharex=True)

    for ii,f in enumerate(files):
        if ii < nfiles/2:
            ax = axarr[0]
        else:
            ax = axarr[1]
        tel = tels[ii]
        dt = epochs[ii]
        wl, flux, ivar = load_spec(f, tel)
        print(tel)
        wl, flux = fluxcal(wl, flux, dt)
        wl, flux = clip_lines(wl, flux, z, tel, dt)
        #plot_lines(ax, z, tel, dt)
        wl, flux = clip_tellurics(wl, flux)
        wl, flux = fluxcal(wl, flux, dt)
        scale = (flux[wl > 4100][0])/2
        plot_spec(ax, wl, flux/scale+nfiles/2-ii%(nfiles/2), tel, dt)
        plot_smoothed_spec(
                ax, wl, flux/scale+nfiles/2-ii%(nfiles/2), ivar, tel, dt)
        ax.tick_params(axis='both', labelsize=14)
    axarr[0].set_ylabel(
            r"Scaled $F_{\lambda}$ + constant",
            fontsize=16)
    axarr[0].set_xlabel(r"Observed Wavelength (\AA)", fontsize=16)
    axarr[1].set_xlabel(r"Observed Wavelength (\AA)", fontsize=16)
    axarr[1].get_yaxis().set_ticks([])
    plt.xlim(3000, 11000)
    #plt.xlim(4900, 5200)
    plt.subplots_adjust(wspace=0)
    axarr[0].set_ylim(0,11)
    #axarr[0].set_ylim(0,5)

    #plt.tight_layout()
    #plt.savefig("spec_sequence.png")
    plt.show()
    #plt.close()
