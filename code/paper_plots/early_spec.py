""" Plot the spectra from the first 5 days """

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15
from astropy.time import Time
from astropy.modeling.blackbody import blackbody_lambda
import glob
from plot_lc import get_lc
import sys
sys.path.append("/Users/annaho/Github/Spectra")
sys.path.append("/Users/annaho/Github/sigfig")
from round_nums import round_sig
from normalize import smooth_spec
from measure_snr import get_snr
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/code")
from load_lum import load_lc
from load_radius import load_radius
from load_temp import load_temp


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


def get_temp(epoch):
    dt_all, temp_all, ltemp, utemp = load_temp()
    temp = np.interp(epoch, dt_all, temp_all)
    return temp


def load_spec(f, tel):
    """ load data from spec file """
    dat = np.loadtxt(f)
    wl = dat[:,0]
    flux = dat[:,1]
    # the DCT spec is shifted by like 15 AA
    if tel == 'DCT':
        wl = wl - 15
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


def plot_smoothed_spec(ax, x, y, ivar, tel, epoch, ls='-', lw=0.5, c='black', label=None):
    """ plot the smoothed spectrum """
    res = get_res(tel)
    choose_x = np.logical_and(x >= 3200, x<= 9300)
    choose = choose_x 
    smoothed = smooth_spec(x, y, ivar, res*3)
    ax.plot(
            x[choose], smoothed[choose], c=c, 
            drawstyle='steps-mid', lw=lw, ls=ls, alpha=1.0, label=label)
    dt_str = r"+%s\,d ($T=%s\,$kK)" %(
            str(np.round(epoch, 1)), (int(round_sig(temp/1000))))
    ax.text(
            x[choose][-1]+100, smoothed[choose][-1],  s=dt_str, 
            horizontalalignment='left', verticalalignment='center', 
            fontsize=12)
    return smoothed


def choose_gal_lines(z, dt):
    """ choose galaxy emission lines given the epoch """
    balmer = np.array([6564.61, 4862.68, 4341.68, 4102.89, 3970.072])
    oiii = np.array([4363, 4932.6, 4960.295, 5008.24]) # O III
    oii = np.array([3727.092, 3729.875])
    nii = np.array([6549.86])
    oi = np.array([6302.046, 6365.536])
    gal_wl = np.hstack((balmer, oiii, oii)) * (z+1)
    return gal_wl


def get_ciii(v):
    return np.array([4649, 5696])*(1-v/3E5)


def plot_ciii(ax, y, v, label=None):
    """ plot ciii ionization lines, shifted by some velocity 
    
    Parameters
    ----------
    v: velocity given in km/s
    """
    wl = get_ciii(v)
    ax.scatter(
        wl, y, marker='v', c='k', label=label)


def get_cii(v):
    return 6580*(1-v/3E5)


def plot_cii(ax, y, v, label=None):
    """ return some ionization lines, shifted by some velocity and redshift
    
    Parameters
    ----------
    v: velocity given in km/s
    """
    wl = get_cii(v)
    ax.scatter(
            wl, y, marker='v', 
            facecolor='white', edgecolor='k', label=label)


def get_siv(v):
    return 4110*(1-v/3E5)


def plot_siv(ax, y, v, label=None):
    """ return some ionization lines, shifted by some velocity and redshift
    
    Parameters
    ----------
    v: velocity given in km/s
    """
    line = get_siv(v)
    ax.scatter(line, y, marker='d', 
            facecolor='black', edgecolor='k', label=label)


def get_oii(v):
    return np.array([4670,4350])*(1-v/3E5)


def plot_oii(ax, y, v, label=None):
    """ return some ionization lines, shifted by some velocity and redshift
    
    Parameters
    ----------
    v: velocity given in km/s
    """
    wl = get_oii(v)
    ax.scatter(wl, y, marker='d', 
            facecolor='white', edgecolor='k', label=label)
    return wl


def plot_gal_lines(ax, z, tel, dt):
    """ Plot galaxy emission lines for a particular redshift """
    res = get_res(tel)
    gal_wl = choose_gal_lines(z, dt)
    for val in gal_wl:
        ax.axvspan(
                val-res/2, val+res/2, ls='--', color='grey', lw=0.5, alpha=0.5)


def clip_lines(wl, flux, z, tel, dt):
    res = get_res(tel)
    gal_wl = choose_gal_lines(z, dt)
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
    scale = flam/flam_meas
    return wl, flux*scale


if __name__=="__main__":
    z = 0.03154

    files, epochs, tels = get_files()
    start = 0
    end = 9
    files = files[start:end]
    epochs = epochs[start:end]
    tels = tels[start:end]
    nfiles = len(files)
    shift = [0, 0.2, 0.4, 0.7, 1.0, 1.1, 1.6, 2.0, 2.2]

    fig,ax = plt.subplots(
            1, 1, figsize=(8,6), sharex=True)

    done = 0
    for ii,f in enumerate(files):
        tel = tels[ii]
        dt = epochs[ii]

        wl, flux, ivar = load_spec(f, tel)
        
        # Only show the spectra up to 7000AA
        choose = wl < 7000
        wl = wl[choose]
        flux = flux[choose]
        ivar = ivar[choose]

        wl, flux = fluxcal(wl, flux, dt)
        wl, flux = clip_lines(wl, flux, z, tel, dt)
        #plot_lines(ax, z, tel, dt)
        wl, flux = clip_tellurics(wl, flux)
        wl, flux = fluxcal(wl, flux, dt)
        temp = get_temp(dt)
        bb = blackbody_lambda(wl, temp).value
        temp = get_temp(dt)
        scale = flux[wl>4100][0]
        shifted = flux/scale-shift[ii]
        plot_spec(ax, wl, shifted, tel, dt)
        if ii == 7:
            smoothed = plot_smoothed_spec(
                    ax, wl, shifted, ivar, tel, dt, lw=1)
        else:
            smoothed = plot_smoothed_spec(
                    ax, wl, shifted, ivar, tel, dt)
        if ii == 7:
            v = 22000
            wl_ciii = get_ciii(v)
            y_ciii = np.array([smoothed[wl<line][-1]+0.1 for line in wl_ciii])
            plot_ciii(ax, y_ciii, v, label="CIII")
            wl_oii = get_oii(v)
            y_oii = np.array([smoothed[wl<line][-1]+0.1 for line in wl_oii])
            plot_oii(ax, y_oii, v, label="OII")
            wl_cii = get_cii(v)
            y_cii = smoothed[wl<wl_cii][-1]+0.05 
            plot_cii(ax, y_cii, v, label="CII")
            wl_siv = get_siv(v)
            y_siv = smoothed[wl<wl_siv][-1]+0.1
            plot_siv(ax, y_siv, v, label="SIV")
        ax.tick_params(axis='both', labelsize=14)
    ax.set_ylabel(
            r"Scaled $F_{\lambda}$ + constant",
            fontsize=16)
    ax.set_xlabel(r"Observed Wavelength (\AA)", fontsize=16)
    ax.set_xlim(3000, 8440)
    ax.set_ylim(-2,1.5)
    ax.legend(loc='upper right', fontsize=12)
    #axarr[0].set_ylim(0,5)

    plt.tight_layout()
    #plt.savefig("early_spectra.png")
    plt.show()
    #plt.close()
