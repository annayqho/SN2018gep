""" Plot the spectra from the first 5 days """

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from collections import OrderedDict
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15
from astropy.time import Time
from extinction import fitzpatrick99
from astropy.modeling.blackbody import blackbody_lambda
import glob
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/code")
from load_lc import get_lc
sys.path.append("/Users/annaho/Github/Spectra")
sys.path.append("/Users/annaho/Github/sigfig")
from round_nums import round_sig
from normalize import smooth_spec
from measure_snr import get_snr
from load_lum import load_lc
from load_radius import load_radius
from load_temp import load_temp


z = 0.03154
SPEC_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/spec"


def get_files(sind, eind):
    """ 
    start_ind: starting index
    end_ind: end index
    """
    files = np.array(glob.glob(SPEC_DIR + "/ZTF18abukavn/*.ascii"))
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
    return files_sorted[sind:eind], dt_sorted[sind:eind], tel_sorted[sind:eind]


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
            x[choose], y[choose], c='lightgrey', 
            drawstyle='steps-mid', lw=0.4, alpha=0.1)
    return ax


def plot_smoothed_spec(ax, x, y, ivar, tel, epoch, ls='-', lw=0.5, c='black', label=None, text=True):
    """ plot the smoothed spectrum """
    res = get_res(tel)
    temp = get_temp(epoch)
    choose_x = np.logical_and(x >= 3200, x<= 9300)
    choose = choose_x 
    smoothed = smooth_spec(x, y, ivar, res*3)
    ax.plot(
            x[choose], smoothed[choose], c=c, 
            drawstyle='steps-mid', lw=lw, ls=ls, alpha=1.0, label=label,
            zorder=10)
    dt_str = r"+%s\,d ($T=%s\,$kK)" %(
            str(np.round(epoch, 1)), (int(round_sig(temp/1000))))
    if text:
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


def get_lines(species, v):
    """ Return wavelengths of lines shifted to some velocity, 
    taking in the species """
    # These wavelengths are in the rest-frame at that redshift
    if species == "CIII":
        return np.array([4649, 5696])*(1-v/3E5+z)
    elif species == "CIV":
        return np.array([4658,5801])*(1-v/3E5+z)
    elif species == "CII":
        return np.array([3919,4267,6580])*(1-v/3E5+z)
    elif species == "SIV":
        return np.array([4110*(1-v/3E5+z)])
    elif species == "OII":
        return np.array([3727,4670,4350])*(1-v/3E5+z)
    elif species == "OIII":
        return np.array([3760])*(1-v/3E5+z)
    elif species == "OIV":
        return np.array([3730])*(1-v/3E5+z)
    else:
        print("I don't recognize this species!")


def plot_lines(ax, y, v, species, offset):
    """ plot ionization lines, shifted by some velocity 
    
    Parameters
    ----------
    v: velocity given in km/s
    """
    markers = {}
    fc = {}

    # Singly ionized species are white
    markers["OII"] = 'd'
    fc["OII"] = 'white'
    markers["CII"] = 'v'
    fc["CII"] = 'white'

    # Doubly ionized species are colored in black
    markers["CIII"] = 'v'
    fc["CIII"] = 'k'
    markers["OIII"] = 'd'
    fc["OIII"] = 'black'

    # Triply ionized species are not filled in
    markers["SIV"] = '+'
    fc["SIV"] = 'black'
    markers["OIV"] = 'x'
    fc["OIV"] = 'black'
    markers["CIV"] = '1'
    fc["CIV"] = 'k'

    wl = get_lines(species, v)
    # offset: sometimes you have to stack them
    # wl are at the redshift of the transient
    ax.scatter(
        wl/(1+z), y*offset, marker=markers[species], 
        facecolor=fc[species], edgecolor='k', label=species, zorder=20)


def plot_species(ax, v, wl, smoothed, sp, offsets):
    """ Plot lines for a list of species """
    for ii,s in enumerate(sp):
        lines = get_lines(s, v)
        y = np.array([smoothed[wl<line][-1]+0.1 for line in lines])
        plot_lines(ax, y, v, s, offsets[ii])


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


def plot_18cow(ax, scale):
    """ Plot 18cow spectrum scaled by some amount """
    wl_cow, flux_cow, ivar_cow = load_spec(
            SPEC_DIR + "/AT2018cow/AT2018cow_20180621_P200_v3.ascii", 'P200')
    plot_smoothed_spec(
            ax, wl_cow/(1+0.014), flux_cow/scale, ivar_cow,
            'P200', 1.0, lw=1, ls='--', 
            label=r"18cow, +5.353d, $T=26\,$kK")


def plot_10vgv(ax, scale):
    dat = np.loadtxt(SPEC_DIR + "/ptf10vgv.txt", delimiter=',')
    x = (dat[:,0])/(1+0.0142)
    y = dat[:,1]
    ax.plot(
            x, y/scale, lw=0.5, c='#84206b', alpha=1, ls='-', 
            label="10vgv, +2d")


def plot_12gzk(ax, scale):
    dat = np.loadtxt(SPEC_DIR + "/ptf12gzk.txt", delimiter=',')
    x = (dat[:,0])/(1+0.0137)
    y = dat[:,1]
    ax.plot(
            x, y/scale-0.1, lw=0.5, c='#f6c746', ls='-', 
            label="12gzk, +3d", zorder=5)


def plot_98bw(ax, scale):
    dat = np.loadtxt(SPEC_DIR + "/98bw.txt", delimiter=',')
    x = dat[:,0] / (1+0.0085)
    y = dat[:,1]
    ax.plot(x, y/scale, lw=0.5, c='#f6c746', ls='-',
            label="98bw, +18d", zorder=5)


def plot_06aj(ax, scale):
    dat = np.loadtxt(SPEC_DIR + "/06aj.txt", delimiter=',')
    x = dat[:,0] / (1+0.033)
    y = dat[:,1]
    ax.plot(x, y/scale+0.5, lw=0.5, c='#84206b', ls='-',
            label="06aj, +18d", zorder=5)


def spec_evol(ax):
    """ Evolution of the early spectra """
    files, epochs, tels = get_files(0, 9)
    nfiles = len(files)
    shift = [0, 0.2, 0.4, 0.7, 1.0, 1.2, 1.6, 2.0, 2.2]
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
        wl, flux = clip_tellurics(wl, flux)
        wl, flux = fluxcal(wl, flux, dt)
        temp = get_temp(dt)
        bb = blackbody_lambda(wl, temp).value
        temp = get_temp(dt)
        scale = flux[wl>4100][0]
        shifted = flux/scale-shift[ii]
        plot_spec(ax, wl/(1+z), shifted, tel, dt)

        smoothed = plot_smoothed_spec(
                ax, wl/(1+z), shifted, ivar, tel, dt)

        # Line identifications for dt=2
        if ii == 1:
            smoothed = plot_smoothed_spec(
                    ax, wl/(1+z), shifted, ivar, tel, dt, lw=1, c='#e55c30')
            sp = ["CIII", "OIII", "CIV", "OIV"]
            offsets = [1, 1, 1.1, 1.1]
            plot_species(ax, 45000, wl, smoothed, sp, offsets) 
        if ii == 3:
            sp = ["CIII", "OIII", "CIV", "OIV"]
            offsets = [0.8, 1, 1.1, 1.1]
            plot_species(ax, 35000, wl, smoothed, sp, offsets)
        if ii == 7:
            smoothed = plot_smoothed_spec(
                    ax, wl/(1+z), shifted, ivar, tel, dt, lw=1, c='#e55c30')
            sp = ["CIII", "OII", "CII", "SIV"]
            offsets = [0.95, 1, 1, 1]
            plot_species(ax, 30000, wl, smoothed, sp, offsets)
    ax.set_xlim(3000, 8440)
    ax.set_ylim(-2,1.5)
    handles, labels = ax.get_legend_handles_labels()
    by_label = OrderedDict(zip(labels,handles))
    ax.legend(
            by_label.values(), by_label.keys(), 
            loc='upper right', fontsize=12, ncol=3)
    ax.set_ylabel(
        r"Scaled $F_{\lambda}$ + cst.",
        fontsize=16)


def early_comparison(ax):
    """ Early spec compared to 18cow, 10vgv, 12gzk """
    files, epochs, tels = get_files(1, 2)
    f = files[0]
    tel = tels[0]
    dt = epochs[0]
    wl, flux, ivar = load_spec(f, tel)
    wl, flux = fluxcal(wl, flux, dt)
    wl, flux = clip_lines(wl, flux, z, tel, dt)
    wl, flux = clip_tellurics(wl, flux)
    scale = flux[wl>4100][0]
    #shifted = flux/scale-shift[ii]
    plot_spec(ax, wl/(1+z), flux/scale, tel, dt)
    smoothed = plot_smoothed_spec(
        ax, wl/(1+z), flux/scale, ivar, tel, dt, lw=1.0, text=False, 
        label='18gep, +1.0d, $T=%s$\,kK' %int(get_temp(1.0)/1000),
        c='#e55c30')
    plot_18cow(ax, 1.6E-14)
    plot_10vgv(ax, 0.6)
    plot_12gzk(ax, 1E-15)
    ax.set_ylim(0.05,1.8)
    ax.set_ylabel(
            r"Scaled $F_{\lambda}$ + cst.",
            fontsize=16)
    ax.legend(fontsize=12, loc='upper right')


def w_comparison(ax):
    """ Compare spec with W feature """
    files, epochs, tels = get_files(7,8)
    f = files[0]
    tel = tels[0]
    dt = epochs[0]
    wl, flux, ivar = load_spec(f, tel)
    wl, flux = fluxcal(wl, flux, dt)
    wl, flux = clip_lines(wl, flux, z, tel, dt)
    wl, flux = clip_tellurics(wl, flux)
    scale = flux[wl>4100][0]
    plot_spec(ax, wl/(1+z), flux/scale, tel, dt)
    smoothed = plot_smoothed_spec(
        ax, wl/(1+z), flux/scale, ivar, tel, dt, lw=1.0, text=False, 
        label='18gep, +4.2d, $T=20$\,kK', c='#e55c30')
    dat = np.loadtxt(SPEC_DIR + "/2008d.txt", delimiter=',')
    x = dat[:,0]
    y = dat[:,1]
    ext = fitzpatrick99(x+100, 0.63)
    ax.plot(
            x-100, y/0.1+ext-2.0, lw=0.5, c='#84206b', 
            label="SN2008D, +1.4d, $T=11$\,kK")
    dat = np.loadtxt(SPEC_DIR + "/ptf12dam.txt", delimiter=',')
    x = dat[:,0]
    y = dat[:,1]
    ax.plot(
            x-750, y/2+0.2, lw=1, c='#f6d746', ls='-', 
            label="12dam, -25d, $T=15$--20\,kK")
    ax.legend(fontsize=12, loc='upper right')
    ax.set_ylim(-0.4,2.5)
    ax.set_ylabel(
        r"Scaled $F_{\lambda}$ + cst.",
        fontsize=16)


def grbsn_comparison(ax):
    """ Compare later spec with GRB-SNe """
    files, epochs, tels = get_files(16,17)
    f = files[0]
    tel = tels[0]
    dt = epochs[0]
    wl, flux, ivar = load_spec(f, tel)
    wl, flux = fluxcal(wl, flux, dt)
    wl, flux = clip_lines(wl, flux, z, tel, dt)
    wl, flux = clip_tellurics(wl, flux)
    scale = flux[wl>4100][0]
    plot_spec(ax, wl/(1+z), flux/scale, tel, dt)
    smoothed = plot_smoothed_spec(
        ax, wl/(1+z), flux/scale, ivar, tel, dt, lw=1.0, text=False, 
        label='18gep, +18d', c='#e55c30')
    plot_98bw(ax, 1E-14)
    plot_06aj(ax, 1E-16)
    ax.legend(fontsize=12, loc='upper left')
    ax.set_xlabel(r"Rest Wavelength (\AA)", fontsize=16)
    ax.set_ylabel(
            r"Scaled $F_{\lambda}$ + cst.",
            fontsize=16)
    

if __name__=="__main__":
    # Top panel as separate figure
    fig,ax = plt.subplots(1,1,figsize=(7,4))
    spec_evol(ax)
    plt.savefig("early_spec_evol.eps", format='eps', dpi=300, bbox_inches='tight')
    plt.close()

    # Bottom panels together
    fig,axarr = plt.subplots(
            3, 1, figsize=(7,8), sharex=True, 
            gridspec_kw={'height_ratios':[1,1,1]})

    ax = axarr[0]
    early_comparison(ax)
     
    ax = axarr[1]
    w_comparison(ax)

    ax = axarr[2]
    grbsn_comparison(ax)

    for ax in axarr:
        ax.tick_params(axis='both', labelsize=14)
        ax.get_yaxis().set_ticks([])

    plt.subplots_adjust(hspace=0.1)
    plt.savefig("early_spec_comparison.eps", format='eps', dpi=500, bbox_inches='tight')
    #plt.show()
    plt.close()
