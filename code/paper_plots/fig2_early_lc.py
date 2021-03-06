""" 
Fit for the explosion epoch by fitting a polynomial 
to the g-band light curve """

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.gridspec as gridspec
from astropy.cosmology import Planck15
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/code")
from load_lc import get_forced_phot


d = Planck15.luminosity_distance(z=0.03154).cgs.value


def mag2lum(mag):
    """ 
    Convert from mag to lum
    AB flux zero point for g-band is 3631 Jy or 1E-23 erg/cm2/s/Hz
    """
    flux = 10**(-(mag + 48.6)/2.5)
    lum = 4 * np.pi * d**2 * flux
    # Uncertainty in magnitude is roughly the fractional uncertainty on the flux
    return lum


def emag2elum(mag,emag):
    lum = mag2lum(mag)
    elum_top = np.abs(mag2lum(mag+emag) - lum)
    elum_bottom = np.abs(mag2lum(mag-emag) - lum)
    elum = np.max([elum_top, elum_bottom], axis=0)
    return elum


def get_t0():
    jd, filt, mag, emag, limmag, code = get_forced_phot()
    # first detection was 2458370.6634
    dt = jd-2458370.6634
    use_filter = 'ztfg'
    # use_filter = 'ztfr' # for referee report
    window = 3
    choose = np.logical_and.reduce(
            (code=='ZTF Camera', filt==use_filter, ~np.isnan(mag)))
    sec = np.logical_and(dt[choose] < window, dt[choose] > 0)

    x = dt[choose][sec]
    y = mag[choose][sec]
    ey = emag[choose][sec]

    npts = len(x)
    ndraw = 10000
    # make 10000 new versions of the y array
    y_new = np.zeros((npts, ndraw))
    for ii,yval in enumerate(y):
        y_new[ii,:] = np.random.normal(loc=yval, scale=ey[ii], size=ndraw)

    # For each version, fit for t0
    t0s = np.zeros(ndraw)
    for ii in np.arange(ndraw):
        yconv = mag2lum(y_new[:,ii])
        out = np.polyfit(x, yconv, deg=2)
        a,b,c = out
        t0s[ii] = (-b + np.sqrt(b**2-4*a*c))/(2*a)

    t0 = np.mean(t0s)
    et0 = np.std(t0s)

    return t0, et0, (a,b,c)


def plot_firstmin_flux(ax):
    """ Plot the first few minutes """
    #dt, mag, emag, instr, filt, det = load_lc()
    jd, filt, mag, emag, limmag, code = get_forced_phot()
    t0 = 2458370.6634
    dt = jd-t0
    choose = np.logical_and.reduce(
            (code=='ZTF Camera', filt=='ztfg', ~np.isnan(mag)))
    lum = mag2lum(mag[choose])
    elum = emag2elum(mag[choose], emag[choose])
    ax.errorbar(
            dt[choose]*24*60, lum/1E26, yerr=elum/1E26, 
            c='k', ms=5, fmt='s', label="P48 $g$", zorder=2)

    # Plot the r-band detections
    choose = np.logical_and.reduce(
            (code=='ZTF Camera', filt=='ztfr', ~np.isnan(mag)))
    lum = mag2lum(mag[choose])
    elum = emag2elum(mag[choose], emag[choose])
    ax.errorbar(
            dt[choose]*24*60, lum/1E26, yerr=elum/1E26, ms=5, fmt='o', 
            mfc='#e55c30', mec='#e55c30', label="P48 $r$", c='#e55c30',
            zorder=0)

    t0,et0,fitparams = get_t0()
    xm = np.linspace(-50,150)
    xd = xm/(60*24)
    yd = fitparams[0]*xd**2 + fitparams[1]*xd + fitparams[2]
    ax.plot(xm, yd/1E26, c='k', lw=0.5, ls='--')
    ax.axvspan(-24.83-2, -24.83+2, lw=0.5, color='lightgrey')
    ax.text(-20, 40, "$t_0=-25\pm2$ min", fontsize=14,
            verticalalignment='top')

    # Format this box
    ax.set_xlim(-50, 150)
    ax.set_ylim(0, 43)
    ax.set_ylabel(r"$L_\nu$ [$10^{26}$ erg/s/Hz]", fontsize=16)
    ax.set_xlabel("Minutes since ZTF discovery", fontsize=16)
    ax.yaxis.set_tick_params(labelsize=14)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.legend(loc='lower right', fontsize=14)


def plot_firstdays_flux(ax):
    """ Plot the first few days """
    #dt, mag, emag, instr, filt, det = load_lc()
    jd, filt, mag, emag, limmag, code = get_forced_phot()
    t0 = 2458370.6634
    dt = jd-t0
    lum = mag2lum(mag)
    elum = emag2elum(mag, emag)
    choose = np.logical_and.reduce(
            (code=='ZTF Camera', filt=='ztfg', ~np.isnan(mag)))
    ax.errorbar(
            dt[choose], lum[choose]/1E28, yerr=elum[choose]/1E28, 
            c='k', ms=5, fmt='s', label="P48 $g$", zorder=2)
    #out,cov = get_fit_func()
    #xlab = np.linspace(-1,2)
    #ylab = out[0]*xlab**2 + out[1]*xlab + out[2]
    #ax.plot(xlab*24*60, ylab/1E28, c='k', ls='--')

    # Plot the r-band detections
    choose = np.logical_and.reduce(
            (code=='ZTF Camera', filt=='ztfr', ~np.isnan(mag)))
    ax.errorbar(
            dt[choose], lum[choose]/1E28, yerr=elum[choose]/1E28, 
            ms=5, fmt='o', mfc='#e55c30', mec='#e55c30', label="P48 $r$", 
            c='#e55c30', zorder=0)

    # plot the fit
    t0,et0,fitparams = get_t0()
    xd = np.linspace(-0.1, 2.1)
    yd = fitparams[0]*xd**2 + fitparams[1]*xd + fitparams[2]
    ax.plot(xd, yd/1E28, c='k', lw=0.5, ls='--')

    # Format this box
    ax.set_xlim(-0.1, 2.1)
    ax.set_ylim(-0.5, 4)
    ax.set_ylabel(r"$L_\nu$ [$10^{28}$ erg/s/Hz]", fontsize=16)
    ax.set_xlabel("Days since ZTF discovery", fontsize=16)
    ax.yaxis.set_tick_params(labelsize=14)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.legend(loc='lower right', fontsize=14)


def plot_firstmin_mag(ax):
    """ Plot the first few minutes """
    jd, filt, mag, emag, limmag, code = get_forced_phot()
    t0 = 2458370.6634
    dt = jd-t0
    choose = np.logical_and.reduce(
            (code=='ZTF Camera', filt=='ztfg', ~np.isnan(mag)))
    ax.errorbar(
            dt[choose]*24*60, mag[choose], yerr=emag[choose], 
            c='k', ms=5, fmt='s', label="P48 $g$", zorder=2)
    ax.text(63, 18.9, "$1.4 \pm 0.1$ mag/hr", fontsize=14,
            verticalalignment='top', horizontalalignment='right')

    # Plot the r-band detections
    choose = np.logical_and.reduce(
            (code=='ZTF Camera', filt=='ztfr', ~np.isnan(mag)))
    ax.errorbar(
            dt[choose]*24*60, mag[choose], yerr=emag[choose], 
            ms=5, fmt='o', mfc='#e55c30', mec='#e55c30', 
            label="P48 $r$", c='#e55c30', zorder=0)

    # Format this box
    ax.set_xlim(-50, 150)
    ax.set_ylim(18, 21.5)
    ax.invert_yaxis()
    ax.set_ylabel(r"Apparent Mag", fontsize=16)
    #ax.set_xlabel("Minutes since first detection", fontsize=16)
    ax.yaxis.set_tick_params(labelsize=14)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.legend(loc='lower right', fontsize=14)


def plot_full_lc(ax):
    """ plot the full LC with pre-bump """
    dt, mag, emag, instr, filt, det = load_lc()
    lum,elum = mag2lum(mag,emag)
    gband = np.logical_and(instr=='P48+ZTF', filt=='g')
    choose = np.logical_and(det, gband)
    ax.errorbar(
            dt[choose], lum[choose]/1E28, yerr=elum[choose]/1E28, 
            c='k', ms=5, fmt='s', label="P48 $g$", zorder=2)
    rband = np.logical_and(instr=='P48+ZTF', filt=='r')
    choose = np.logical_and(det, rband)
    ax.errorbar(
            dt[choose], lum[choose]/1E28, yerr=elum[choose]/1E28, 
            ms=5, fmt='o', mfc='white', mec='grey', label="P48 $r$", c='grey',
            zorder=0)

    out,cov = get_fit_func()
    xlab = np.linspace(-1,2)
    ylab = out[0]*xlab**2 + out[1]*xlab + out[2]
    ax.plot(xlab, ylab/1E28, c='k', ls='--')

    # Vertical line for the first UVOT epoch
    textor = 'vertical' # textorientation
    ax.axvline(x=0.48, lw=2, c='lightblue') # UVOT
    ax.text(
            0.48, 4.2, 'Swift', fontsize=14, horizontalalignment='right',
            rotation=textor)
    ax.axvline(x=0.7, lw=2, c='lightblue') # LT
    ax.text(
            0.7, 4.2, 'LT', fontsize=14, horizontalalignment='right',
            rotation=textor)
    ax.axvline(x=1.0, lw=2, c='lightblue') # P200/P60
    ax.text(
            1.0, 0.8, 'P200', fontsize=14, horizontalalignment='right',
            rotation=textor)
    ax.text(1.0, 4.2, 'P60', fontsize=14, horizontalalignment='right',
            rotation=textor)
    ax.axvline(x=1.7, lw=2, c='lightblue') # LT
    ax.text(1.7, 4.2, 'LT', fontsize=14, horizontalalignment='right',
            rotation=textor)
    ax.axvline(x=2.0, lw=2, c='lightblue') # P200
    ax.text(2.0, 2.0, 'P200', fontsize=14, horizontalalignment='right',
            rotation=textor)
    ax.axvline(x=2.7, lw=2, c='lightblue') # LT
    ax.text(2.7, 4.2, 'LT', fontsize=14, horizontalalignment='right',
            rotation=textor)

    ax.legend(loc='lower right', fontsize=14)
    ax.set_ylabel(r"$L_\nu$ [$10^{28}$ erg/s/Hz]", fontsize=16)
    ax.set_xlabel("Days since ZTF discovery", fontsize=16)
    ax.yaxis.set_tick_params(labelsize=14)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.set_xlim(-0.2, 2.2)
    ax.set_ylim(-0.5,4.7)


if __name__=="__main__":
    # Initialize the figure
    # fig = plt.figure(figsize=(8,5))
    # gs = fig.add_gridspec(
    #         nrows=2, ncols=2, wspace=0.05)
    # plot_firstmin_mag(fig.add_subplot(gs[0,0]))
    # plot_firstmin_flux(fig.add_subplot(gs[1,0]))
    # plot_firstdays_flux(fig.add_subplot(gs[:,1]))

    # gs.tight_layout(fig)
    # #plt.savefig("first_days.eps", format='eps', dpi=500)
    # plt.show()

    fig,ax = plt.subplots(1,1,figsize=(8,5))
    plot_full_lc(ax)
