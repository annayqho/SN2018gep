""" 
Plot the full early LC, including progenitor stuff
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.cosmology import Planck15
from astropy.io import ascii

z = 0.03154
d = Planck15.luminosity_distance(z=0.03154).cgs.value
# JD of first (r-band) detection
t0 = 2458370.6634

def mag_to_flux(mag, emag):
    """ Convert from mag to flux """
    # AB flux zero point for g-band is 3631 Jy or 1E-23 erg/cm2/s/Hz
    flux = 10**(-(mag + 48.6)/2.5)
    lum = 4 * np.pi * d**2 * flux
    # Uncertainty in magnitude is roughly the fractional uncertainty on the flux
    eflux = emag*flux
    elum = 4 * np.pi * d**2 * eflux
    return np.array(lum), np.array(elum)


def get_data():
    DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/phot"

    # Full light curve from Danny
    f = DATA_DIR + "/ZTF18abukavn.csv"
    dat = ascii.read(f)
    mjd = dat['mjd']
    jd = mjd + 2400000.5
    dt = jd - t0
    filt = dat['filter']
    mag = dat['mag']
    emag = dat['magerr']
    limmag = dat['lim_mag']
    code = dat['instrument']
    sn_det = np.logical_and.reduce(
            (code=='ZTF Camera', dt > -1, ~np.isnan(mag)))
    prog_det = np.logical_and.reduce(
            (code=='final photometry', dt < 0, ~np.isnan(mag)))
    prog_nondet = np.logical_and(code=='final photometry', np.isnan(mag))
    return dt, filt, mag, emag, limmag, sn_det, prog_det, prog_nondet


def full_lc(figx, figy, xmin, xmax, ymin, ymax, figname, annotations=False):
    """ Plot the full LC in g and r, showing all the progenitor detections """
    dt, filt, mag, emag, limmag, sn_det, prog_det, prog_nondet = get_data()

    # Initialize the figure
    fig,ax = plt.subplots(1,1,figsize=(figx, figy))

    # Plot the g-band LC
    choose = np.logical_and(sn_det, filt=='ztfg')
    gcol = '#140b34'
    ax.errorbar(
            dt[choose], mag[choose], yerr=emag[choose], 
            c=gcol, ms=5, fmt='s', label="P48 $g$, single images", 
            zorder=3, lw=0.5)

    # Plot the r-band LC
    choose = np.logical_and(sn_det, filt=='ztfr')
    rcol = '#e55c30'
    ax.errorbar(
            dt[choose], mag[choose], yerr=emag[choose], 
            ms=5, fmt='o', mfc=rcol, mec=rcol, label="P48 $r$, single images", 
            c=rcol, zorder=2, lw=0.5)

    # Plot the r-band non-detection, 32.5 min before the first det
    # rband = np.logical_and(instr=='P48+ZTF', filt=='r')
    # choose = np.logical_and(~det, rband)
    # xnondet = -32.5/60/24
    # ynondet = 21.25
    # ax.scatter(
    #         xnondet, ynondet, color='grey', marker='_', label=None, s=20)
    # ax.arrow(
    #         xnondet, ynondet, 0, 0.7, length_includes_head=True,
    #         head_width=1, head_length=0.2, fc='grey', ec='grey')

    # Plot the g-band prog LC
    choose = np.logical_and(prog_det, filt=='ztfg')
    ax.errorbar(
            dt[choose], mag[choose], yerr=emag[choose], 
            c=gcol, mec=gcol, mfc='white', ms=5, fmt='s', zorder=1,
            label="P48 $g$, 3-day stacks")

    # Plot the g-band prog non-detections
    choose = np.logical_and(prog_nondet, filt=='ztfg')
    ax.scatter(
            dt[choose], limmag[choose], 
            color=gcol, 
            s=30, marker='_', label=None, zorder=0, lw=1.0)

    # Plot the r-band prog LC
    choose = np.logical_and(prog_det, filt=='ztfr')
    ax.errorbar(
            dt[choose], mag[choose], yerr=emag[choose], 
            ms=5, fmt='o', mfc='white', mec=rcol, c=rcol,
            zorder=1, lw=0.5,
            label="P48 $r$, 3-day stacks")

    # Plot the r-band prog non-detections
    choose = np.logical_and(prog_nondet, filt=='ztfr')
    ax.scatter(
            dt[choose], limmag[choose], 
            color=rcol, 
            s=30, marker='_', label=None, zorder=0, lw=1.0)

    ax.set_ylim(ymin, ymax)
    # Add an axis on the right-hand side showing the absolute mag
    ax2 = ax.twinx()
    ax2.set_ylabel(
            "Absolute Mag",
            fontsize=16, rotation=270, labelpad=15.0)
    y_f = lambda y_i: y_i - Planck15.distmod(z=z).value
    ymin, ymax = ax.get_ylim()
    ax2.set_ylim((y_f(ymin), y_f(ymax)))
    ax2.plot([],[])
    ax2.tick_params(axis='both', labelsize=14)
    ax2.invert_yaxis()

    if annotations:
        # a bunch of stuff showing our follow-up timeline
        textor = 'vertical' # textorientation
        ax.axvline(x=0.48, ymin=0, ymax=0.2, lw=0.5, c='grey') # UVOT
        ax.text(
                0.48, 20.2, 'Swift', fontsize=12, horizontalalignment='center',
                rotation=textor)
        ax.axvline(x=0.7, ymin=0, ymax=0.05, lw=0.5, c='grey') # LT
        ax.text(
                0.73, 22, 'LT', fontsize=12, horizontalalignment='center',
                rotation=textor)
        ax.axvline(x=1.0, ymin=0, ymax=0.1, lw=0.5, c='grey') # P200/P60
        ax.text(
                1.0, 19.7, 'P200, P60', fontsize=12, horizontalalignment='center',
                rotation=textor) 
        ax.axvline(x=1.7, ymin=0, ymax=0.05, lw=0.5, c='grey') # LT
        ax.text(1.7, 22, 'LT', fontsize=12, horizontalalignment='center',
                rotation=textor)
        ax.axvline(x=2.0, ymin=0, ymax=0.2, lw=0.5, c='grey') # P200
        ax.text(2.0, 20.3, 'P200', fontsize=12, horizontalalignment='center',
                rotation=textor)

        # and an inset
        axins = inset_axes(
                ax, 2, 1, loc=1,
                bbox_to_anchor=(0.57,0.98),
                bbox_transform=ax.transAxes) 

        # Plot the g-band LC
        choose = np.logical_and(sn_det, filt=='ztfg')
        gcol = '#140b34'
        axins.errorbar(
                dt[choose]*24*60, mag[choose], yerr=emag[choose], 
                c=gcol, ms=5, fmt='s', label="P48 $g$, single images", 
                zorder=3, lw=0.5)

        # Plot the r-band LC
        choose = np.logical_and(sn_det, filt=='ztfr')
        rcol = '#e55c30'
        axins.errorbar(
                dt[choose]*24*60, mag[choose], yerr=emag[choose], 
                ms=5, fmt='o', mfc=rcol, mec=rcol, label="P48 $r$, single images", 
                c=rcol, zorder=2, lw=0.5)

        # Format the inset
        axins.set_xlim(-50,160)
        axins.set_ylim(18,21.5)
        axins.invert_yaxis()
        axins.tick_params(axis='both', labelsize=12)
        axins.set_ylabel(r"Apparent Mag", fontsize=12)
        axins.set_xlabel(r"Minutes since ZTF discovery", fontsize=12)

        axins.text(
                0.05, 0.8, "1.3 mag/hr", fontsize=12, 
                transform=axins.transAxes)
        axins.axvline(x=-23, ls='--', lw=0.5)
    
    else:
        ax.legend(loc='upper left', fontsize=14)

    # Format this box
    ax.set_xlim(xmin, xmax)
    ax.set_ylabel(r"Apparent Mag", fontsize=16)
    ax.set_xlabel("Days since ZTF discovery", fontsize=16)
    ax.yaxis.set_tick_params(labelsize=14)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.invert_yaxis()

    plt.tight_layout()
    plt.show()
    #plt.savefig(figname)


def lc_fit():
    """ Plot the zoomed LC in luminosity """
    dt, filt, mag, emag, limmag, sn_det, prog_det, prog_nondet = get_data()
    lum, elum = mag_to_flux(mag, emag)
    progf, progef = mag_to_flux(prog_mag, prog_emag)

    # Initialize the figure
    fig,ax = plt.subplots(1,1,figsize=(8,3))

    # Plot the g-band LC
    gband = np.logical_and(instr=='P48+ZTF', filt=='g')
    choose = np.logical_and(det, gband)
    ax.errorbar(
            dt[choose], lum[choose]/1E28, yerr=elum[choose]/1E28, 
            c='k', ms=5, fmt='s', label="P48 $g$", zorder=2)

    # Plot the r-band detections
    rband = np.logical_and(instr=='P48+ZTF', filt=='r')
    choose = np.logical_and(det, rband)
    ax.errorbar(
            dt[choose], lum[choose]/1E28, yerr=elum[choose]/1E28, 
            ms=5, fmt='o', mfc='white', mec='grey', label="P48 $r$", c='grey',
            zorder=2)

    # Plot the g-band prog LC
    choose = np.logical_and(prog_det, prog_filter=='ztfg')
    ax.errorbar(
            prog_dt[choose], progf[choose]/1E28, yerr=progef[choose]/1E28, 
            c='k', ms=5, fmt='s', label=None, zorder=2)

    # Plot the g-band prog non-detections
    # choose = np.logical_and(prog_nondet, prog_filter=='ztfg')
    # ax.scatter(
    #         prog_dt[choose], prog_limmag[choose], 
    #         color='k', 
    #         s=20, marker='_', label=None, zorder=2)

    # Plot the r-band prog LC
    #choose = np.logical_and(prog_det, prog_filter=='ztfr')
    #ax.errorbar(
    #        prog_dt[choose], progf[choose], yerr=progef[choose], 
    #        ms=5, fmt='o', mfc='white', mec='grey', label=None, c='grey',
    #        zorder=0)
# 
    # Plot the r-band prog non-detections
#     choose = np.logical_and(prog_nondet, prog_filter=='ztfr')
#     ax.scatter(
#             prog_dt[choose], prog_limmag[choose], 
#             color='grey', 
#             s=10, marker='_', label=None, zorder=2)
# 
    # Format this box
    ax.set_xlim(-10, 3)
    ax.set_ylabel(r"$L_\nu$ [$10^{28}$ erg/s/Hz]", fontsize=16)
    ax.set_xlabel("Days since ZTF discovery", fontsize=16)
    ax.set_yscale('log')
    ax.yaxis.set_tick_params(labelsize=14)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.legend(loc='lower right', fontsize=14)

    plt.tight_layout()
    plt.show()
 

if __name__=="__main__":
    #full_lc(10, 3, -175, 35, 15.5, 23, "full_gr.png")
    full_lc(10, 3, -18, 2.2, 15.5, 22.5, "zoom_gr.png", annotations=True)
