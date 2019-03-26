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


def full_lc():
    """ Plot the full LC in g and r, showing all the progenitor detections """
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
            (code=='final photometry', dt < -1, ~np.isnan(mag)))
    prog_nondet = np.logical_and(code=='final photometry', np.isnan(mag))

    # Initialize the figure
    fig,ax = plt.subplots(1,1,figsize=(10,3))

    # Plot the g-band LC
    choose = np.logical_and(sn_det, filt=='ztfg')
    ax.errorbar(
            dt[choose], mag[choose], yerr=emag[choose], 
            c='k', ms=5, fmt='s', label="P48 $g$", zorder=2)

    # Plot the r-band detections
    choose = np.logical_and(sn_det, filt=='ztfr')
    ax.errorbar(
            dt[choose], mag[choose], yerr=emag[choose], 
            ms=5, fmt='o', mfc='white', mec='grey', label="P48 $r$", c='grey',
            zorder=0)

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
            c='k', ms=5, fmt='s', label=None, zorder=2)

    # Plot the g-band prog non-detections
    choose = np.logical_and(prog_nondet, filt=='ztfg')
    ax.scatter(
            dt[choose], limmag[choose], 
            color='k', 
            s=20, marker='_', label=None, zorder=2)

    # Plot the r-band prog LC
    choose = np.logical_and(prog_det, filt=='ztfr')
    ax.errorbar(
            dt[choose], mag[choose], yerr=emag[choose], 
            ms=5, fmt='o', mfc='white', mec='grey', label=None, c='grey',
            zorder=0)

    # Plot the r-band prog non-detections
    choose = np.logical_and(prog_nondet, filt=='ztfr')
    ax.scatter(
            dt[choose], limmag[choose], 
            color='grey', 
            s=10, marker='_', label=None, zorder=2)

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

    # Format this box
    ax.set_xlim(-175, 33)
    ax.set_ylabel(r"Apparent Mag", fontsize=16)
    ax.set_xlabel("Days since $t_0$", fontsize=16)
    ax.yaxis.set_tick_params(labelsize=14)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.invert_yaxis()
    ax.legend(loc='lower right', fontsize=14)

    plt.tight_layout()
    plt.show()
    #plt.savefig("lc_full.png")


def lc_zoom():
    """ Plot the zoomed LC in luminosity """
    DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/phot"
    f = DATA_DIR + "/ZTF18abukavn_opt_phot.dat"
    dat = np.loadtxt(f, dtype=str, delimiter=' ')
    instr = dat[:,0]
    jd = dat[:,1].astype(float)
    dt = jd-t0
    filt = dat[:,2]
    mag = dat[:,3].astype(float)
    emag = dat[:,4].astype(float)
    det = np.logical_and(mag<99, ~np.isnan(mag))
    lum, elum = mag_to_flux(mag, emag)

    # Progenitor data
    f = DATA_DIR + "/precursor.csv"
    prog = ascii.read(f)
    prog_mjd = prog['mjd']
    prog_jd = prog_mjd + 2400000.5
    prog_dt = prog_jd - t0
    prog_filter = prog['filter']
    prog_mag = prog['mag']
    prog_emag = prog['magerr']
    prog_limmag = prog['lim_mag']
    code = prog['instrument']
    prog_det = np.logical_and(prog_dt < 0, ~np.isnan(prog_mag))
    prog_nondet = np.logical_and(code=='ZTF Deep Stacks', np.isnan(prog_mag))
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
            zorder=0)

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
    ax.set_xlabel("Days since $t_0$", fontsize=16)
    ax.set_yscale('log')
    ax.yaxis.set_tick_params(labelsize=14)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.legend(loc='lower right', fontsize=14)

    plt.tight_layout()
    plt.show()
 

if __name__=="__main__":
    full_lc()
