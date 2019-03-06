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

d = Planck15.luminosity_distance(z=0.033).cgs.value
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
    return lum, elum


def full_lc():
    """ Plot the full LC in g and r, showing all the progenitor detections """
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
    prog_nondet = np.logical_and(prog_dt < 0, np.isnan(prog_mag))

    # Initialize the figure
    fig,ax = plt.subplots(1,1,figsize=(8,4))

    # Plot the g-band LC
    gband = np.logical_and(instr=='P48+ZTF', filt=='g')
    choose = np.logical_and(det, gband)
    ax.errorbar(
            dt[choose], mag[choose], yerr=emag[choose], 
            c='k', ms=5, fmt='s', label="P48 $g$", zorder=2)

    # Plot the r-band detections
    rband = np.logical_and(instr=='P48+ZTF', filt=='r')
    choose = np.logical_and(det, rband)
    ax.errorbar(
            dt[choose], mag[choose], yerr=emag[choose], 
            ms=5, fmt='o', mfc='white', mec='grey', label="P48 $r$", c='grey',
            zorder=0)

    # Plot the g-band prog LC
    choose = np.logical_and(prog_det, prog_filter=='ztfg')
    ax.errorbar(
            prog_dt[choose], prog_mag[choose], yerr=prog_emag[choose], 
            c='k', ms=5, fmt='s', label=None, zorder=2)

    # Plot the g-band prog non-detections
    choose = np.logical_and(prog_nondet, prog_filter=='ztfg')
    for ii,dt_val in enumerate(prog_dt[choose]):
        ax.arrow(
            dt_val, prog_limmag[choose][ii], 0, 0.5, length_includes_head=True,
            head_width=0.01, head_length=0.1, fc='k', ec='k') 

    # Plot the r-band prog LC
    choose = np.logical_and(prog_det, prog_filter=='ztfr')
    ax.errorbar(
            prog_dt[choose], prog_mag[choose], yerr=prog_emag[choose], 
            ms=5, fmt='o', mfc='white', mec='grey', label=None, c='grey',
            zorder=0)



# # Show the last non-detection
# ax.axvline(x=-32.5, ls=':', c='grey')
# ax.text(-32.5, 19, 'ND', fontsize=14,
#         horizontalalignment='center',
#         verticalalignment='center')
# ax.axvline(x=-23, ls='--', c='k')
# ax.text(-18, 19, '$t_0$', fontsize=16, horizontalalignment='center',
#         verticalalignment='center')

    # Format this box
    #ax.set_xlim(-40, 80)
    #ax.set_ylim(18.5,21)
    ax.set_ylabel(r"Apparent Mag", fontsize=16)
    ax.set_xlabel("Days since $t_0$", fontsize=16)
    ax.yaxis.set_tick_params(labelsize=14)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.invert_yaxis()
    ax.legend(loc='lower right', fontsize=14)

    plt.tight_layout()
    plt.show()


if __name__=="__main__":
    full_lc()
