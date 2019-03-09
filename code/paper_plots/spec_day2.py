""" Comparison of early spectrum to early 18cow spectrum
to make a point about the broad absorption feature """

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15
import glob
from spec_sequence import *


fig,ax = plt.subplots(1,1,figsize=(8,4))

gepfile = "ZTF18abukavn_20180911_P200_v2.ascii"
#gepfile = "ZTF18abukavn_20180909_LT_v1.ascii"
gep_dir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/spec/ZTF18abukavn"
wl, flux, ivar = load_spec(gep_dir + "/%s" %gepfile, 'P200')
wl, flux = fluxcal(wl, flux, 2.07)
wl, flux = clip_lines(wl, flux, 0.03154, 'P200', 2.07)
wl, flux = clip_tellurics(wl, flux)
plot_spec(ax, wl, flux/1E-15, 'P200', 2.07)
plot_smoothed_spec(
        ax, wl, flux/1E-15, ivar, 'P200', 2.07)

ax.text(0.9, 0.9, r"SN2018gep at $\Delta t=2.07$d (P200)", fontsize=14,
        transform=ax.transAxes, horizontalalignment='right',
        verticalalignment='top')
ax.text(0.9, 0.8, r"$T_\mathrm{eff}\approx30,000$\,K", fontsize=14,
        transform=ax.transAxes, horizontalalignment='right',
        verticalalignment='top')
ax.text(0.9, 0.7, r"$v\approx35,000\,$km/s", fontsize=14,
        transform=ax.transAxes, horizontalalignment='right',
        verticalalignment='top')

plt.tick_params(labelsize=14)
plt.xlim(3200, 9320)
plt.ylim(0.015, 8.7)
plt.xlabel(r"Observed Wavelength (\AA)", fontsize=16)
plt.ylabel(
    r"Flux $f_{\lambda}$ ($10^{-15}$ erg cm$^{-2}$ s$^{-1}$ \AA$^{-1}$)", 
    fontsize=16)
plt.legend(fontsize=14)

plt.tight_layout()
#plt.savefig("spec_day2.png")
plt.show()

