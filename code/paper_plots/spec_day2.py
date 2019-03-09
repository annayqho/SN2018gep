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

#plt.plot(
#        wl, flux/1E-15, c='k', drawstyle='steps-mid', lw=0.5)

plt.tick_params(labelsize=14)
plt.xlim(3700, 7000)
plt.ylim(0.3, 4.1)
plt.xlabel(r"Observed Wavelength (\AA)", fontsize=16)
plt.ylabel(
    r"Flux $f_{\lambda}$ ($10^{-15}$ erg cm$^{-2}$ s$^{-1}$ \AA$^{-1}$)", 
    fontsize=16)
plt.legend(fontsize=14)

plt.tight_layout()
#plt.savefig("spec_day2.png")
plt.show()

