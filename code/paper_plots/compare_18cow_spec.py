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

cow_dir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/spec/AT2018cow"
dt = 5.353
gepfile = "ZTF18abukavn_20180910_P200_v6.ascii"
gep_dir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/spec/ZTF18abukavn"
tel = 'P200'
dt = 1.0
wl_gep, flux_gep, ivar_gep = load_spec(gep_dir + "/%s" %gepfile, tel)
wl, flux = fluxcal(wl_gep, flux_gep, dt)
wl, flux = clip_lines(wl, flux, 0.03154, 'P200', dt)
wl, flux = clip_tellurics(wl, flux)
plot_spec(ax, wl, flux/1E-15, 'P200', dt)
plot_smoothed_spec(
        ax, wl, flux/1E-15, ivar_gep, 'P200', dt, lw=1.0, 
        label="SN2018gep, $\Delta t=1.0$")
ax.legend(fontsize=14, loc='upper right')

wl_cow, flux_cow, ivar_cow = load_spec(
        cow_dir + "/AT2018cow_20180621_P200_v3.ascii", 'P200')
plot_smoothed_spec(
        ax, wl_cow, flux_cow/1E-15/4.7, ivar_cow, 
        'P200', dt, lw=1, ls='--', label="AT2018cow/4.7, $\Delta t=5.353$")

# Add suggestions for lines
ax.scatter(3490, 5.2, marker='v', c='k')
ax.text(
        3490, 5.3, "CII", fontsize=12, 
        verticalalignment='bottom', horizontalalignment='center')
ax.scatter(3800, 4.2, marker='v', c='k')
ax.text(
        3800, 4.3, "CII", fontsize=12, 
        verticalalignment='bottom', horizontalalignment='center')
ax.scatter(4141, 3.5, marker='v', c='k')
ax.text(
        4141, 3.6, "CIII", fontsize=12, 
        verticalalignment='bottom', horizontalalignment='center')
ax.scatter(5072, 1.6, marker='v', c='k')
ax.text(
        5072, 1.7, "CIII", fontsize=12, 
        verticalalignment='bottom', horizontalalignment='center')
ax.scatter(3319, 4.6, marker='v', c='k')
ax.text(
        3319, 4.7, "OII", fontsize=12, 
        verticalalignment='bottom', horizontalalignment='center')

plt.tick_params(labelsize=14)
plt.xlim(3190, 9000)
plt.ylim(0.25, 5.7)
plt.xlabel(r"Observed Wavelength (\AA)", fontsize=16)
plt.ylabel(
    r"Flux $f_{\lambda}$ ($10^{-15}$ erg cm$^{-2}$ s$^{-1}$ \AA$^{-1}$)", 
    fontsize=16)
plt.legend(fontsize=14)

plt.tight_layout()
#plt.savefig("spec_comp_18cow.png")
plt.show()

