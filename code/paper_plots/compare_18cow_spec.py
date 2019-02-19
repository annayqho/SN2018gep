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


fig = plt.figure(figsize=(8,4))

cow_dir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/spec/AT2018cow"
wl_cow, flux_cow, ivar_cow = load_spec(
        cow_dir + "/AT2018cow_20180621_P200_v3.ascii", 'P200')

gepfile = "ZTF18abukavn_20180910_P200_v6.ascii"
#gepfile = "ZTF18abukavn_20180909_LT_v1.ascii"
gep_dir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/spec/ZTF18abukavn"
wl_gep, flux_gep, ivar_gep = load_spec(gep_dir + "/%s" %gepfile, 'P200')
wl_gep, flux_gep = clip_lines(wl_gep, flux_gep, 0.03154, 'P200', 1.0)

plt.plot(
        wl_cow, flux_cow/4.7/1E-15,
        c='grey', drawstyle='steps-mid', lw=0.5, 
        label="AT2018cow/4.7, $\Delta t=5.353$")
plt.plot(
        wl_gep, flux_gep/1E-15,
        c='k', drawstyle='steps-mid', lw=0.5,
        label="AT2018gep, $\Delta t=1.0$")

plt.tick_params(labelsize=14)
plt.xlim(3700, 7000)
plt.ylim(0.3, 4.1)
plt.xlabel(r"Observed Wavelength (\AA)", fontsize=16)
plt.ylabel(
    r"Flux $f_{\lambda}$ ($10^{-15}$ erg cm$^{-2}$ s$^{-1}$ \AA$^{-1}$)", 
    fontsize=16)
plt.legend(fontsize=14)

plt.tight_layout()
plt.savefig("spec_comp_18cow.png")
#plt.show()

