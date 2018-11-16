""" Plot the spectral sequence """

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15
import glob

# We have to apply offsets, because the spectra are too close to each other

fig,ax = plt.subplots(1,1,figsize=(6,10))

files = glob.glob("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/spec/ZTF18abukavn/*.ascii")

t0 = 2458370.6634 # in JD
for f in files:
    # In Dan's 18cow paper, he interpolates over host narrow features
    dat = np.loadtxt(f)
    wl = dat[:,0]
    flux_raw = dat[:,1]

    # If the spectra are from LT, then you multiply by 1E-15?
    tel = f.split("_")[2]
    if tel == 'LT':
        flux = flux_raw*1E-15
    elif tel == 'P200':
        flux = flux_raw
    elif tel == 'Keck1':
        flux = flux_raw
    elif tel == 'DCT':
        flux = flux_raw
    elif tel == 'NOT':
        flux = flux_raw
    elif tel == 'P60':
        flux = flux_raw

    # Plot the spectrum
    ax.plot(wl, flux, c='k', alpha=0.7, drawstyle='steps-mid')

ax.set_ylabel(
    r"Flux $f_{\lambda}$ (erg cm${}^{-2}$ s$^{-1}$ \AA$^{-1}$)", fontsize=16)
ax.set_xlabel(
        r"Observed Wavelength (\AA)", fontsize=16)
ax.yaxis.set_tick_params(labelsize=14)
ax.xaxis.set_tick_params(labelsize=14)
#ax.legend(loc='upper right', fontsize=12)
#ax.set_xscale('log')

plt.tight_layout()
#plt.savefig("lc.png")
plt.show()
