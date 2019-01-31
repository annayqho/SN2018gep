""" Plot the UVOT grism spectrum """

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')
from astropy.io import fits as pyfits

ddir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/uv"
dat = ddir + "/sumpha2.fit"
spec = pyfits.open(dat)[1].data
wl = spec['wave'] # every Angstrom
flux = spec['flux']
err = spec['sqrtvariance']
err[err==0] = 0.0000001
ivar = 1/err**2
choose = wl <= 4100

# Make a grid of wavelengths, from 2300 to 4000
binsize = 1
wl_grid = np.arange(2305, 4000, binsize)
flux_grid = np.zeros(len(wl_grid))
eflux_grid = np.zeros(len(wl_grid))

for ii,wl_val in enumerate(wl_grid):
    choose = np.abs(wl-wl_val) <= 5
    if sum(err[choose] == 0)>0:
        print(wl_val)
    w = ivar[choose]
    mean,wsum = np.average(
            flux[choose], weights=w, returned=True)
    flux_grid[ii] = mean
    emean = np.sqrt(1/np.sum(w))
    eflux_grid[ii] = emean

fig,ax = plt.subplots(1,1,figsize=(6,4))

ax.plot(wl_grid, flux_grid, drawstyle='steps-mid', lw=0.5, c='k')
ax.fill_between(
        wl_grid,
        flux_grid-eflux_grid, 
        flux_grid+eflux_grid, 
        color='grey', alpha=0.3)

# The Time of the UVOT grism spectrum: 2018-09-15T13:29:46.795880
ax.text(
        0.9, 0.9, r"$\Delta t = 6.4$ days", 
        transform=ax.transAxes, fontsize=14,
        horizontalalignment='right',
        verticalalignment='top')

ax.set_ylim(3E-16, 1E-14)
ax.set_ylabel(r"Flux (erg cm$^{-2}$ s$^{-1}$ \AA$^{-1}$)", fontsize=16)
ax.set_yscale('log')
ax.set_xlabel(r"Rest-Frame Wavelength (\AA)", fontsize=16)
ax.tick_params(axis='both', labelsize=14)

plt.tight_layout()

#plt.show()
plt.savefig("uv_grism_spec.png")
