""" Plot the UVOT grism spectrum """

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')
from astropy.io import fits as pyfits
from astropy.cosmology import Planck15
from uv_lc import get_uv_lc

d = Planck15.luminosity_distance(z=0.3154).cgs.value

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

# Scale to the UV photometry.
# The U-band filter is roughly 3000-3900
# http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=Swift/UVOT.U
# First, interpolate the UV LC for U-filter to 6.4 days
dt, filt, nu, lnu, elnu = get_uv_lc()
u_dt = dt[filt == 'U']
u_lnu = lnu[filt == 'U']
u_lnu_val = np.interp(6.4, u_dt, u_lnu)
u_lum = u_lnu_val * nu[filt == 'U'][0] / (4 * np.pi * d**2)
# Then, integrate the spectrum from 3000 to 3900 AA
choose = np.logical_and(wl_grid > 3000, wl_grid < 3900)
spec_int = np.trapz(flux_grid[choose], wl_grid[choose]) # 3.2E-12
# spec_int is slightly higher than u_lum
# Scaling factor
scale = u_lum / spec_int

fig,ax = plt.subplots(1,1,figsize=(6,4))

ax.plot(wl_grid, flux_grid*scale, drawstyle='steps-mid', lw=0.5, c='k')
ax.fill_between(
        wl_grid,
        (flux_grid-eflux_grid)*scale, 
        (flux_grid+eflux_grid)*scale, 
        color='grey', alpha=0.3)

# The Time of the UVOT grism spectrum: 2018-09-15T13:29:46.795880
ax.text(
        0.9, 0.9, r"$\Delta t = 6.4$ days", 
        transform=ax.transAxes, fontsize=14,
        horizontalalignment='right',
        verticalalignment='top')

ax.set_ylim(1E-17, 3E-16)
ax.set_ylabel(r"Flux (erg cm$^{-2}$ s$^{-1}$ \AA$^{-1}$)", fontsize=16)
ax.set_yscale('log')
ax.set_xlabel(r"Rest-Frame Wavelength (\AA)", fontsize=16)
ax.tick_params(axis='both', labelsize=14)

plt.tight_layout()

#plt.show()
plt.savefig("uv_grism_spec.png")
