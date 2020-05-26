""" Plot the UVOT grism spectrum """

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')
from astropy.io import fits as pyfits
from astropy.cosmology import Planck15
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/code")
from load_lc import get_uv_lc


# Dictionary mapping the filter to the wavelength
bands = {}
bands['V'] = 5468
bands['B'] = 4392
bands['U'] = 3465
bands['UVW1'] = 2600
bands['UVM2'] = 2246
bands['UVW2'] = 1928


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

# flux_grid is in units erg/cm2/s/A.
# now we want to scale to the UV photometry

# get the U-filter measurement at 6.4 days

# in mJy
dt, filt, lnu, elnu = get_uv_lc()
# transform to cgs
lnu_cgs = lnu * 1E-3 * 1E-23  # erg/cm2/s/Hz
elnu_cgs = elnu * 1E-3 * 1E-23  # erg/cm2/s/Hz

# effective frequency of U-band filter
nu = 3E18 / bands['U'] # Hz
u_dt = dt[filt == 'U']
u_lnu = lnu_cgs[filt == 'U']

# interpolate to epoch of UV grism spectrum
u_lnu_val = np.interp(6.4, u_dt, u_lnu)

# Convert from per Hz to per AA
u_lwl_val = u_lnu_val * (bands['U']*1E8)**2 / 3E10

# Then, integrate the spectrum from 3000 to 3900 AA
choose = np.logical_and(wl_grid > 3000, wl_grid < 3900)
spec_int = np.trapz(flux_grid[choose], wl_grid[choose]) # 3.2E-12

# spec_int is slightly higher than u_lum
# Scaling factor
scale = u_lwl_val / spec_int

fig,ax = plt.subplots(1,1,figsize=(6,4))

flux_scaled = flux_grid * scale
ax.plot(wl_grid, flux_scaled, drawstyle='steps-mid', lw=0.5, c='k')
ax.fill_between(
        wl_grid,
        (flux_grid-eflux_grid)*scale, 
        (flux_grid+eflux_grid)*scale, 
        color='lightgrey', alpha=0.3)

# The Time of the UVOT grism spectrum: 2018-09-15T13:29:46.795880
ax.text(
        0.9, 0.9, r"$\Delta t = 6.4$ days", 
        transform=ax.transAxes, fontsize=14,
        horizontalalignment='right',
        verticalalignment='top')

ax.set_ylim(9E-18, 2E-16)
ax.set_ylabel(r"Flux (erg cm$^{-2}$ s$^{-1}$ \AA$^{-1}$)", fontsize=16)
ax.set_yscale('log')
ax.set_xlabel(r"Rest-Frame Wavelength (\AA)", fontsize=16)
ax.tick_params(axis='both', labelsize=14)

plt.tight_layout()

#plt.show()
plt.savefig("uv_grism_spec.eps", dpi=300)
