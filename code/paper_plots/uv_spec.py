""" Plot the UVOT grism spectrum """

import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')
from astropy.io import fits as pyfits

ddir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/uv"
dat = ddir + "/sumpha2.fit"
spec = pyfits.open(dat)[1].data
wl = spec['wave']
flux = spec['flux']
err = spec['sqrtvariance']
choose = wl <= 4100

fig,ax = plt.subplots(1,1,figsize=(6,4))

ax.plot(wl[choose], flux[choose], drawstyle='steps-mid', lw=0.5, c='k')
ax.fill_between(
        wl[choose], flux[choose]-err[choose], 
        flux[choose]+err[choose],
        color='grey', alpha=0.3)

# The Time of the UVOT grism spectrum: 2018-09-15T13:29:46.795880
ax.text(
        0.9, 0.9, r"$\Delta t = 6.4$ days", 
        transform=ax.transAxes, fontsize=12,
        horizontalalignment='right',
        verticalalignment='top')

ax.set_ylim(3E-16, 4E-14)
ax.set_ylabel(r"Flux (erg cm$^{-2}$ s$^{-1}$ \AA$^{-1}$", fontsize=14)
ax.set_yscale('log')
ax.set_xlabel(r"Rest-Frame Wavelength (\AA)", fontsize=14)
ax.tick_params(axis='both', labelsize=12)

plt.tight_layout()

#plt.show()
plt.savefig("uv_grism_spec.png")
