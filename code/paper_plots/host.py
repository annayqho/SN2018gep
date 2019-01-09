""" Plot a host image from SDSS """

import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.wcs
from astroquery.sdss import SDSS
from astropy import coordinates as coords

ra = 250.950926 
dec = 41.045380

ddir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/host"
gim = fits.open(ddir + "/frame-g-002335-3-0028.fits.bz2")

# Figure out pos from header
head = gim[0].header
wcs = astropy.wcs.WCS(head)
target_pix = wcs.wcs_world2pix([(np.array([ra,dec], np.float_))], 1)[0]
xpos = target_pix[0]
ypos = target_pix[1]

# Plot 600x600 cutout
im = gim[0].data[int(xpos-300):int(xpos+300),int(ypos-300):int(ypos+300)]

plt.imshow(im)

plt.show()
