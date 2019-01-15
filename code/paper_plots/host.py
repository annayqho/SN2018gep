""" Plot a host image from SDSS """

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.wcs
from astroquery.sdss import SDSS
from astropy import coordinates as coords

ra = 250.950926 
dec = 41.045380

def combine_images(band1, band2):
    ddir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/host"
    im1 = fits.open(ddir + "/frame-%s-002335-3-0028.fits" %band1)
    im2 = fits.open(ddir + "/frame-%s-002335-3-0028.fits" %band2)
    data1 = im1[0].data
    data2 = im2[0].data
    comb = (data1 + data2)/2
    hdu = fits.PrimaryHDU(comb)
    hdu.writeto('%s%sim.fits' %(band1,band2), overwrite=True)

combine_images('g', 'r')
combine_images('i', 'z')

# Figure out pos from header
# head = gim[0].header
# Figure out pos from header
# head = gim[0].header
# wcs = astropy.wcs.WCS(head)
# target_pix = wcs.wcs_world2pix([(np.array([ra,dec], np.float_))], 1)[0]
# xpos = target_pix[0]
# ypos = target_pix[1]
# 
# # Plot 600x600 cutout
# #im = gim[0].data[int(xpos-300):int(xpos+300),int(ypos-300):int(ypos+300)]
# im = gim[0].data
# 
# plt.imshow(im)
# 
# plt.show()
