""" Plot a host image from SDSS

This is Fig 1 in the paper.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.io import fits
import astropy.wcs
from astropy import coordinates as coords
from astropy.visualization import make_lupton_rgb

ra = 250.950926 
dec = 41.045380

ddir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/host"

u = fits.open(
        ddir + "/G019.250.801+41.089.U.8087_8342_10281_10536.fits")[0].data
r = fits.open(
        ddir + "/G019.250.801+41.089.R.8087_8342_10281_10536.fits")[0].data
g = fits.open(
        ddir + "/G019.250.801+41.089.G.8087_8342_10281_10536.fits")[0].data
z = fits.open(
        ddir + "/G019.250.801+41.089.Z.8087_8342_10281_10536.fits")[0].data

rgb = make_lupton_rgb(z/8, r/4, u, minimum=[2,2,2], stretch=20, Q=10)

fig,ax = plt.subplots()

ax.text(0.05, 0.9, "CFHTLS/$urz$", fontsize=24, transform=ax.transAxes,
        horizontalalignment='left', color='white')

ax.imshow(rgb[35:240,35:240,:], origin='lower')

# mark the transient location
imsize = 240-35+1
ax.plot([imsize/2, imsize/2], [imsize/2, imsize/2-10], c='white')
ax.plot([imsize/2, imsize/2+10], [imsize/2, imsize/2], c='white')

# Mark compass
ax.plot((imsize-10,imsize-10), (imsize-10,imsize-20), color='white', lw=2)
ax.text(
        imsize-10, imsize-23, "S", color='white', fontsize=16,
        horizontalalignment='center', verticalalignment='top')
ax.plot((imsize-10,imsize-20), (imsize-10,imsize-10), color='white', lw=2)
ax.text(
        imsize-23, imsize-10, "E", color='white', fontsize=16,
        horizontalalignment='right', verticalalignment='center')
ax.axis('off')

# Mark image scale
# I think it's 0.186 arcsec
x = 10 
y = 20
x2 = x + 5/0.186
ax.plot((x,x2), (y,y), c='white', lw=2)
ax.text((x2+x)/2, y*1.1, "5''", color='white', fontsize=16,
        verticalalignment='bottom', horizontalalignment='center')
ax.text((x2+x)/2, y/1.1, "(3.3 kpc)", color='white', fontsize=16,
        verticalalignment='top', horizontalalignment='center')

plt.savefig("host.eps", dpi=500, bbox_inches='tight')

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
