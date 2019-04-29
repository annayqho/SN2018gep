""" Plot a nice Ic-BL spectrum. Let's use the NOT one. """

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15
from astropy.time import Time
import glob

DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/spec/ZTF18abukavn/"
f = DIR + "ZTF18abukavn_20180917_NOT_v1.ascii"

dat = np.loadtxt(f)
wl = dat[:,0]
flux = dat[:,1]

x = wl
y = flux 

fig,ax= plt.subplots(1,1,figsize=(8,3))
ax.plot(x, y/1E-15, c='k', drawstyle='steps-mid', lw=0.5)

ax.set_ylabel(
    r"Flux ($10^{-15}$ erg cm$^{-2}$ s$^{-1}$ \AA$^{-1}$)", 
    fontsize=16)
ax.set_xlim(3400,9100)
ax.set_ylim(0,3)
ax.set_xlabel(
        r"Observed Wavelength (\AA)", fontsize=16)
ax.xaxis.set_tick_params(labelsize=14)
ax.yaxis.set_tick_params(labelsize=14)

plt.tight_layout()
#plt.show()
plt.savefig("not_icbl.png")
