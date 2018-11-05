""" Show the X-ray limits compared to the iPTF16asu X-ray light curves """

import glob
import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15

def plot_lc(f, name=None):
    dt = []
    lum = []
    with open(f, "r") as inputf:
        for line in inputf.readlines():
            if len(line) > 40:
                dt_s = float(line.split('\t')[0])
                dt.append(dt_s/86400)
                flux = float(line.split('\t')[3])
                lum.append(flux * 4 * np.pi * d**2)
    if name:
        plt.plot(dt, lum, c='grey', alpha=1.0, lw=2, zorder=1)
    else:
        # just thin grey
        plt.plot(dt, lum, c='grey', alpha=0.1, lw=1)
    return dt, lum


# ZTF18abukavn
d = Planck15.luminosity_distance(z=0.033).cgs.value
ratio = 3.27E-11 # count-to-flux rate, erg/cm2/ct

xrtlc = Table.read(
    "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/from_brad.dat",
    format='ascii')
dt_sec = xrtlc['col1']
ct = xrtlc['col4']
lum = ct * ratio * 4 * np.pi * d**2 
plt.scatter(dt_sec/86400, lum, marker='v', c='k', zorder=5)

# All the GRBs
data_dir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/xrtlc"
# flist = glob.glob(data_dir + "/*_xrt_bin.txt")
# for f in flist:
#     plot_lc(f)

# Individual LLGRBs
f = data_dir + "/060218_xrt_bin.txt" 
dt, lum = plot_lc(f, name="060218")
plt.text(
    0.06, 1.26E46, "060218", 
    horizontalalignment='center',
    verticalalignment='bottom', fontsize=14)

f = data_dir + "/100316d_xrt_bin.txt" 
dt, lum = plot_lc(f, name="100316D")
plt.text(
    dt[-1], lum[-1], '100316D', 
    horizontalalignment='right',
    verticalalignment='top', fontsize=14)

f = data_dir + "/030329_xray.dat"
# this one has a different file
dat = np.loadtxt(f)
dt = dat[:,0] 
flux = dat[:,1] * 1E-12
lum = flux * 4 * np.pi * d**2
plt.plot(dt, lum, c='grey', lw=2, zorder=1)
plt.text(
    dt[0], lum[0], '030329',
    horizontalalignment='center',
    verticalalignment='bottom', fontsize=14)

# need 2009bb
# need 1998bw

# Formatting
plt.xlim(1E-3,1E2)
plt.xlabel(r"$\Delta t$ [days]", fontsize=14)
plt.ylabel(
r"X-ray luminosity [0.3-10 keV, erg\,s${}^{-1}$]", fontsize=14)
plt.yscale('log')
plt.xscale('log')
plt.ylim(1E40, 1E47)
plt.tick_params(axis='both', labelsize=14)

plt.tight_layout()
#plt.show()
plt.savefig("xray_lc.png")
