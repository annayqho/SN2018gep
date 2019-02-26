""" Collage of early light curves for all LLGRB-SNe """

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.io import ascii

fig, axarr = plt.subplots(2,2,sharex=True,sharey=True)

datadir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/lc"

# Plot the GRB 980425 early LC
dat = ascii.read(
    datadir + "/lc_980425.txt", format="fixed_width", delimiter="&")
jd_raw = np.array(dat['JD'])
delta = jd_raw[0] - 17/24 # first obs is 17 hrs after burst
jd = jd_raw - delta
bands = ['B','Vc','Rc']
shapes = ['o', 'o', 's']
colors = ['black', 'white', 'black']
ls = ['-', '--', ':']
#colors = ['blue', 'green', 'red', 'lightsalmon']
for ii,band in enumerate(bands):
    lc = dat[band]
    bad = np.array(['nodata' in val for val in lc])
    mag = np.array([float(val.split("$\\pm$")[0]) for val in lc[~bad]])
    mag_err = np.array(
            [float(val.split("$\\pm$")[1]) for val in lc[~bad]])
    t = jd[~bad]
    for ax in axarr.reshape(-1):
        ax.plot(t,mag,c='k', linestyle=ls[ii])
        ax.errorbar(
                t,mag,yerr=mag_err,ecolor='k', fmt=shapes[ii],
                mfc=colors[ii],mec='k',label=band,ms=5)
        ax.legend()

ax.set_xlim(0,4)
ax.set_ylim(15,16.5)
ax.invert_yaxis()
plt.subplots_adjust(wspace=0, hspace=0)
plt.show()
