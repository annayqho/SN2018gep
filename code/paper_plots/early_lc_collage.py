""" Collage of early light curves for all LLGRB-SNe """

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.io import ascii

fig, axarr = plt.subplots(2,2,sharex=True,sharey=True)

datadir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/lc"


def get_98bw_lc(band):
    dat = ascii.read(
        datadir + "/lc_980425.txt", format="fixed_width", delimiter="&")
    jd_raw = np.array(dat['JD'])
    delta = jd_raw[0] - 17/24 # first obs is 17 hrs after burst
    jd = jd_raw - delta
    lc = dat[band]
    bad = np.array(['nodata' in val for val in lc])
    mag = np.array([float(val.split("$\\pm$")[0]) for val in lc[~bad]])
    mag_err = np.array(
            [float(val.split("$\\pm$")[1]) for val in lc[~bad]])
    t = jd[~bad]
    return t, mag, mag_err


# Plot the GRB 980425 r-band early LC in each panel,
# color-coded by the Vc-Rc color
t, mag, mag_err = get_98bw_lc('Rc')
vt, vmag, vmag_err = get_98bw_lc('Vc')
for ax in axarr.reshape(-1):
    ax.plot(t, mag, c='lightgrey', alpha=0.7)
ax = axarr.reshape(-1)[0]
ax.errorbar(
    t, mag, yerr=mag_err, c='k', fmt='o', ms=3, alpha=1.0)
ax.text(
        0.9, 0.9, "1998bw, r-band", 
        fontsize=14, transform=ax.transAxes,
        horizontalalignment='right',
        verticalalignment='top',
        bbox=dict(
            boxstyle="round", fc='white', ec='k', 
            lw=0.5, alpha=1.0, pad=0.4))

ax.set_xlim(0,4)
ax.set_ylim(15,16.5)
ax.invert_yaxis()
plt.subplots_adjust(wspace=0, hspace=0)
plt.show()
