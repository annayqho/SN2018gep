""" Fit Nickel light curves / Arnett models """

import numpy as np
import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from astropy.table import Table


def az(y, z):
    return 2*z*np.exp(-2*z*y+z**2)


def bz(y, s, z):
    return 2*z*np.exp(-2*z*y+2*z*s+z**2)


def lph(x, y, s):
    eps_ni = 3.90E10
    eps_co = 6.78E9
    return mni * np.eps(-x**2) * \
            ((eps_ni-eps_co)*INTEGRAL(az(z)) * eps_co * INTEGRAL(bz(z)))


def load_lc():
    DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data"
    dat = np.loadtxt("%s/physevol.dat" %DATA_DIR, dtype=str)
    mjd = dat[:,0].astype(float)
    dt = mjd-58370.6473
    lum = dat[:,8].astype(float) * 3.839E33 # original units L_sun
    plt.plot(dt, lum)


if __name__=="__main__":
    # Load the bolometric light curve
    load_lc()
    plt.show()
