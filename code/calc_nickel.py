""" See whether the LC is consistent with being powered by Ni decay at
late times """

import numpy as np
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/code")
from load_lum import load_lc


def qgamma(t, mni):
    """
    t: time in days
    mni: nickel mass in solar masses
    """
    return (mni)*(6.45*np.exp(-t/8.76) + 1.38*np.exp(-t/111.4))*1E43


def qpos(t, mni):
    """
    t: time in days
    mni: nickel mass in solar masses
    """
    return 4.64 * (mni) * (-np.exp(-t/8.76) + np.exp(-t/111.4)) * 1E41


def qdep(t, t0, mni):
    """ The total energy deposition rate """
    qdep = qgamma(t, mni) * (1-np.exp(-(t0/t)**2)) + qpos(t, mni)


if __name__=="__main__":
    # load the bolometric light curve
    dt, lum, llum, ulum = load_lc()

    # calculate the ratio of the bolometric luminosity
    # to the integrated bol LC

    # since you integrate from the first point,
    # you can't integrate the first point from itself
    lum_ratio = np.zeros(len(dt)-1)

    # for each point after the first one, evaluate...
    for ii in np.arange(1, len(dt)-1):
        # the integral up until and including that point
        integral = np.trapz((lum*dt)[0:ii+1], dt[0:ii+1])
        lum_ratio[ii-1] = lum[ii] / integral

    # now, calculate the energy deposition integral
    qdep(t, t0, mni)
