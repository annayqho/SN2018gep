""" Integrate blackbody across a ZTF filter """

import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling.blackbody import blackbody_lambda, blackbody_nu
from astropy import units as u
from astropy.cosmology import Planck15


def get_flux(mag):
    return 1E-23 * 3631 * 10**(-mag/2.5)


def make_bb(T, R):
    """ 
    Given a certain temperature, 
    simulate measurements
    in ZTF g- and r-band 
    """
    temp = teff * u.K

    d = Planck15.luminosity_distance(z=0.03154).cgs.value
    datadir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data"
    gfilt = np.loadtxt(datadir + "/Palomar_ZTF.g.dat")
    rfilt = np.loadtxt(datadir + "/Palomar_ZTF.r.dat")

    gwl = gfilt[:,0] * u.AA
    gtrans = gfilt[:,1]
    g_flam = blackbody_lambda(gwl, temp)
    g_f = g_flam*gtrans
    g_l = np.pi * g_f * 4 * np.pi * d**2 

    rwl = rfilt[:,0] * u.AA
    rtrans = rfilt[:,1]
    r_flam = blackbody_lambda(rwl, temp)
    r_f = r_flam*rtrans

    # int_0^infty Blam(T)dlam = sig*T**4 / pi

    # Spectral flux
    h = 6.63E-27
    c = 3E10
    k = 1.38E-16
    blam = (2 * h * c**2 / wl**5) / (np.exp(h*c/(wl*k*T))-1)
    spflux = np.pi * blam


if __name__=="__main__":
    rmag = 19.745
    rmeas = get_flux(rmag)
    gmag = 19.395
    gmeas = get_flux(gmag)
    #make_bb(5000)
