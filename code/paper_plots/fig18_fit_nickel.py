""" Fit Nickel light curves / Arnett models """

import numpy as np
import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/code")
from load_lum import load_lc
from scipy.integrate import quad
from scipy.optimize import curve_fit
from astropy.table import Table

# CONSTANTS
eps_ni = 3.90E10 # cgs
eps_co = 6.78E9 # cgs
t_ni = 8.80 * 86400 # days to sec
t_co = 111.3 * 86400 # days

def az(z, y):
    return 2*z*np.exp(-2*z*y+z**2)


def bz(z, s, y):
    return 2*z*np.exp(-2*z*y+2*z*s+z**2)


def lph(dt, mni, tdiff):
    # dt comes in days, need to convert to seconds
    t = dt * 86400

    # mni comes in solar masses, need to convert to grams
    mni_g = mni * 1.989E33

    # tdiff comes in days, need to convert to seconds
    t_m = tdiff * 86400

    x = t/t_m
    y = t_m/(2*t_ni)
    s = t_m*(t_co-t_ni) / (2*t_co*t_ni)

    # integrate A(z)dz from 0 to x
    a = np.array([quad(az, 0, xval, args=(y,))[0] for xval in x])
    aerr = np.array([quad(az, 0, xval, args=(y,))[1] for xval in x])

    # integrate B(z)dz from 0 to x
    b = np.array([quad(bz, 0, xval, args=(s, y))[0] for xval in x])
    berr = np.array([quad(bz, 0, xval, args=(s, y))[1] for xval in x])

    lum = mni_g * np.exp(-x**2) * ((eps_ni-eps_co)*a + eps_co*b)
    return lum


if __name__=="__main__":
    # Load the bolometric light curve
    dt, lum, llum, ulum = load_lc()

    # Plot a sample Ni light curve, say M_Ni = 0.1 solar masses
    mni0 = 0.3 # solar masses
    tdiff0 = 10 # days

    #kappa = 0.1 # cm2/g

    # Plot
    fig,ax = plt.subplots(1,1,figsize=(7,5))
    ax.errorbar(
            dt, lum, yerr=[llum, ulum], fmt='o', 
            mec='k', mfc='k', ms=5, c='k')

    # Various models
    dt_plot = np.linspace(0.1, 40, 1000)

    # same as in Whitesides for iPTF16asu
    ni_lum = lph(dt_plot, 0.55, 1.5)
    ax.plot(dt_plot, ni_lum, ls='--', c='k')
    ax.text(
            30, 1E43, 
            r"$M_\mathrm{Ni}=0.55 M_\odot, \tau_m=1.5\,\mathrm{d}$",
            horizontalalignment='center')

    # Best-fit to 10+ days
    choose = dt >= 10
    popt, pcov = curve_fit(
            lph, dt[choose], lum[choose], p0=[mni0, tdiff0])#, bounds=(0, [1.0, 100]))
    mni = popt[0]
    tdiff = popt[1]
    ni_lum = lph(dt_plot, mni, tdiff)
    ax.plot(dt_plot, ni_lum, ls='--', c='orange')
    ax.text(
            33, 4E42, 
            r"$M_\mathrm{Ni}= 0.27 M_\odot, \tau_m=5.5\,\mathrm{d}$",
            horizontalalignment='center')

    # Best-fit to the last two days
    choose = dt > 28
    popt, pcov = curve_fit(
            lph, dt[choose], lum[choose], p0=[mni0, tdiff0])
    mni = popt[0]
    tdiff = popt[1]
    dt_plot = np.linspace(0.1, 40, 1000)
    ni_lum = lph(dt_plot, mni, tdiff)
    ax.plot(dt_plot, ni_lum, ls='--', c='purple')
    ax.text(
            38, 2E42, 
            r"$M_\mathrm{Ni}= 0.13 M_\odot, \tau_m=15\,\mathrm{d}$",
            horizontalalignment='center')
    
    ax.set_ylabel("Bolometric Luminosity (erg/s)", fontsize=16)
    ax.set_xlabel(
        r"Days since $t_0=$JD 2458370.6473 (UT 2018 Sept 09.15)", fontsize=16)
    ax.yaxis.set_tick_params(labelsize=14)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.set_yscale('log')
    ax.set_ylim(1E42, 4E44)
    ax.set_xlim(-3, 45)
    plt.show()
    #plt.savefig("nifit.png")
