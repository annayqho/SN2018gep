""" The various calculations for shock cooling,
using analytical formulas from Piro (2015) """

import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import numpy as np
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/code")
from astropy.time import Time
from load_lum import load_lc


def plot_14gqr():
    """ 
    Plot the bolometric light curve of iPTF14gqr
    """
    datadir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data"
    dat = np.loadtxt("%s/bol_lc/iPTF14gqr_bbLCRT.txt" %datadir)
    dt = dat[:,0]
    lum = dat[:,1]
    elum = dat[:,2]
    plt.errorbar(
            dt, lum, yerr=elum, fmt='s', 
            c='#e55c30', mec='#e55c30', mfc='#e55c30', label="iPTF14gqr")


def fit(t, kappa, Mc, E51, Re, Me):
    """ Using a set of parameters, produce a light curve according to
    the model from Tony Piro's 2015 paper
    
    t = time in days
    kappa: opacity in cgs
    Mc: core mass in solar masses
    Re: envelope radius in cm
    Me: envelope mass in solar masses
    """
    # The time to peak luminosity, from Nakar & Piro 2014
    tp_day = 0.9 * (kappa/0.34)**(0.5) * \
            E51**(-0.25) * (Mc)**(0.17) * (Me/0.01)**0.57
    tp = tp_day * 86400

    # The velocity in the extended material (cm/s)
    ve = 2E9 * E51**0.5 * Mc**(-0.35) * (Me/0.01)**(-0.15)

    # the expansion timescale in seconds
    te = Re/ve

    # The amount of energy passed into the extended material (erg)
    Ee = 4E49 * E51 * Mc**(-0.7) * (Me/0.01)**0.7

    # The light curve
    t_sec = t*86400
    L = (te*Ee/tp**2) * np.exp(-(t_sec*(t_sec+2*te))/(2*tp**2))

    # Return time in days
    return t, L


    # set parameters
    B = 1E14
    Pms = 20

    tm = 4 * B14**(-2) * Pms**2
    Lm = (Em/tm)/(1+t/tm)**2


def plot_18gep():
    """ Plot the bolometric LC of 2018gep """
    dt, lum, llum, ulum = load_lc()
    plt.errorbar(
            dt, lum, yerr=[llum,ulum], 
            fmt='o', c='#140b34', mfc='#140b34', mec='#140b34',
            lw=0.5, label="SN2018gep")


def plot_16asu():
    """ Plot the bolometric LC of 2018gep """
    datadir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data"
    t0 = Time("2016-05-10").jd + 0.53
    et0 = 0.17
    dat = np.loadtxt(datadir + "/bol_lc/iPTF16asu_bolometric_luminosity.txt",
        delimiter=" ")
    dt = dat[:,0]-t0
    lum = dat[:,1]
    elum = dat[:,2]
    c = '#84206b'
    plt.errorbar(
            dt, lum, xerr=0.17, yerr=elum,
            fmt='v', mec=c, mfc='white', c=c,
            label="iPTF16asu", zorder=0, lw=0.5)


if __name__=="__main__":
    fig = plt.figure(figsize=(8,6))
    t = np.linspace(0,10,1000) 

    # Parameters for 14gqr
    kappa = 0.2
    Mc = 0.23 # solar mass
    E = 1.38E50
    Re = 3E13
    Me = 8.8E-3

    plot_14gqr()
    x,y = fit(t, kappa, Mc, E/1E51, Re, Me)
    plt.plot(
            x, y, c='#e55c30', 
            label=r"$M_\mathrm{sh}=8.8 \times 10^{-3}\,M_\odot$,$R_\mathrm{sh}=3 \times 10^{13}$\,cm,$E_K=1.38 \times 10^{50}\,$erg")

    # Parameters for 18gep
    kappa = 0.2
    Mc = 0.23 # solar mass
    E = 1E52
    Re = 3E13
    Me = 8.8E-3
    plot_18gep()
    x,y = fit(t, kappa, Mc, E/1E51, Re, Me)
    plt.plot(
            x, y, c='#140b34', ls='--',
            label=r"$M_\mathrm{sh}=8.8 \times 10^{-3}\,M_\odot$,$R_\mathrm{sh}=3 \times 10^{13}$\,cm,$E_K=10^{52}\,$erg")

    # Make the shell way more massive
    kappa = 0.2
    Mc = 0.23 # solar mass
    E = 1E52
    Re = 7E12
    Me = 0.5
    x, y = fit(t, kappa, Mc, E/1E51, Re, Me)
    plt.plot(
            x, y, c='#140b34', ls='-',
            label=r"$M_\mathrm{sh}=0.5\,M_\odot$,$R_\mathrm{sh}=7 \times 10^{12}$\,cm,$E_K=10^{52}\,$erg")

    # What about iPTF16asu?
    # leave this out for now...
    kappa = 0.2
    Mc = 0.23
    E = 3.8E51
    Re = 1.7E12
    Me = 0.45
    x, y = fit(t, kappa, Mc, E/1E51, Re, Me)

    plt.xlim(-0.2,8)
    plt.ylim(5E41,1E46)
    plt.gca().tick_params(labelsize=14)
    plt.xlabel("$\Delta t$ (days)", fontsize=16)
    plt.ylabel("$L_\mathrm{bol}$ (erg/s)", fontsize=16)
    plt.yscale('log')
    plt.legend(loc='upper right', fontsize=12)
    plt.show()
    #plt.savefig("piro_models.png")
