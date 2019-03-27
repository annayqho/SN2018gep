""" The various calculations for shock cooling,
using analytical formulas from Piro (2015) """

import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import numpy as np
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/code")
from load_lum import load_lc


def plot_14gqr():
    """ 
    Plot the bolometric light curve of iPTF14gqr
    """
    datadir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data"
    dat = np.loadtxt("%s/data/bol_lc/iPTF14gqr_bbLCRT.txt" %datadir)
    dt = dat[:,0]
    lum = dat[:,1]
    elum = dat[:,2]
    plt.errorbar(
            dt, lum, yerr=elum, fmt='s', 
            c='#e55c30', mec='#e55c30', mfc='#e55c30', label="iPTF14gqr")


def fit(t, kappa, Mc, E51, Re, Me):
    # Produce a bolometric light curve
    tp_day = 0.9*(kappa/0.34)**(0.5) * E51**(-0.25) * (Mc)**(0.17) * (Me/0.01)**0.57
    tp = tp_day * 86400
    ve = 2E9 * E51**0.5 * Mc**(-0.35) * (Me/0.01)**(-0.15)
    te = Re/ve
    Ee = 4E49 * E51 * Mc**(-0.7) * (Me/0.01)**0.7

    L = (te*Ee/tp**2) * np.exp(-(t*(t+2*te))/(2*tp**2))
    return t/86400, L


def magnetar():
    # set parameters
    B = 1E14
    Pms = 20

    tm = 4 * B14**(-2) * Pms**2
    Lm = (Em/tm)/(1+t/tm)**2


def plot_18gep():
    # Plot the data
    dt, lum, llum, ulum = load_lc()
    plt.errorbar(
            dt, lum, yerr=[llum,ulum], 
            fmt='o', c='#140b34', mfc='#140b34', mec='#140b34',
            lw=0.5, label="SN2018gep")


if __name__=="__main__":
    fig = plt.figure(figsize=(8,6))
    t = np.linspace(0,10,1000) * 86400

    plot_14gqr()
    # order: t, kappa, Mc, E51, Re, Me
    x,y = fit(t, 0.2, 0.23, 1.38E50/1E51, 3E13, 8.8E-3)
    plt.plot(
            x, y, c='#e55c30', 
            label=r"$M_\mathrm{sh}=8.8 \times 10^{-3}\,M_\odot$,$R_\mathrm{sh}=3 \times 10^{13}$\,cm,$E_K=2 \times 10^{50}\,$erg")

    plot_18gep()
    x, y = fit(t, 0.2, 0.23, 10, 3E13, 8.8E-3)
    plt.plot(
            x, y, c='#140b34', ls='--',
            label=r"$M_\mathrm{sh}=8.8 \times 10^{-3}\,M_\odot$,$R_\mathrm{sh}=3 \times 10^{13}$\,cm,$E_K=10^{52}\,$erg")

    x, y = fit(t, 0.2, 0.23, 10, 7E12, 0.5)
    plt.plot(
            x, y, c='#140b34', ls='-',
            label=r"$M_\mathrm{sh}=0.5\,M_\odot$,$R_\mathrm{sh}=7 \times 10^{12}$\,cm,$E_K=10^{52}\,$erg")

    plt.xlim(-0.2,8)
    plt.ylim(5E41,1E46)
    plt.gca().tick_params(labelsize=14)
    plt.xlabel("$\Delta t$ (days)", fontsize=16)
    plt.ylabel("$L_\mathrm{bol}$ (erg/s)", fontsize=16)
    plt.yscale('log')
    plt.legend(loc='upper right', fontsize=12)
    plt.show()
