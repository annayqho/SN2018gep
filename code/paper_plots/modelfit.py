""" Plot the physical evolution compared to David's model """

import numpy as np
from math import floor, log10
import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/code")
from scipy.optimize import curve_fit
from astropy.io import ascii
from astropy.table import Table
from load_lum import load_lc
from load_radius import load_radius
from load_temp import load_temp
from load_model import load


def round_sig(x, sig=2):
    if x < 0:
        return -round(-x, sig-int(floor(log10(-x)))-1)
    return round(x, sig-int(floor(log10(x)))-1)


def ndec(num):
    dec = str(num).split('.')[-1]
    return len(dec)


def powlaw(x, b, m):
    """ Fitting function for self-absorbed part """
    y = 10**(m*np.log10(x)+b)
    return y


def fit_pow(x, y, y_err=[], m=-5/3):
    """ If you know the power-law index already (e.g. -5/3)
    find the best-fit y-intercept """
    if len(y_err) > 0:
        popt, pcov = curve_fit(
                lambda x, b: powlaw(x, b, m), x, y,
                sigma = y_err, absolute_sigma = True, p0=[45])
    else:       
        popt, pcov = curve_fit(
                lambda x, b: powlaw(x, b, m), x, y, p0=[45])
    return popt, pcov


def print_table():
    # Load the bolometric light curve
    dt, lum, llum, ulum = load_lc()

    # Load the radius
    dt, rad, lrad, urad = load_radius()

    # Load the temperature
    dt, temp, ltemp, utemp = load_temp()

    # Table of measurements
    dtprint = np.array([round_sig(val,2) for val in dt])

    lprint = np.array([round_sig(val,2) for val in lum/3.839E43])
    ulprint = np.array(
            [np.round(val,ndec(lprint[ii])) \
            for ii,val in enumerate(ulum/3.839E43)])
    llprint = np.array(
            [np.round(val,ndec(lprint[ii])) \
            for ii,val in enumerate(llum/3.839E43)])

    rprint = np.array([round_sig(val,2) for val in rad/1.496E13])
    urprint = np.array(
            [np.round(val,ndec(rprint[ii])) \
            for ii,val in enumerate(urad/1.496E13)])
    lrprint = np.array(
            [np.round(val,ndec(rprint[ii])) \
            for ii,val in enumerate(lrad/1.496E13)])

    tprint = np.array([round_sig(val,2) for val in temp/1E3])
    utprint = np.array(
            [np.round(val,ndec(tprint[ii])) \
            for ii,val in enumerate(utemp/1E3)])
    ltprint = np.array(
            [np.round(val,ndec(tprint[ii])) \
            for ii,val in enumerate(ltemp/1E3)])

    outputf = open("physevol_tab.txt", "w")
    outputf.write("\\begin{table}[] \n")
    outputf.write("\centering \n")
    outputf.write(
            "\caption{Physical evolution of AT2018gep from blackbody fits.\
            Uncertainties represent the 16-to-84 percentile range from a\
            Monte Carlo simulation with 600 trials.} \n")
    outputf.write("\\begin{tabular}{lrrr} \n")
    outputf.write("\hline \n")
    outputf.write(
            "$\Delta t$ & $L (10^{10} L_\odot)$ & $R$ (AU) & $T$ (kK) \\\ \n")
    outputf.write("\hline")

    for ii,l in enumerate(lprint):
        linestr = "$%s$ & $%s^{%s}_{%s}$ & $%s^{%s}_{%s}$ & $%s^{%s}_{%s}$ \\\ \n" %(
                dtprint[ii], 
                l, "+%s" %ulprint[ii], "-%s" %llprint[ii],
                rprint[ii], "+%s" %urprint[ii], "-%s" %lrprint[ii],
                tprint[ii], "+%s" %utprint[ii], "-%s" %ltprint[ii])
        outputf.write(linestr)

    outputf.write("\hline \n")
    outputf.write("\end{tabular} \n")
    outputf.write("\label{tab:physevol} \n")
    outputf.write("\end{table} \n")

    outputf.close()


def plot():
    # Load the bolometric light curve
    dt, lum, llum, ulum = load_lc()

    # Load the radius
    dt, rad, lrad, urad = load_radius()

    # Load the temperature
    dt, temp, ltemp, utemp = load_temp()

    # Load David's model
    mdt, ml, mr, mt = load()     

    # Initialize the figure
    fig,axarr = plt.subplots(3,1, figsize=(6,8), sharex=True)

    # Luminosity panel
    axarr[0].errorbar(
            dt, lum, yerr=[llum,ulum], fmt='o', c='k', label="SN2018gep")
    axarr[0].plot(mdt, ml, c='k', lw=1, ls="--", label="Model")

    # Plot the SN
    ddir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/bol_lc"
    dat = Table.read(ddir + "/sn2010bh.dat", format='ascii.fast_no_header')
    sndt = dat['col1']
    snlum = dat['col2']
    xsn = sndt
    ysn = snlum*2
    axarr[0].plot(xsn, ysn, c='k', ls=':', lw=1, label='SN2010bh, x2')

    # Plot the total
    ysn_grid = np.interp(mdt, xsn, ysn) # put SN LC onto model grid
    tot = ysn_grid + ml
    axarr[0].plot(
            mdt, tot, c='k', ls='-', lw=1, label="Model + SN2010bhx2")

    # Legend
    axarr[0].legend(fontsize=14)

    # Radius panel
    axarr[1].errorbar(
            dt, rad, yerr=[lrad,urad], fmt='o', c='k', label="SN2018gep")
    axarr[1].plot(mdt, mr, c='k', lw=1, ls='--', label="Model")
    axarr[1].legend(fontsize=14, loc='lower right')

    # Temperature panel
    axarr[2].errorbar(
            dt, temp, yerr=[ltemp,utemp], fmt='o', c='k', label="SN2018gep")
    tsn = (tot/(4*np.pi*mr**2*5.67E-5))**0.25
    axarr[2].plot(mdt, mt, c='k', lw=1, ls='--', label="Model")
    axarr[2].plot(
            mdt, tsn, c='k', lw=1, ls='-', 
            label=r"$L_\mathrm{bol} / (4\pi R^2 \sigma)$")
    axarr[2].legend(fontsize=14, loc='upper right')

    # Formatting
    axarr[0].xaxis.label.set_visible(False)
    axarr[1].xaxis.label.set_visible(False)

    axarr[2].tick_params(axis='both', labelsize=16)
    axarr[0].tick_params(axis='y', labelsize=16)
    axarr[1].tick_params(axis='y', labelsize=16)

    for ii in np.arange(0,3):
        axarr[ii].set_yscale('log')

    axarr[1].set_ylim(1E14, 1E16)
    axarr[2].set_xlim(0, 18)
    #axarr[2].set_xscale('log')
    axarr[2].set_xlabel(r'Days since $t_0$', fontsize=16)
    axarr[0].set_ylabel(r'$L_\mathrm{bol}$ (erg/s)', fontsize=16)
    axarr[1].set_ylabel(r'$R_\mathrm{ph}$ (cm)', fontsize=16)
    axarr[2].set_ylabel(r'$T_\mathrm{eff}$ (K)', fontsize=16)

    plt.subplots_adjust(hspace=0)
    plt.tight_layout()
    #plt.show()
    plt.savefig("modelfit.png")


if __name__=="__main__":
    plot()
