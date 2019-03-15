""" Plot the physical evolution as a power law """

import numpy as np
from math import floor, log10
import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/code")
from scipy.optimize import curve_fit
from astropy.io import ascii
from load_lum import load_lc
from load_radius import load_radius
from load_temp import load_temp


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
    rprint = np.array([int(r) if r>=10 for r in rprint])
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


    # Initialize the figure
    fig,axarr = plt.subplots(3,1, figsize=(6,8), sharex=True)

    # Luminosity panel
    axarr[0].errorbar(dt, lum, yerr=[llum,ulum], fmt='o', c='k')
    # Fit power law to all points after 1 day
    m = -5/3
    mstr = '-5/3'
    choose = dt > 1
    b,berr = fit_pow(
            dt[choose], lum[choose], np.min((llum, ulum), axis=0)[choose], m=m)
    xfit = np.linspace(1, max(dt))
    yfit = powlaw(xfit, b, m)
    axarr[0].plot(xfit, yfit, ls='--', c='#f98e09')
    axarr[0].text(25, 5E42, '$t^{%s}$' %mstr,
            horizontalalignment='left', verticalalignment='center', 
            fontsize=14)
    axarr[0].axvline(x=3.22, ls='-', c='lightgrey', lw=3, zorder=0)
    axarr[0].text(
            3, 1E43, "$t_\mathrm{rise}<3.2\,$d", fontsize=14,
            horizontalalignment='right')

    # Fit another power law
    m = -2
    mstr = '-2'
    choose = dt > 1
    b,berr = fit_pow(
            dt[choose], lum[choose], np.min((llum, ulum), axis=0)[choose], m=m)
    xfit = np.linspace(1, max(dt))
    yfit = powlaw(xfit, b, m)
    axarr[0].plot(xfit, yfit, ls='--', c='#57106e')
    axarr[0].text(2, 5E44, '$t^{%s}$' %mstr,
            horizontalalignment='left', verticalalignment='center', fontsize=14)

    # Fit a power law to points before 3.5 days
    choose = dt <= 3.5
    m,b = np.polyfit(np.log10(dt[choose]), np.log10(lum[choose]), deg=1)
    xfit = np.linspace(0.3, 5)
    yfit = 10**(m*np.log10(xfit)+b)
    mstr = str(np.round(m, 1))
    axarr[0].plot(xfit, yfit, ls='--', c='grey')
    axarr[0].text(1, 1E44, '$t^{%s}$' %mstr,
            horizontalalignment='left', verticalalignment='center', fontsize=14)

    # Formatting

    # Radius panel
    axarr[1].errorbar(dt, rad, yerr=[lrad,urad], fmt='o', c='k')

    # Plot lines of constant velocity
    xvals = np.linspace(1E-3, 1E2, 1000)

    # v = 0.1c
    yvals = 0.1 * (3E10) * xvals * 86400
    axarr[1].plot(xvals, yvals, ls='--', c='grey')
    axarr[1].text(1, 2E14, 'v=0.1c', fontsize=14)

    # v = 0.26c
    #yvals = 0.26 * (3E10) * xvals * 86400
    #axarr[1].plot(xvals, yvals, ls='--', c='grey')
    #axarr[1].text(0.4, 7E14, 'v=0.26c', fontsize=14)

    # Temperature panel
    choose = np.logical_and(dt>1, dt<19)
    axarr[2].errorbar(dt, temp, yerr=[ltemp,utemp], fmt='o', c='k')
    m = -0.92
    b,berr = fit_pow(
            dt[choose], temp[choose], 
            np.min((ltemp, utemp), axis=0)[choose], m=m)
    xfit = np.linspace(1, max(dt))
    yfit = 10**(m*np.log10(xfit)+b)
    mstr = str(np.round(m, 1))
    axarr[2].plot(xfit, yfit, ls='--', c='#f98e09')
    axarr[2].text(10, 1E4, '$t^{%s}$' %mstr,
            horizontalalignment='left', verticalalignment='center', 
            fontsize=14)
    choose = dt <= 3.5
    m,b = np.polyfit(np.log10(dt[choose]), np.log10(temp[choose]), deg=1)
    xfit = np.linspace(0.3, 5)
    yfit = 10**(m*np.log10(xfit)+b)
    mstr = str(np.round(m, 1))
    axarr[2].plot(xfit, yfit, ls='--', c='grey')
    axarr[2].text(1, 2E4, '$t^{%s}$' %mstr,
            horizontalalignment='left', verticalalignment='center', 
            fontsize=14)

    # Showing where 15,000 K is
    axarr[2].axhline(y=15000)

    # Formatting
    axarr[0].xaxis.label.set_visible(False)
    axarr[1].xaxis.label.set_visible(False)

    axarr[2].tick_params(axis='both', labelsize=16)
    axarr[0].tick_params(axis='y', labelsize=16)
    axarr[1].tick_params(axis='y', labelsize=16)

    for ii in np.arange(0,3):
        axarr[ii].set_yscale('log')

    axarr[1].set_ylim(1E14, 1E16)
    axarr[2].set_xlim(2E-1, 50)
    axarr[2].set_xscale('log')
    axarr[2].set_xlabel(r'Days since $t_0$', fontsize=16)
    axarr[0].set_ylabel(r'$L_\mathrm{bol}$ (erg/s)', fontsize=16)
    axarr[1].set_ylabel(r'$R_\mathrm{ph}$ (cm)', fontsize=16)
    axarr[2].set_ylabel(r'$T_\mathrm{eff}$ (K)', fontsize=16)

    plt.subplots_adjust(hspace=0)
    plt.tight_layout()
    plt.show()
    #plt.savefig("bbfit_log.png")


if __name__=="__main__":
    plot()
