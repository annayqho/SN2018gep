""" Plot the physical evolution as a power law """

import numpy as np
from math import floor, log10
import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/code")
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from scipy.optimize import curve_fit
from astropy.io import ascii
from astropy.time import Time
from load_lum import load_lc
from load_radius import load_radius
from load_temp import load_temp

datadir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data"


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
        t = tprint[ii]
        ut = utprint[ii]
        lt = ltprint[ii]
        if tprint[ii] >= 10:
            t = int(t)
            ut = int(ut)
            lt = int(lt)
        linestr = "$%s$ & $%s^{%s}_{%s}$ & $%s^{%s}_{%s}$ & $%s^{%s}_{%s}$ \\\ \n" %(
                dtprint[ii], 
                l, "+%s" %ulprint[ii], "-%s" %llprint[ii],
                int(rprint[ii]), "+%s" %int(urprint[ii]), 
                "-%s" %int(lrprint[ii]),
                t, "+%s" %ut, "-%s" %lt)
        outputf.write(linestr)

    outputf.write("\hline \n")
    outputf.write("\end{tabular} \n")
    outputf.write("\label{tab:physevol} \n")
    outputf.write("\end{table} \n")

    outputf.close()


def lum_panel(ax):
    """ Panel showing the luminosity evolution """
    dt, lum, llum, ulum = load_lc()

    # choose just the optical points
    opt = np.logical_or(dt < 0.1, np.logical_and(dt > 0.5, dt < 3.22))
    ax.errorbar(
            dt[opt], lum[opt], yerr=[llum[opt],ulum[opt]], 
            fmt='o', c='lightgrey', lw=0.5)

    # UV + opt points
    ax.errorbar(
            dt[~opt], lum[~opt], yerr=[llum[~opt],ulum[~opt]], ms=8,
            fmt='o', mec='k', mfc='lightgrey', lw=0.5, c='k', 
            label="SN2018gep")
    ax.plot(dt, lum, lw=1, c='k')
    ax.set_ylabel(r'$L_\mathrm{bol}$ (erg/s)', fontsize=16)
    ax.set_yscale('log')
    ax2 = ax.twinx()
    ax2.set_ylabel(
            r"$(L_\odot$)", 
            fontsize=16, rotation=270, labelpad=15.0)
    y_f = lambda y_i: y_i/3.839E33
    ymin, ymax = ax.get_ylim()
    ax2.set_ylim((y_f(ymin), y_f(ymax)))
    ax2.plot([],[])
    ax2.set_yscale('log')
    ax.tick_params(axis='both', labelsize=16)
    ax2.tick_params(axis='both', labelsize=16)


def rad_panel(ax):
    """ Panel showing the radius evolution """
    dt, rad, lrad, urad = load_radius()
    opt = np.logical_or(dt < 0.1, np.logical_and(dt > 0.5, dt < 3.22))
    ax.errorbar(
            dt[opt], rad[opt]/1E15, yerr=[lrad[opt]/1E15,urad[opt]/1E15], 
            fmt='o', mfc='lightgrey', c='lightgrey', lw=0.5)
    ax.errorbar(
            dt[~opt], rad[~opt]/1E15, yerr=[lrad[~opt]/1E15,urad[~opt]/1E15], 
            fmt='o', c='k', lw=0.5, ms=8, mec='k', mfc='lightgrey')
    ax.plot(
            dt, rad/1E15, lw=1, c='k')

    # Plot lines of constant velocity
    xvals = np.linspace(-1, 1E2, 1000)

    # v = 0.1c
    yvals = 3E14 + 0.1 * (3E10) * xvals * 86400
    ax.plot(xvals, yvals/1E15, ls='--', lw=0.5, c='grey')
    ax.text(26, 6.8, 'v=0.1c', fontsize=14, rotation=20)

    # Inset showing the first few days
    axins = inset_axes(
            ax, 1.5, 1, loc=2)
    axins.errorbar(
            dt[opt], rad[opt]/1.496E13, yerr=[lrad[opt]/1E15,urad[opt]/1E15], 
            fmt='o', c='lightgrey', lw=0.5)
    axins.errorbar(
            dt[~opt], rad[~opt]/1.496E13, 
            yerr=[lrad[~opt]/1E15,urad[~opt]/1E15], 
            fmt='o', mec='k', mfc='lightgrey', ms=8, lw=0.5)
    axins.plot(
            dt, rad/1.496E13, lw=1, c='k')
    axins.plot(xvals, yvals/1.496E13, ls='--', lw=0.5, c='grey')
    axins.set_xlim(-0.5,3)
    axins.set_ylim(10,90)
    axins.tick_params(labelsize=12)
    axins.yaxis.tick_right()
    #axins.set_ylabel("($10^{14}$\,cm)", rotation=270, fontsize=12)
    axins.yaxis.set_label_position("right")

    ax.set_ylabel(r'$R_\mathrm{ph}$ ($10^{15}$ cm)', fontsize=16)
    ax.set_ylim(0,10)

    ax2 = ax.twinx()
    ax2.set_ylabel(
            r"(AU)", 
            fontsize=16, rotation=270, labelpad=15.0)
    y_f = lambda y_i: y_i*1E15/1.496E13
    ymin, ymax = ax.get_ylim()
    ax2.set_ylim((y_f(ymin), y_f(ymax)))
    ax2.plot([],[])
    ax.tick_params(axis='both', labelsize=16)
    ax2.tick_params(axis='both', labelsize=16)


def temp_panel(ax):
    """ Panel showing the temp evolution """
    dt, temp, ltemp, utemp = load_temp()
    opt = np.logical_or(dt < 0.1, np.logical_and(dt > 0.5, dt < 3.22))
    ax.errorbar(
            dt[opt], temp[opt], yerr=[ltemp[opt],utemp[opt]], 
            fmt='o', c='lightgrey', lw=0.5)
    ax.errorbar(
            dt[~opt], temp[~opt], yerr=[ltemp[~opt],utemp[~opt]], 
            fmt='o', mec='k', c='k', mfc='lightgrey', lw=0.5)
    ax.plot(dt, temp, c='k', lw=1)
    ax.tick_params(axis='both', labelsize=16)
    ax.axhline(y=5000, c='grey', ls='--', lw=0.5)
    ax.text(0.5, 4500, "5000 K", fontsize=14, verticalalignment='top')


def lum_16asu(ax):
    """ 
    Plot the lum evolution of iPTF16asu for comparison 
    
    Best-fit explosion time:
    2016 May 10.53 +/- 0.17 days
    """
    t0 = Time("2016-05-10").jd + 0.53
    et0 = 0.17
    dat = np.loadtxt(datadir + "/bol_lc/iPTF16asu_bolometric_luminosity.txt",
        delimiter=" ")
    dt = dat[:,0]-t0 
    lum = dat[:,1]
    elum = dat[:,2]
    ax.errorbar(
            dt, lum, xerr=0.17, yerr=elum, 
            marker='s', mec='#e55c30', mfc='white', c='#e55c30', 
            label="iPTF16asu", zorder=0, ls='--', lw=0.5)
    ax.plot(
            dt, lum, c='#e55c30', zorder=0, ls='--', lw=1)


def temp_16asu(ax):
    """ Plot the lum evolution of iPTF16asu for comparison """
    t0 = Time("2016-05-10").jd + 0.53
    et0 = 0.17
    dat = np.loadtxt(datadir + "/bol_lc/iPTF16asu_Temp_Radius.txt",
        delimiter=" ")
    jd = dat[:,0]
    dt = jd-t0
    temp = dat[:,1]
    etemp = dat[:,2]
    choose = etemp > 0
    ax.errorbar(
            dt[choose], temp[choose],
            xerr=[0.17]*len(dt[choose]), yerr=temp[choose], 
            marker='s', mec='#e55c30', mfc='white', c='#e55c30', 
            zorder=0, ls='--', lw=0.5)
    ax.plot(
            dt[choose], temp[choose], c='#e55c30', zorder=0, ls='--', lw=1)


def rad_16asu(ax):
    """ Plot the radius evolution of iPTF16asu for comparison """
    t0 = Time("2016-05-10").jd + 0.53
    et0 = 0.17
    dat = np.loadtxt(datadir + "/bol_lc/iPTF16asu_Temp_Radius.txt",
        delimiter=" ")
    jd = dat[:,0]
    dt = jd-t0
    rad = dat[:,3] # to cm
    erad =dat[:,4]
    choose = erad > 0
    ax.errorbar(
            dt[choose], rad[choose]/1E15, xerr=0.17, yerr=erad[choose]/1E15, 
            marker='s', mec='#e55c30', mfc='white', c='#e55c30', 
            zorder=0, ls='--', lw=0.5)
    ax.plot(
            dt[choose], rad[choose]/1E15, c='#e55c30', zorder=0, ls='--', lw=1)

    # make an inset showing the first few days



def plot(scale='loglinear'):
    # Initialize the figure
    fig,axarr = plt.subplots(3,1, figsize=(7,8), sharex=True)

    # Luminosity panel
    lum_panel(axarr[0])
    lum_16asu(axarr[0])
    axarr[0].legend(fontsize=12) # wait until after 16asu
    if scale=='loglinear':
        axarr[0].set_ylim(1E42, 1E45)
    elif scale=='loglog':
        axarr[0].set_ylim(5E41, 1E45)

    # Radius panel
    rad_panel(axarr[1])
    rad_16asu(axarr[1])

    # Temperature panel
    temp_panel(axarr[2])
    temp_16asu(axarr[2])
    axarr[2].set_yscale('log')
    if scale=='loglog':
        axarr[2].set_xscale('log')
        axarr[2].set_xlim(0.03, 100)
    elif scale=='loglinear':
        axarr[2].set_xlim(-1, 40)

    axarr[0].xaxis.label.set_visible(False)
    axarr[1].xaxis.label.set_visible(False)

    axarr[0].tick_params(axis='y', labelsize=16)
    axarr[1].tick_params(axis='y', labelsize=16)

    axarr[2].set_xlabel(r'Days since $t_0$', fontsize=16)
    axarr[2].set_ylabel(r'$T_\mathrm{eff}$ (K)', fontsize=16)

    #plt.subplots_adjust(hspace=0)
    plt.tight_layout()
    #plt.show()
    plt.savefig("bbfit_%s.eps" %scale, format='eps', dpi=1000)


if __name__=="__main__":
    plot(scale='loglinear')
