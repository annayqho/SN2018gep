""" Fit an engine model to the bolometric light curve """

import numpy as np
import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from scipy.optimize import curve_fit
from load_lum import load_lc


def func(t, Ee, te, alpha):
    """ Takes t in sec, Ee in erg, te in sec """
    lum = Ee * (alpha-1) / (te * (1+t/te)**alpha)
    return lum



def logfunc(t, Ee, te, alpha):
    """ Takes t in 10^4 sec, Ee in 10^51 erg, te in 10^4 sec """
    loglum = 51 + np.log10(Ee) + np.log10(alpha-1) \
             - 4 - np.log10(te) - alpha*np.log10(1+t/te)
    return loglum


if __name__=="__main__":
    dt_all, lum_all, llum, ulum = load_lc()

    # Ignore the first point
    dt = dt_all[1:]
    lum = lum_all[1:]

    dt_sec = dt * 86400

    # sample lum from its distribution where sigma = 20% of lum
    npts = len(dt)
    ndraw = 10

    # Make ndraw new versions of a new lum array, repopulated with the lum
    # from the Gaussian dist with STD of ulum

    lum_new = np.zeros((npts, ndraw))
    # For each version, do a new fit
    # and save the resulting Ee, te, alpha
    for ii,lumval in enumerate(lum):
        # log uncertainty is dx/x
        lum_new[ii,:] = np.random.normal(
                loc=lumval, scale=0.20*lumval, size=ndraw)

    Ees = np.zeros(ndraw)
    tes = np.zeros(ndraw)
    alphas = np.zeros(ndraw)

    for ii in np.arange(ndraw):
        Ee0 = 1.48E50
        te0 = 0.82*86400
        alpha0 = 1.9
        p0 = np.array([Ee0,te0,alpha0])
        logp0 = np.array([Ee0/1E51,te0/1E4,alpha0])
        popt, pcov = curve_fit(
                logfunc, dt_sec/1E4, np.log10(lum_new[:,ii]), p0=logp0,
                sigma = [0.2]*len(dt_sec), 
                absolute_sigma = True,
                maxfev=10000,
                bounds=([1E49/1E51,0.1*86400/1E4,1], [1E52/1E51,2*86400/1E4,3]))
        Ees[ii] = popt[0]
        tes[ii] = popt[1] 
        alphas[ii] = popt[2]
        plt.errorbar(dt, lum_new[:,ii], yerr=0.2*lum_new[:,ii], c='k', fmt='s')
        tplot = np.linspace(3,40,1000)
        yplot = 10**logfunc(tplot*86400/1E4, popt[0], popt[1], popt[2])
        plt.plot(tplot, yplot, c='k', ls='--')
        plt.yscale('log')
        plt.savefig("iter%s.png" %ii)
        plt.close()


    # Measure the median and standard deviation 
    plt.hist(Ees)
    plt.savefig("E_hist.png")
    plt.close()
    plt.hist(tes)
    plt.savefig("t_hist.png")
    plt.close()
    plt.hist(alphas)
    plt.savefig("alpha_hist.png")
    plt.close()
    Ee = np.median(Ees) 
    te = np.median(tes) 
    alpha = np.median(alphas)

    Eestr = r"%s \times 10^{50}" %str(np.round(Ee*1E51/1E50,1))
    testr = r"%s" %str(np.round(te*1E4/86400,2))
    alphastr = r"%s" %str(np.round(alpha,1))
    str1 = r"$L_e(t) = \frac{E_e}{t_e} \frac{(\alpha-1)}{(1+t/t_e)^\alpha}$"
    str2 = r"$E_e = %s$ erg, $t_e = %s$ days, $\alpha=%s$" %(
            Eestr,testr,alphastr)

    plt.errorbar(dt_all, lum_all, yerr=[llum,ulum], c='k', fmt='s')
    tplot = np.linspace(3,40,1000)
    lplot = 10**logfunc(tplot*86400/1E4, Ee, te, alpha)
    lplot = func(tplot*86400, 1.5E50, 0.82*86400, 1.9)
    plt.plot(tplot, lplot, c='k', ls='--')
    plt.yscale('log')
    plt.tick_params(axis='both', labelsize=12)
    plt.xlabel(r'Days since $t_0$', fontsize=16)
    plt.ylabel(r'$L_\mathrm{bol}$ (erg/s)', fontsize=16)
    plt.text(
            35, 5E44, str1, horizontalalignment='right', fontsize=12)
    plt.text(
            35, 3E44, str2, horizontalalignment='right', fontsize=12)


