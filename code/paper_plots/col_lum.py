""" Color evolution, luminosity evolution phase space """

import numpy as np
import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from astropy.cosmology import Planck15


def at2018gep(ax):
    distmod = Planck15.distmod(z=0.033).value
    DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/phot"

    f = DATA_DIR + "/ZTF18abukavn_opt_phot.dat"
    dat = np.loadtxt(f, dtype=str, delimiter=' ')
    instr = dat[:,0]
    jd = dat[:,1].astype(float)
    filt = dat[:,2]
    mag = dat[:,3].astype(float) 
    emag = dat[:,4].astype(float)

    det = np.logical_and(mag<99, ~np.isnan(mag))
    nondet = np.logical_or(mag==99, np.isnan(mag))
    zp = 2458370.6473
    dt = jd-zp

    dt_grid = np.arange(0,30,1)

    band = filt=='g'
    choose = np.logical_and(det, band)
    order = np.argsort(dt[choose])
    g = np.interp(dt_grid, dt[choose][order], mag[choose][order])

    # To get the uncertainties, for each dt point, choose the 5 points
    # closest in time. Take the median uncertainty.
    gerr = np.zeros(len(g))
    for ii,gval in enumerate(g):
        pts = np.argsort(np.abs(dt[choose]-dt_grid[ii]))[0:5]
        gerr[ii] = np.median(emag[choose][pts])

    # Later I should probably do a proper perturbation -> actual
    # measurement of the uncertainty with a Monte Carlo,
    # but I don't feel like doing that right now.

    band = filt=='r'
    choose = np.logical_and(det, band)
    order = np.argsort(dt[choose])
    r = np.interp(dt_grid, dt[choose][order], mag[choose][order])

    rerr = np.zeros(len(r))
    for ii,rval in enumerate(r):
        pts = np.argsort(np.abs(dt[choose]-dt_grid[ii]))[0:5]
        rerr[ii] = np.median(emag[choose][pts])

    gr = g - r
    gr_err = np.sqrt(gerr**2+rerr**2)
    cb = ax.errorbar(
            gr, g-distmod, xerr=gr_err, yerr=gerr, c='k', fmt='.')
    cb = ax.scatter(
            gr, g-distmod, c=dt_grid, cmap='inferno', marker='o', zorder=5)
    ax.scatter(
            0, 0, marker='o', c='k', label="AT2018gep")
    return cb


def at2018cow(ax):
    distmod = Planck15.distmod(z=0.014).value
    DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/lc"
    f = DATA_DIR + "/18cow.txt"
    dat = np.loadtxt(f, dtype=str)
    jd = dat[:,0].astype(float)
    filt = dat[:,2]
    mag = dat[:,4].astype(float)
    emag = dat[:,5].astype(float)

    zp = 58285
    dt = jd-zp

    dt_grid = np.arange(0,30,1)

    band = filt=='g'
    choose = band
    order = np.argsort(dt[choose])
    g = np.interp(dt_grid, dt[choose][order], mag[choose][order])

    # To get the uncertainties, for each dt point, choose the 5 points
    # closest in time. Take the median uncertainty.
    gerr = np.zeros(len(g))
    for ii,gval in enumerate(g):
        pts = np.argsort(np.abs(dt[choose]-dt_grid[ii]))[0:5]
        gerr[ii] = np.median(emag[choose][pts])

    # Later I should probably do a proper perturbation -> actual
    # measurement of the uncertainty with a Monte Carlo,
    # but I don't feel like doing that right now.

    band = filt=='r'
    choose = band
    order = np.argsort(dt[choose])
    r = np.interp(dt_grid, dt[choose][order], mag[choose][order])

    rerr = np.zeros(len(r))
    for ii,rval in enumerate(r):
        pts = np.argsort(np.abs(dt[choose]-dt_grid[ii]))[0:5]
        rerr[ii] = np.median(emag[choose][pts])

    gr = g - r
    gr_err = np.sqrt(gerr**2+rerr**2)
    ax.errorbar(
            gr, g-distmod, xerr=gr_err, yerr=gerr, c='k', fmt='.', zorder=0)
    ax.scatter(
            gr, g-distmod, c=dt_grid, cmap='inferno', marker='s', zorder=2)
    ax.scatter(
            0, 0, marker='s', c='k', label="AT2018cow")


def drout(ax):
    ddir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/lc"
    z = {'10ah':0.074, '10bjp':0.113, '11qr':0.324, '12bb':0.101, '12bv':0.405,
         '12brf':0.275, '11bbq': 0.646, 
         '13duy':0.270, '13dwm':0.245, '13ess':0.296}
    inputf = ddir + "/drout14.txt"
    dat = np.loadtxt(inputf,delimiter=';',dtype=str)
    names = np.array([val.strip() for val in dat[:,0]])
    unames = np.unique(names)
    for ii,name in enumerate(unames):
        print(name)
        choose = names == name
        redshift = z[name]
        # only do this if there is a known redshift
        filt = np.array(dat[:,1][choose]).astype(str)
        dt = dat[:,2][choose].astype(float)
        islim = dat[:,3][choose]
        mag = dat[:,4][choose].astype(float)
        emag = dat[:,5][choose].astype(str)
        distmod = Planck15.distmod(z=redshift).value
        isdet = np.array([val == " " for val in islim])
        isfilt = np.array(['g' in val for val in filt])
        pts = np.logical_and(isdet, isfilt)
        gdt = dt[pts]
        g = mag[pts]
        eg = emag[pts].astype(float)
        isfilt = np.array(['r' in val for val in filt])
        pts = np.logical_and(isdet, isfilt)
        rdt = dt[pts]
        r = mag[pts]
        er = emag[pts].astype(float)

        dt_grid = np.sort(np.intersect1d(rdt,gdt))

        if len(dt_grid) > 1:
            # only do this if npts > 1
            keep = np.array([np.where(rdt==val)[0][0] for val in dt_grid])
            user = r[keep]
            useer = er[keep]
            keep = np.array([np.where(gdt==val)[0][0] for val in dt_grid])
            useg = g[keep]
            useeg = eg[keep]
            gr = useg-user
            egr = np.sqrt(useeg**2+useer**2)

            ax.scatter(gr[0], useg[0]-distmod, c='grey')
            ax.scatter(gr[-1], useg[-1]-distmod, c='grey')
            ax.text(gr[0], useg[0]-distmod, name)
            if ii == 0:
                ax.plot(gr, useg-distmod, ls='-', c='grey', zorder=0, lw=3, 
                        alpha=0.5, label="PS-1 gold (Drout+14)")
                # ax.fill_between(gr, useg-distmod-useeg,
                #     useg-distmod+useeg, color='grey', alpha=0.3, zorder=0,
                #     label="PS-1 gold (Drout+14)")
                #ax.fill_betweenx(useg-distmod, gr-egr, gr+egr,
                #    color='grey', alpha=0.3, zorder=0,
                #    label="PS-1 gold (Drout+14)")
            else:
                ax.plot(gr, useg-distmod, ls='-', c='grey', zorder=0, lw=3, 
                        alpha=0.5)
                # ax.fill_between(gr, useg-distmod-useeg,
                #     useg-distmod+useeg, color='grey', alpha=0.3, zorder=0)
                #ax.fill_betweenx(useg-distmod, gr-egr, gr+egr,
                #    color='grey', alpha=0.3, zorder=0)


if __name__=="__main__":
    fig,ax = plt.subplots(1,1,figsize=(8,5))
    cb = at2018gep(ax)
    cbar = plt.colorbar(cb)
    at2018cow(ax)
    drout(ax)
    ax.tick_params(axis='both', labelsize=14)
    cbar.ax.set_ylabel("Days from some $t_0$", fontsize=12)
    cbar.ax.tick_params(labelsize=12)
    ax.set_xlabel("$g-r$, observer frame", fontsize=16)
    ax.set_ylabel("Absolute $g$-band mag, observer frame", fontsize=16)
    plt.xlim(-1, 1.2)
    plt.ylim(-15.3, -21)
    plt.legend(prop={'size':12})

    plt.tight_layout()

    plt.show()
    #plt.savefig("g_gr.png")
