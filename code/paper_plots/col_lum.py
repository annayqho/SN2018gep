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

    dt_grid = np.arange(0,60,1)

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
            0, 0, marker='o', c='k', label="ZTF18abukavn (AT2018gep)")
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

    dt_grid = np.arange(0,60,1)

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

            i = 0
            ax.errorbar(
                    gr[i], useg[i]-distmod, 
                    xerr=egr[i], yerr=useeg[i], c='grey', marker='+',
                    alpha=0.3)
            i = -1
            ax.errorbar(
                    gr[i], useg[i]-distmod, 
                    xerr=egr[i], yerr=useeg[i], c='grey', marker='+',
                    alpha=0.3)
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


def arcavi(ax):
    ddir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/lc"
    inputf = ddir + "/arcavi_info.txt"
    info = np.loadtxt(inputf,delimiter=';',dtype=str)
    unames = np.array([val.strip() for val in info[:,0]])
    z_key = info[:,3].astype(float)

    inputf = ddir + "/arcavi_lc.txt"
    dat = np.loadtxt(inputf,delimiter=';',dtype=str)
    names = np.array([val.strip() for val in dat[:,0]])
    for ii,name in enumerate(unames):
        choose = names == name
        redshift = z_key[unames==name]
        filt = np.array(dat[:,2][choose]).astype(str)
        mjd = dat[:,3][choose].astype(float)
        dt = mjd-mjd[0]
        islim = dat[:,4][choose]
        mag = dat[:,5][choose].astype(float)
        emag = dat[:,6][choose].astype(str)
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

        dt_grid = np.sort(np.intersect1d(rdt.astype(int),gdt.astype(int)))

        if len(dt_grid) > 1:
            print(dt_grid[-1])
            # only do this if npts > 1
            keep = np.array(
                    [np.where(rdt.astype(int)==val)[0][0] for val in dt_grid])
            user = r[keep]
            useer = er[keep]
            keep = np.array(
                    [np.where(gdt.astype(int)==val)[0][0] for val in dt_grid])
            useg = g[keep]
            useeg = eg[keep]
            gr = useg-user
            egr = np.sqrt(useeg**2+useer**2)

            col = 'orange'
            i = 0
            ax.errorbar(
                    gr[i], useg[i]-distmod, 
                    xerr=egr[i], yerr=useeg[i], c=col, marker='+',
                    alpha=0.3)
            i = -1
            ax.errorbar(
                    gr[i], useg[i]-distmod, 
                    xerr=egr[i], yerr=useeg[i], c=col, marker='+',
                    alpha=0.3)
            ax.text(gr[0], useg[0]-distmod, name)
            if ii == 2:
                ax.plot(gr, useg-distmod, ls='-', c=col, zorder=0, lw=3, 
                        alpha=0.3, label="SNLS (Arcavi+16)")
                # ax.fill_between(gr, useg-distmod-useeg,
                #     useg-distmod+useeg, color='grey', alpha=0.3, zorder=0,
                #     label="PS-1 gold (Drout+14)")
                #ax.fill_betweenx(useg-distmod, gr-egr, gr+egr,
                #    color='grey', alpha=0.3, zorder=0,
                #    label="PS-1 gold (Drout+14)")
            else:
                ax.plot(gr, useg-distmod, ls='-', c=col, zorder=0, lw=3, 
                        alpha=0.3)


def sn2002bj(ax):
    distmod = Planck15.distmod(z=0.012108).value
    dt = np.array([1, 3, 4, 5, 9, 11, 15])
    B = np.array([14.79, 14.97, 15.03, 15.15, 15.75, 16.27, 17.92])
    V = np.array([14.91, 15.03, 15.10, 15.20, 15.71, 16.04, 17.93])
    R = np.array([14.97, 15.06, 15.08, 15.21, 15.65, 16.06, 17.32])
    ax.plot(B-R, B-distmod, ls='-', c='purple', zorder=0, lw=3, 
            alpha=0.3, label="SN2002bj")


def fbot():
    fig,ax = plt.subplots(1,1,figsize=(9,7))
    cb = at2018gep(ax)
    cbar = plt.colorbar(cb)
    at2018cow(ax)
    drout(ax)
    arcavi(ax)
    #sn2002bj(ax) # it's a spec. classified SN, not an FBOT
    ax.tick_params(axis='both', labelsize=14)
    cbar.ax.set_ylabel("Days from some $t_0$", fontsize=12)
    cbar.ax.tick_params(labelsize=12)
    ax.set_xlabel("$g-r$, observer frame", fontsize=16)
    ax.set_ylabel("Absolute $g$-band mag, observer frame", fontsize=16)
    plt.xlim(-1, 1.5)
    plt.ylim(-12.5, -21)
    plt.legend(prop={'size':12})

    plt.tight_layout()
    #plt.show()
    plt.savefig("fbot_g_gr.png")


def asu(ax):
    ddir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/lc"
    dat = np.loadtxt(ddir + "/lc_16asu.txt", delimiter='&', dtype='str')
    t = dat[:,0].astype(float)
    dt = t-t[0]
    band = np.array([val.strip() for val in dat[:,2]])
    mag_raw = dat[:,3]
    mag = np.zeros(len(mag_raw))
    emag = np.zeros(len(mag))
    for ii,val in enumerate(mag_raw):
        if '>' not in val:
            mag[ii] = float(val.split('$pm$')[0])
            emag[ii] = float(val.split('$pm$')[1])

    choose = np.logical_and(mag > 0, band == 'g')
    g = mag[choose]
    gdt = dt[choose]
    choose = np.logical_and(mag > 0, band == 'r')
    r = mag[choose]
    rdt = dt[choose]
      
    dt_grid = np.arange(20, 36, 1)
    gplt = np.interp(dt_grid, gdt, g) 
    rplt = np.interp(dt_grid, rdt, r)
    gr = gplt-rplt
    ax.plot(
            gr, gplt-39.862, c='grey', 
            lw=3, alpha=0.3, label="iPTF16asu")
    ax.errorbar(


if __name__=="__main__":
    fig,ax = plt.subplots(1,1,figsize=(9,7))
    cb = at2018gep(ax)
    asu(ax)

    # Formatting
    cbar = plt.colorbar(cb)
    ax.tick_params(axis='both', labelsize=14)
    cbar.ax.set_ylabel("Days from some $t_0$", fontsize=12)
    cbar.ax.tick_params(labelsize=12)
    ax.set_xlabel("$g-r$, observer frame", fontsize=16)
    ax.set_ylabel("Absolute $g$-band mag, observer frame", fontsize=16)
    plt.xlim(-1, 1.5)
    plt.ylim(-12.5, -21)
    plt.legend(prop={'size':12})

    plt.tight_layout()

    plt.show()
    #plt.savefig("icbl_g_gr.png")

