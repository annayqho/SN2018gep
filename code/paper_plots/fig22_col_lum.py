""" Color evolution, luminosity evolution phase space """

import numpy as np
import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from astropy.io import ascii
from astropy.cosmology import Planck15
from astropy.time import Time


def at2018gep(ax):
    distmod = Planck15.distmod(z=0.03154).value
    DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/phot"
    t0 = 2458370.6473

    # For this, only use the P48 g & r forced photometry
    f = DATA_DIR + "/precursor.csv"
    prog = ascii.read(f)
    mjd = prog['mjd']
    jd = mjd + 2400000.5
    dt = jd - t0
    filt = prog['filter']
    mag = prog['mag']
    emag = prog['magerr']
    code = prog['instrument']
    det = np.logical_and.reduce((dt > 0, ~np.isnan(mag), code=='ZTF Camera'))

    # Interpolate everything onto the r-band light curve
    band = filt == 'ztfr'
    choose = np.logical_and(det, band)
    order = np.argsort(dt[choose])
    dt_r = dt[choose][order]
    mag_r = mag[choose][order]

    band = filt=='ztfg'
    choose = np.logical_and(det, band)
    order = np.argsort(dt[choose])
    dt_g = dt[choose][order]
    mag_g = mag[choose][order]

    tgrid = np.copy(dt_g)
    # keep all the points below 1 days, and then only choose 1 pt per day
    tgrid = np.hstack((tgrid[tgrid<1], np.arange(1,30,1)))
    #tgrid = np.logspace(np.log10(tgrid[0]), np.log10(tgrid[-1]), 20)

    r = np.interp(tgrid, dt_r, mag_r)
    g = np.interp(tgrid, dt_g, mag_g)

    gr = g - r
    xdata = gr
    ydata = g-distmod
    ax.plot(xdata, ydata, c='k', zorder=5, lw=2)

    markt = np.array([1/24, 1, 10, 20])
    labs = np.array(["1 hour", "1 day", "10 days", "20 days"])
    cols = np.array(['#f6d746', '#e55c30', '#84206b', '#140b34'])
    for ii,t in enumerate(markt):
        prevptx = np.interp(t-0.1, tgrid, xdata)#xdata[tgrid<t][-1]
        prevpty = np.interp(t-0.1, tgrid, ydata)#ydata[tgrid<t][-1]
        newptx = np.interp(t,tgrid,xdata)
        newpty = np.interp(t,tgrid,ydata)
        ax.annotate('', 
            xytext=(prevptx, prevpty),
            xy=(newptx, newpty),
            arrowprops=dict(color=cols[ii], width=1, headlength=10),
            label=labs[ii], zorder=10)
        ax.scatter(0,0,marker='^',c=cols[ii],s=100, label=labs[ii])
    ax.legend(loc='upper right')
        # ax.text(newptx, newpty, "$\Delta t$= %s" %labs[ii], fontsize=12,
        #         horizontalalignment='left', verticalalignment='top')
    ax.text(
            -0.3, -16, "SN2018gep", fontsize=14,
            horizontalalignment='right')


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
    print(dt)

    tgrid = np.arange(dt[0],30,1)

    band = filt=='g'
    choose = band
    order = np.argsort(dt[choose])
    g = np.interp(tgrid, dt[choose][order], mag[choose][order])

    band = filt=='r'
    choose = band
    order = np.argsort(dt[choose])
    r = np.interp(tgrid, dt[choose][order], mag[choose][order])

    gr = g - r
    xdata = gr
    ydata = g-distmod
    ax.plot(
            xdata, ydata, c='k', ls='--', lw=2, zorder=2)
    markt = np.array([1.8, 10, 20])
    labs = np.array(["1 day", "10 days", "20 days"])
    cols = np.array(['#e55c30', '#84206b', '#140b34'])
    for ii,t in enumerate(markt):
        prevptx = np.interp(t-0.1, tgrid, xdata)#xdata[tgrid<t][-1]
        prevpty = np.interp(t-0.1, tgrid, ydata)#ydata[tgrid<t][-1]
        newptx = np.interp(t,tgrid,xdata)
        newpty = np.interp(t,tgrid,ydata)
        ax.annotate('', 
            xytext=(prevptx, prevpty),
            xy=(newptx, newpty),
            arrowprops=dict(color=cols[ii], width=1, headlength=10))
        # ax.text(newptx, newpty, "$\Delta t$= %s" %labs[ii], fontsize=12,
        #         horizontalalignment='left', verticalalignment='top')
    ax.text(-0.37, -20.4, "AT2018cow", fontsize=14)


def drout_all(ax):
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


def ps1_10ah(ax):
    """ the only Drout gold transient with z < 0.1 """
    gr = np.array([0, -0.2, 0.05])
    phase = np.array([-3, 0, 14])
    g = np.array([22.18, 19.95, 21.64])
    distmod = Planck15.distmod(z=0.074).value
    G = g-distmod
    ax.plot(gr, G, ls='-', c='grey', zorder=0, lw=3, alpha=0.5)
    ax.text(gr[0], G[0], 'PS1-10ah', fontsize=12)


def ksn2015k(ax):
    diff = Planck15.distmod(z=0.09).value
    gr = -0.17
    G = 20.01-diff
    ax.errorbar(
            gr, G, xerr=0.20, yerr=0.12, 
            fmt='v', c='#84206b', mfc='#84206b', ms=12, label=None)
    ax.text(gr, G, "KSN2015K", fontsize=14,
            horizontalalignment='right', verticalalignment='bottom')


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
    fig,ax = plt.subplots(1,1,figsize=(7,5))
    at2018gep(ax)
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
            gr, gplt-39.862, c='black', 
            lw=3, alpha=0.3, label="No radio, no GRB")
    ax.text(gr[0], gplt[0]-39.862, "iPTF16asu")


def sn2009bb(ax):
    """ got this from the open SN catalog """
    ddir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/lc"
    dat = np.loadtxt(ddir + "/sn2009bb.txt", delimiter=',', dtype='str')
    mjd = dat[:,1].astype(float)
    mag = dat[:,2].astype(float)
    #emag = dat[:,3].astype(float)
    band = dat[:,5]
    dt = mjd-mjd[0]
    choose = band == 'g'
    g = mag[choose]
    gdt = dt[choose]

    choose = np.logical_or(band == 'r', band == 'R')
    r = mag[choose]
    rdt = dt[choose]
      
    distmod = Planck15.distmod(z=0.0104).value
    dt_grid = np.arange(0, 50, 1)
    gplt = np.interp(dt_grid, gdt, g) 
    rplt = np.interp(dt_grid, rdt, r)
    gr = gplt-rplt
    ax.plot(
            gr, gplt-distmod, c='orange', 
            lw=3, alpha=0.3, label="Radio, no GRB")
    ax.text(gr[0], gplt[0]-distmod, "SN2009bb")


def sn2012ap(ax):
    """ got this from the open SN catalog """
    ddir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/lc"
    dat = np.loadtxt(ddir + "/sn2012ap.txt", delimiter=',', dtype='str')
    mjd = dat[:,1].astype(float)
    mag = dat[:,2].astype(float)
    band = dat[:,5]
    dt = mjd-mjd[0]

    choose = band == 'B'
    g = mag[choose]
    gdt = dt[choose]

    choose = band == 'R'
    r = mag[choose]
    rdt = dt[choose]
      
    distmod = Planck15.distmod(z=0.01224).value
    dt_grid = np.arange(0, 18, 1)
    gplt = np.interp(dt_grid, gdt, g) 
    rplt = np.interp(dt_grid, rdt, r)
    gr = gplt-rplt
    ax.plot(
            gr, gplt-distmod, c='orange', 
            lw=3, alpha=0.3)
    ax.text(gr[0], gplt[0]-distmod, "SN2012ap")


def iptf17cw(ax):
    # R-band 
    g = np.array([19.041, 19.235, 19.294, 20.93])
    r = np.array([18.703, 18.721, 18.793, 19.998])
    distmod = Planck15.distmod(z=0.093).value
    gr = g-r
    ax.plot(
            gr, g-distmod, c='orange', 
            lw=3, alpha=0.3)
    ax.text(gr[0], g[0]-distmod, "iPTF17cw")


def sn2006aj(ax):
    t0 = Time('2006-02-18').mjd + 0.149
    t = np.array(
            [Time('2006-02-20').mjd + 0.162,
             Time('2006-02-21').mjd + 0.196,
             Time('2006-02-23').mjd + 0.196,
             Time('2006-02-25').mjd + 0.112,
             Time('2006-02-27').mjd + 0.167,
             Time('2006-03-04').mjd + 0.097])
    dt = t-t0
    V = np.array([17.84, 17.73, 17.30, 17.10, 17.0, 17.09])
    R = np.array([17.76, 17.22, 17.22, 16.96, 16.84, 16.88])
    gr = V-R
    distmod = Planck15.distmod(z=0.033023).value
    plt.scatter(
            gr, V-distmod, marker='o', edgecolor='orange', 
            facecolor='white', label='Ic-BL')
    plt.plot(gr, V-distmod, c='orange')


def sn1998bw(ax):
    """ got this from the open SN catalog """
    ddir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/lc"
    dat = np.loadtxt(ddir + "/sn1998bw.txt", delimiter=',', dtype='str')
    mjd = dat[:,1].astype(float)
    mag = dat[:,2].astype(float)
    band = dat[:,5]
    dt = mjd-mjd[0]

    choose = band == 'V'
    g = mag[choose]
    gdt = dt[choose]

    choose = np.logical_or(band == 'R', band == "r'")
    r = mag[choose]
    rdt = dt[choose]
      
    distmod = Planck15.distmod(z=0.0085).value
    dt_grid = np.arange(0, 50, 1)
    gplt = np.interp(dt_grid, gdt, g) 
    rplt = np.interp(dt_grid, rdt, r)
    gr = gplt-rplt
    ax.plot(
            gr, gplt-distmod, c='purple', 
            lw=3, alpha=0.3)
    ax.text(gr[1], gplt[1]-distmod, "SN1998bw",
            horizontalalignment='right')


def sn2010bh(ax):
    """ got this from the open SN catalog """
    ddir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/lc"
    dat = np.loadtxt(ddir + "/sn2010bh.txt", delimiter=',', dtype='str')
    mjd = dat[:,1].astype(float)
    mag = dat[:,2].astype(float)
    band = dat[:,5]
    dt = mjd-mjd[0]

    choose = band == 'V'
    g = mag[choose]
    gdt = dt[choose]

    choose = np.logical_or(band == 'R_c', band == "r'")
    r = mag[choose]
    rdt = dt[choose]
      
    distmod = Planck15.distmod(z=0.0593).value
    dt_grid = np.arange(5, 34, 1)
    gplt = np.interp(dt_grid, gdt, g) 
    rplt = np.interp(dt_grid, rdt, r)
    gr = gplt-rplt
    ax.plot(
            gr, gplt-distmod, c='purple', 
            lw=3, alpha=0.3)
    ax.text(gr[1], gplt[1]-distmod, "SN2010bh",
            horizontalalignment='right')


def sn2007ru(ax):
    """ got this from the open SN catalog """
    ddir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/lc"
    dat = np.loadtxt(ddir + "/sn2007ru.txt", delimiter=',', dtype='str')
    mjd = dat[:,1].astype(float)
    mag = dat[:,2].astype(float)
    band = dat[:,5]
    dt = mjd-mjd[0]

    choose = band == 'V'
    g = mag[choose]
    gdt = dt[choose]

    choose = np.logical_or(band == 'R_c', band == "r'")
    r = mag[choose]
    rdt = dt[choose]
      
    distmod = Planck15.distmod(z=0.0593).value
    dt_grid = np.arange(5, 34, 1)
    gplt = np.interp(dt_grid, gdt, g) 
    rplt = np.interp(dt_grid, rdt, r)
    gr = gplt-rplt
    ax.plot(
            gr, gplt-distmod, c='purple', 
            lw=3, alpha=0.3)
    ax.text(gr[1], gplt[1]-distmod, "SN2010bh",
            horizontalalignment='right')


if __name__=="__main__":
    fig,ax = plt.subplots(1,1,figsize=(7,5))
    cb = at2018gep(ax)
    at2018cow(ax)
    ksn2015k(ax)

    # Formatting
    ax.tick_params(axis='both', labelsize=14)
    ax.set_xlabel("$g-r$, observer frame", fontsize=16)
    ax.set_ylabel("Absolute $g$-band mag, observer frame", fontsize=16)
    plt.xlim(-0.7, 1.0)
    plt.ylim(-15.8, -20.7)
    plt.legend(prop={'size':14})

    plt.tight_layout()

    #plt.show()
    plt.savefig("g_gr.eps", format="eps", dpi=1000)
