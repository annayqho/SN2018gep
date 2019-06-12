""" Plot full light curves, one panel per band
Advice on aesthetic from Erik Petigura

this is Fig 3 in the paper """


import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15
import glob
import extinction
from uv_lc import get_uv_lc

zp = 2458370.6473


def plot_inset():
    # zoomed-in window showing the earliest non-detection and detection
    axins = inset_axes(
            ax, 2, 1, loc=1,
            bbox_to_anchor=(0.87,0.98),
            bbox_transform=ax.transAxes)
    choose = np.logical_and(det, band)
    axins.errorbar(
        dt[choose]*24, mag[choose], emag[choose], fmt='s', ms=6,
        mec=rcol, mfc=rcol, c=rcol, label='r', zorder=9)
    choose = np.logical_and(nondet, band)
    axins.arrow(
            2458370.6408-zp, 19.97, 0, 0.5, length_includes_head=True,
            head_width=0.2, head_length=0.3, fc='k', ec='k')
    band = filt=='g'
    choose = np.logical_and(np.logical_and(det, band), dt*24 < 3)
    axins.errorbar(
            dt[choose]*24, mag[choose], emag[choose],
            fmt='o', ms=5, mec='#57106e', mfc='white', c='#57106e', label='g')

    # fit a line to this early g-band data
    out = np.polyfit(dt[choose]*24, mag[choose], deg=1, w=1/emag[choose])
    m,b = out
    dt_plt = np.linspace(-1,3)
    y_plt = m*dt_plt + b
    axins.plot(dt_plt, y_plt, ls='--', c='k', lw=0.5)
    axins.text(0.5, 0.5, "31.2 mag/day", fontsize=12, transform=axins.transAxes,
            verticalalignment='top')

    axins.set_xlim(-0.1,3)
    axins.set_ylim(18,21)
    axins.tick_params(axis='both', labelsize=12)
    axins.set_xlabel(r"Hours since $t_0$", fontsize=12)
    axins.invert_yaxis()
    ax.plot([-1, -1], [21, 18], c='k', lw=0.5)
    ax.plot([1, 1], [21, 18], c='k', lw=0.5)
    ax.plot([-1, 1], [18, 18], c='k', lw=0.5)
    ax.plot([-1, 1], [21, 21], c='k', lw=0.5)
    #mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
     

def get_lc():
    # get optical light curves 
    DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/phot"
    f = DATA_DIR + "/ZTF18abukavn_opt_phot.dat"
    dat = np.loadtxt(f, dtype=str, delimiter=' ')
    instr = dat[:,0]
    jd = dat[:,1].astype(float)
    filt = dat[:,2]
    mag = dat[:,3].astype(float)
    emag = dat[:,4].astype(float)
    dt = jd-zp

    # add the UV light curves
    add_dt, add_filt, fnu_mjy, efnu_mjy = get_uv_lc()
    # convert to AB mag
    add_mag = -2.5 * np.log10(fnu_mjy*1E-3) + 8.90
    add_emag = (efnu_mjy/fnu_mjy) # I think it's just the ratio
    choose = add_emag < 50
    dt = np.append(dt, add_dt[choose])
    filt = np.append(filt, add_filt[choose])
    mag = np.append(mag, add_mag[choose])
    emag = np.append(emag, add_emag[choose])

    return dt, filt, mag, emag


def plot_lc():
    dt, filt, mag, emag = get_lc()
    det = np.logical_and(mag<99, ~np.isnan(mag))
    nondet = np.logical_or(mag==99, np.isnan(mag))

    fig,axarr = plt.subplots(
            4, 3, figsize=(8,8), sharex=True, sharey=True)

    for ii,use_f in enumerate(bands):
        ax = axarr.reshape(-1)[ii]
        choose = np.logical_and(det, filt == use_f)
        order = np.argsort(dt[choose])
        ax.errorbar(
                dt[choose][order], mag[choose][order]-ext[use_f], 
                emag[choose][order], c='black', fmt='o', ms=3,
                alpha=1.0, zorder=5)
        # for each panel, plot all of them as a grey background
        for f in bands:
            choose = np.logical_and(det, filt == f)
            order = np.argsort(dt[choose])
            ax.plot(
                    dt[choose][order], mag[choose][order]-ext[f], 
                    c='lightgrey', alpha=0.7, zorder=0)

        # for each panel, show the last non-detection (which was r-band)
        ax.arrow(
                2458370.6408-zp, 20.47, 0, 0.5, length_includes_head=True,
                head_width=2, head_length=0.2, fc='r', ec='r', zorder=10)
        # ax.arrow(
        #         2458370.6408-zp, 19.97, 0, 0.5, length_includes_head=True,
        #         head_width=0.05, head_length=0.2, fc='k', ec='k')

        # for each panel, show the r-band peak with a cross
        print(2458374.65-zp)
        ax.scatter(
                2458374.65-zp, 16.3, marker='+', c='r', zorder=10) 

        ax.yaxis.set_tick_params(labelsize=14)
        ax.xaxis.set_tick_params(labelsize=14)

        # for each panel, also show absolute mag
        # if ii % 3 == 2:
        #     ax2 = ax.twinx()
        #     y_f = lambda y_i: y_i-Planck15.distmod(z=0.03154).value
        #     ymin, ymax = ax.get_ylim()
        #     ax2.set_ylim((y_f(ymin), y_f(ymax)))
        #     ax2.plot([],[])
        #     ax2.tick_params(axis='both', labelsize=14)

        # label with the band
        if use_f == 'r':
            usecol = 'red'
            wid = 1.0
        else:
            usecol = 'k'
            wid = 0.5
        ax.text(
                0.9, 0.9, "$%s$" %use_f, 
                fontsize=14, transform=ax.transAxes,
                horizontalalignment='right',
                verticalalignment='top',
                bbox=dict(
                    boxstyle="round", fc='white', ec=usecol, 
                    lw=wid, alpha=1.0, pad=0.4))


    # Final reconfiguring
    axarr.reshape(-1)[-1].set_visible(False)
    plt.subplots_adjust(hspace=0, wspace=0)
    #ax.set_xlim(1,1)
    ax.set_xlim(-5, 70)
    ax.set_ylim(15, 21)
    ax.invert_yaxis()
    fig.text(0.5, 0.04, 
        r"Days since $t_0=$JD 2458370.6473 (UT 2018 Sept 09.15)", 
        ha='center', fontsize=16) 
    fig.text(0.04, 0.5, 'Apparent Mag (AB)', fontsize=16, rotation='vertical')
    #fig.text(0.9, 0.5, 'Absolute Mag', fontsize=16, rotation=270)

    plt.savefig("lc.eps", format='eps', dpi=1000)
    #plt.show()


if __name__=="__main__":
    plot_lc()
