""" Christoffer's host-subtracted photometry
Plot the light curve! """

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15
import glob
from uv_lc import get_uv_lc

zp = 2458370.6473
# extinction values
bands = ['UVW2', 'UVM2', 'UVW1', 'U', 'u', 'B', 'g', 'V', 'r', 'i', 'z']
ext = {}
ext['UVW2'] = 0.05 # placeholder; 1928
ext['UVM2'] = 0.05 # placeholder; 2246
ext['UVW1'] = 0.05 # placeholder; 2600
ext['U'] = 0.05 # placeholder; 3465
ext['u'] = 0.045 # 3543
ext['B'] - 0.04 # placeholder; 4392
ext['g'] = 0.035 # 4770
ext['V'] = 0.03 # placeholder; 5468
ext['r'] = 0.024 # 6231
ext['i'] = 0.018 # 7625
ext['z'] = 0.013 # 9134


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

    axins.set_xlim(-0.3,3)
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
    dt = np.append(dt, add_dt)
    filt = np.append(filt, add_filt)
    add_mag = -2.5 * np.log10(fnu_mjy*1E-3) + 8.90
    mag = np.append(mag, add_mag)
    add_emag = (efnu_mjy/fnu_mjy) # I think it's just the ratio
    emag = np.append(emag, add_emag)

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
                emag[choose][order], c='black', fmt='o', ms=4,
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
                2458370.6408-zp, 19.97, 0, 0.5, length_includes_head=True,
                head_width=2, head_length=0.2, fc='k', ec='k')
        ax.yaxis.set_tick_params(labelsize=14)
        ax.xaxis.set_tick_params(labelsize=14)

        # for each panel, also show absolute mag
        if ii % 2 != 0:
            ax2 = ax.twinx()
            ax2.set_ylabel(
                    "Absolute Magnitude",
                    fontsize=14, rotation=270, labelpad=15.0)
            y_f = lambda y_i: y_i-Planck15.distmod(z=0.033).value
            ymin, ymax = ax.get_ylim()
            ax2.set_ylim((y_f(ymin), y_f(ymax)))
            ax2.plot([],[])
            ax2.tick_params(axis='both', labelsize=14)

        # label with the band
        ax.text(
                0.9, 0.9, "$%s$-band" %use_f, 
                fontsize=14, transform=ax.transAxes,
                horizontalalignment='right')


    # Final reconfiguring
    plt.subplots_adjust(hspace=0, wspace=0)
    ax.invert_yaxis()
    fig.text(0.5, 0.04, 
        r"Days since $t_0=$JD 2458370.6473 (UT 2018 Sept 09.15)", 
        ha='center', fontsize=16) 
    fig.text(0.04, 0.5, 'Apparent Mag', fontsize=16, rotation='vertical')

    #plt.savefig("lc.png")
    plt.show()


if __name__=="__main__":
    plot_lc()
