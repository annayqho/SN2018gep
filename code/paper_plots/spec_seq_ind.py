""" Plot the spectral sequence, but save as individual files to make a .gif """

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15
from astropy.time import Time
import glob

files = glob.glob(
"/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/spec/ZTF18abukavn/*.ascii")
files = np.array(files[12:])
dt = np.zeros(len(files))
cols = np.array([""]*len(dt), dtype='U10')

# Read in all of the files, pull out the corresponding dates, and sort by date
t0 = 2458370.6634 # in JD
for ii,f in enumerate(files):
    tel = f.split("_")[2]
    alldat = open(f).readlines()
    if tel == 'LT':
        for line in alldat:
            if 'DATE-OBS' in line:
                obsdate = line[13:36]
                t = Time(obsdate, format='isot').jd
                dt[ii] = t-t0
        cols[ii] = 'magenta'
    elif tel == 'P200':
        for line in alldat:
            if 'UT shutter open' in line:
                obsdate = line[12:35]
                print(obsdate)
                t = Time(obsdate, format='isot').jd
                dt[ii] = t-t0
        cols[ii] = 'lightblue'
    elif tel == 'Keck1':
        for line in alldat:
            if 'DATE_BEG' in line:
                obsdate = line[13:32]
                t = Time(obsdate, format='isot').jd
                dt[ii] = t-t0
        cols[ii] = 'red'
    elif tel == 'DCT':
        obsdate = '2018-09-14T00:00:00' # temporary
        t = Time(obsdate, format='isot').jd
        dt[ii] = t-t0
        cols[ii] = 'yellow'
    elif tel == 'NOT':
        obsdate = '2018-09-17T00:00:00' # temporary
        t = Time(obsdate, format='isot').jd
        dt[ii] = t-t0
        cols[ii] = 'green'
    elif tel == 'P60':
        for line in alldat:
            if 'MJD_OBS' in line:
                obsdate = float(line[11:25])
                t = Time(obsdate, format='mjd').jd
                dt[ii] = t-t0
        cols[ii] = 'black'
    else:
        print("couldn't find telescope")
        print(tel)
order = np.argsort(dt)
files_sorted = files[order]
dt_sorted = dt[order]
cols = cols[order]

nfiles = len(files_sorted)

fig,ax = plt.subplots(1,1,figsize=(8,5))

# Loop through the sorted files, and plot the spectra
for ii,f in enumerate(files_sorted):
    # In Dan's 18cow paper, he interpolates over host narrow features
    print(f)
    tel = f.split("_")[2]
    dat = np.loadtxt(f)
    wl = dat[:,0]
    flux = dat[:,1]

    # Plot the spectrum
    x = wl
    y = flux 
    ax.plot(x, y, c='k', drawstyle='steps-mid', lw=0.5)

    # Label the dt
    #t_raw = f.split("_")[1]
    #t = Time(t_raw[:4] + "-" + t_raw[4:6] + "-" + t_raw[6:8], format='iso').jd
    #dt = t-t0
    dt_str = r"$\Delta t$=%s d" %str(np.round(dt_sorted[ii], 1))
    ax.text(
            0.98, 0.9, s=dt_str, 
            horizontalalignment='right', verticalalignment='center', 
            fontsize=14, transform=ax.transAxes)
    ax.text(
            0.98, 0.7, s=tel, 
            horizontalalignment='right', verticalalignment='center', 
            fontsize=14, transform=ax.transAxes)
    ax.set_ylabel(
        r"Flux $f_{\lambda}$", fontsize=16)
    ax.yaxis.set_ticks([])
    if tel == 'NOT':
        ax.set_ylim(min(y)/15, max(y)/2) # get rid of noisy end

    ax.set_xlabel(
            r"Observed Wavelength (\AA)", fontsize=16)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.set_xlim(3000,10500)
    #ax.set_ylim(-25, 50)
    #ax.legend(loc='upper right', fontsize=12)
    #ax.set_xscale('log')

    plt.tight_layout()
    plt.savefig("spec_%s.png" %ii)
    plt.close()
#plt.show()
