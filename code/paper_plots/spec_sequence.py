""" Plot the spectral sequence """

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

# We have to apply offsets, because the spectra are too close to each other

fig,ax = plt.subplots(1,1,figsize=(8,10))

files = glob.glob(
"/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/spec/ZTF18abukavn/*.ascii")
files = np.array(files)
dt = np.zeros(len(files))

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
    elif tel == 'P200':
        for line in alldat:
            if 'UT shutter open' in line:
                obsdate = line[12:35]
                print(obsdate)
                t = Time(obsdate, format='isot').jd
                dt[ii] = t-t0
    elif tel == 'Keck1':
        for line in alldat:
            if 'DATE_BEG' in line:
                obsdate = line[13:32]
                t = Time(obsdate, format='isot').jd
                dt[ii] = t-t0
    elif tel == 'DCT':
        obsdate = '2018-09-14T00:00:00' # temporary
        t = Time(obsdate, format='isot').jd
        dt[ii] = t-t0
    elif tel == 'NOT':
        obsdate = '2018-09-17T00:00:00' # temporary
        t = Time(obsdate, format='isot').jd
        dt[ii] = t-t0
    elif tel == 'P60':
        for line in alldat:
            if 'MJD_OBS' in line:
                obsdate = float(line[11:25])
                t = Time(obsdate, format='mjd').jd
                dt[ii] = t-t0
    else:
        print("couldn't find telescope")
        print(tel)
order = np.argsort(dt)
files_sorted = files[order]
dt_sorted = dt[order]

# Loop through the sorted files, and plot the spectra
for ii,f in enumerate(files_sorted):
    # In Dan's 18cow paper, he interpolates over host narrow features
    tel = f.split("_")[2]
    dat = np.loadtxt(f)
    wl = dat[:,0]
    flux = dat[:,1]

    # Plot the spectrum
    x = wl
    y = flux / flux[-1] - ii*3
    ax.plot(x, y, c='k', alpha=0.7, drawstyle='steps-mid', lw=0.5)
    
    # Label the dt
    #t_raw = f.split("_")[1]
    #t = Time(t_raw[:4] + "-" + t_raw[4:6] + "-" + t_raw[6:8], format='iso').jd
    #dt = t-t0
    dt_str = r"$\Delta t$=%s d" %str(np.round(dt_sorted[ii], 1))
    ax.text(
            x[-1], y[-1], s=dt_str, 
            horizontalalignment='left', verticalalignment='center', 
            fontsize=14)

ax.set_ylabel(
    r"Flux $f_{\lambda}$ (arbitrary units)", fontsize=16)

ax.set_xlabel(
        r"Observed Wavelength (\AA)", fontsize=16)
ax.yaxis.set_ticks([])
ax.xaxis.set_tick_params(labelsize=14)
ax.set_xlim(3000,12000)
ax.set_ylim(-60, 20)
#ax.legend(loc='upper right', fontsize=12)
#ax.set_xscale('log')

#plt.tight_layout()
#plt.savefig("lc.png")
plt.show()
