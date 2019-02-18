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
from plot_lc import get_lc

files = glob.glob(
"/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/spec/ZTF18abukavn/*.ascii")
files = np.array(files[0:10])
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

colors = np.array(['lightpink', 'pink', 'r', 
    'orange', 'darkorange', 'yellow', 'lightgreen', 'green', 
    'darkgreen', 'lightblue', 'blue', 'darkblue', 'purple', 
    'black'])
ncolors = len(colors)

nfiles = len(files_sorted)

fig,ax = plt.subplots(1,1,figsize=(8,10))

# Loop through the sorted files, and plot the spectra
for ii,f in enumerate(files_sorted):
    print(f)
    tel = f.split("_")[2]
    dat = np.loadtxt(f)
    wl = dat[:,0]
    flux = dat[:,1]
    dt_spec = dt_sorted[ii]

    # FLUX CAL TO R-BAND LIGHT CURVE    
    # get r-band LC
    dt, filt, det, mag, emag = get_lc()
    choose = np.logical_and(det, filt=='r')
    # interpolate to this epoch
    rval = np.interp(dt_spec, dt[choose], mag[choose])
    # TEMP: r is roughly 658nm +/- 138nm
    # TEMP: assume AB mag
    lam = 6580 # in angstroms
    c = 3E18 # angstrom/s
    fnu = 1E-23 * 3631 * 10**(rval/(-2.5)) # erg/s/cm2/Hz
    flam = fnu * (c/lam**2) # should be erg/s/cm2/AA
    print("The flux at this wavelength is:")
    print(flam)
    # scale factor
    flam_meas = np.interp(lam, wl, flux)
    scale = (flam/flam_meas)/1E-15

    # Plot the spectrum, scaled to this flux value
    x = wl
    y = flux 
    ax.plot(
            x, y*scale+(nfiles-ii), drawstyle='steps-mid', 
            lw=0.5, c=colors[ii%ncolors])

    # Label the dt
    #t_raw = f.split("_")[1]
    #t = Time(t_raw[:4] + "-" + t_raw[4:6] + "-" + t_raw[6:8], format='iso').jd
    #dt = t-t0
    dt_str = r"$\Delta t$=%s d" %str(np.round(dt_spec, 1))
    ax.text(
            8000, (y*scale+(nfiles-ii))[-1], s=dt_str, 
            horizontalalignment='right', verticalalignment='bottom', 
            fontsize=14)#, transform=ax.transAxes)
    # ax.text(
    #         0.98, 0.7, s=tel, 
    #         horizontalalignment='right', verticalalignment='center', 
    #         fontsize=14, transform=ax.transAxes)
#ax.set_ylabel(
#    r"$f_{\lambda}$ ($10^{-15}$ erg cm$^{-2}$ s$^{-1}$ \AA$^{-1}$)", 
#    fontsize=16)
ax.set_ylabel("$f_{\lambda}$ + offset", fontsize=16)
#ax.yaxis.set_ticks([])
#if tel == 'NOT':
#    ax.set_ylim(min(y)/15, max(y)/2) # get rid of noisy end

ax.set_xlabel(
        r"Observed Wavelength (\AA)", fontsize=16)
ax.xaxis.set_tick_params(labelsize=14)
ax.yaxis.set_tick_params(labelsize=14)
ax.set_xlim(3800,8100)
ax.set_ylim(1, 13)
#ax.legend(loc='upper right', fontsize=12)
#ax.set_xscale('log')

plt.tight_layout()
#plt.savefig("spec_%s.png" %ii)
#plt.show()
#plt.close()
plt.show()
