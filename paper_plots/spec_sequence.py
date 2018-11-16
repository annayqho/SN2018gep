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

files = glob.glob("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/spec/ZTF18abukavn/*.ascii")

t0 = 2458370.6634 # in JD
for ii,f in enumerate(files[0:5]):
    # In Dan's 18cow paper, he interpolates over host narrow features
    dat = np.loadtxt(f)
    wl = dat[:,0]
    flux_raw = dat[:,1]

    # If the spectra are from LT, then you multiply by 1E-15?
    tel = f.split("_")[2]
    if tel == 'LT':
        flux = flux_raw*1E-15
        # and read in the observation time
        alldat = open(f).readlines()
        for line in alldat:
            if 'DATE-OBS' in line:
                obsdate = line[13:36]
                print(obsdate)
                t = Time(obsdate, format='isot').jd
                dt = t-t0
    elif tel == 'P200':
        flux = flux_raw
        for line in alldat:
            if 'DATE-OBS' in line:
                print(line)
                obsdate = line[13:36]
                print(obsdate)
                t = Time(obsdate, format='isot').jd
                dt = t-t0
    elif tel == 'Keck1':
        flux = flux_raw
    elif tel == 'DCT':
        flux = flux_raw
    elif tel == 'NOT':
        flux = flux_raw
    elif tel == 'P60':
        flux = flux_raw

    # Plot the spectrum
    ax.plot(wl, flux-ii*1E-15, c='k', alpha=0.7, drawstyle='steps-mid')
    
    # Label the dt
    #t_raw = f.split("_")[1]
    #t = Time(t_raw[:4] + "-" + t_raw[4:6] + "-" + t_raw[6:8], format='iso').jd
    #dt = t-t0
    dt_str = r"$\Delta t$=%s d" %str(np.round(dt, 1))
    ax.text(
            wl[-1], flux[-1]-ii*1E-15, s=dt_str, 
            horizontalalignment='left', verticalalignment='center', 
            fontsize=14)

ax.set_ylabel(
    r"Flux $f_{\lambda}$ (arbitrary units)", fontsize=16)

ax.set_xlabel(
        r"Observed Wavelength (\AA)", fontsize=16)
ax.yaxis.set_ticks([])
ax.xaxis.set_tick_params(labelsize=14)
ax.set_xlim(3700,12000)
ax.set_ylim(-5E-15,3E-15)
#ax.legend(loc='upper right', fontsize=12)
#ax.set_xscale('log')

#plt.tight_layout()
#plt.savefig("lc.png")
plt.show()
