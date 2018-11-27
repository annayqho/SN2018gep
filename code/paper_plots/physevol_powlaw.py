""" Plot the physical evolution as a power law """

import numpy as np
import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/code")
from fit_arnett import load_lc
from extrapolate_size import load_radius
from load_temp import load_temp

# Load the bolometric light curve
dt, lum, elum = load_lc()

# Load the radius
dt, rad, erad = load_radius()

# Load the temperature
dt, temp, etemp = load_temp()

# Initialize the figure
fig,axarr = plt.subplots(3,1, figsize=(6,8), sharex=True)

# Plot the data

# Luminosity panel
axarr[0].errorbar(dt, lum, yerr=elum, fmt='.', c='k')
choose = dt < 2.5
m,b = np.polyfit(np.log10(dt[choose]), np.log10(lum[choose]), deg=1)
xfit = np.linspace(min(dt[choose]), max(dt[choose]))
yfit = 10**(m*np.log10(xfit)+b)
axarr[0].plot(xfit, yfit, ls='--', c='grey')
axarr[0].text(1E-1, 1E44, '$t^{-0.17}$',
        horizontalalignment='left', verticalalignment='center', fontsize=14)

choose = dt > 2.5
m,b = np.polyfit(np.log10(dt[choose]), np.log10(lum[choose]), deg=1)
m = -5/3
xfit = np.linspace(min(dt[choose]), max(dt[choose]))
yfit = 10**(m*np.log10(xfit)+b)
axarr[0].plot(xfit, yfit, ls='--', c='grey')
axarr[0].text(10, 2E43, '$t^{-5/3}$',
        horizontalalignment='left', verticalalignment='center', fontsize=14)

# Radius panel

# Radius panel
axarr[1].errorbar(dt, rad, yerr=erad, fmt='.', c='k')

# Plot lines of constant velocity
xvals = np.linspace(1E-3, 1E2, 1000)
# v = c
yvals = (3E10) * xvals * 86400
axarr[1].plot(xvals, yvals, ls='--', c='grey')
axarr[1].text(1E-1, 8E14, 'v=c', fontsize=12)

# v = 0.1c
yvals = 0.1 * (3E10) * xvals * 86400
axarr[1].plot(xvals, yvals, ls='--', c='grey')
axarr[1].text(1, 2E14, 'v=0.1c', fontsize=12)

# Temperature panel
axarr[2].errorbar(dt, temp, yerr=etemp, fmt='.', c='k')
choose = dt < 2.5
m,b = np.polyfit(np.log10(dt[choose]), np.log10(temp[choose]), deg=1)
xfit = np.linspace(min(dt[choose]), max(dt[choose]))
yfit = 10**(m*np.log10(xfit)+b)
axarr[2].plot(xfit, yfit, ls='--', c='grey')
axarr[2].text(1E-1, 2E4, '$t^{-0.11}$',
        horizontalalignment='left', verticalalignment='center', fontsize=14)

choose = np.logical_and(dt > 2.5, dt < 11)
m,b = np.polyfit(np.log10(dt[choose]), np.log10(temp[choose]), deg=1)
print(m)
xfit = np.linspace(min(dt[choose]), max(dt[choose]))
yfit = 10**(m*np.log10(xfit)+b)
axarr[2].plot(xfit, yfit, ls='--', c='grey')
axarr[2].text(10, 1E4, '$t^{-0.89}$',
        horizontalalignment='left', verticalalignment='center', fontsize=14)

# Formatting
axarr[0].xaxis.label.set_visible(False)
axarr[1].xaxis.label.set_visible(False)

axarr[2].tick_params(axis='both', labelsize=12)
axarr[0].tick_params(axis='y', labelsize=12)
axarr[1].tick_params(axis='y', labelsize=12)

for ii in np.arange(0,3):
    axarr[ii].set_yscale('log')

axarr[1].set_ylim(1E14, 1E16)
axarr[2].set_xlim(2E-3, 50)
axarr[2].set_xscale('log')
axarr[2].set_xlabel(r'Days since $t_0$', fontsize=16)
axarr[0].set_ylabel(r'$L_\mathrm{bol}$ (erg/s)', fontsize=16)
axarr[1].set_ylabel(r'$R_\mathrm{ph}$ (cm)', fontsize=16)
axarr[2].set_ylabel(r'$T_\mathrm{eff}$ (K)', fontsize=16)


plt.subplots_adjust(hspace=0)
plt.tight_layout()
#plt.show()
plt.savefig("bbfit_log.png")
