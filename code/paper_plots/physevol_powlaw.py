""" Plot the physical evolution as a power law """

import numpy as np
import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/code")
from load_lum import load_lc
from load_radius import load_radius
from load_temp import load_temp

# Load the bolometric light curve
dt, lum, llum, ulum = load_lc()

# Load the radius
dt, rad, lrad, urad = load_radius()

# Load the temperature
dt, temp, ltemp, utemp = load_temp()

# Initialize the figure
fig,axarr = plt.subplots(3,1, figsize=(6,8), sharex=True)

# Plot the data

# Luminosity panel
axarr[0].errorbar(dt, lum, yerr=[llum,ulum], fmt='.', c='k')
m,b = np.polyfit(np.log10(dt), np.log10(lum), deg=1)
m = -5/3
xfit = np.linspace(min(dt), max(dt))
yfit = 2.1*10**(m*np.log10(xfit)+b)  # chi by eye
axarr[0].plot(xfit, yfit, ls='--', c='grey')
axarr[0].text(7, 1E44, '$t^{-5/3}$',
        horizontalalignment='left', verticalalignment='center', fontsize=14)

# Radius panel

# Radius panel
axarr[1].errorbar(dt, rad, yerr=[lrad,urad], fmt='.', c='k')

# Plot lines of constant velocity
xvals = np.linspace(1E-3, 1E2, 1000)

# v = 0.1c
yvals = 0.1 * (3E10) * xvals * 86400
axarr[1].plot(xvals, yvals, ls='--', c='grey')
axarr[1].text(1, 2E14, 'v=0.1c', fontsize=12)

# v = 0.26c
yvals = 0.26 * (3E10) * xvals * 86400
axarr[1].plot(xvals, yvals, ls='--', c='grey')
axarr[1].text(0.4, 7E14, 'v=0.26c', fontsize=12)

# Temperature panel
choose = np.logical_and(dt>1, dt<15)
axarr[2].errorbar(dt, temp, yerr=[ltemp,utemp], fmt='.', c='k')
m,b = np.polyfit(np.log10(dt[choose]), np.log10(temp[choose]), deg=1)
xfit = np.linspace(min(dt), max(dt))
yfit = 10**(m*np.log10(xfit)+b)
mstr = str(np.round(m, 2))
axarr[2].plot(xfit, yfit, ls='--', c='grey')
axarr[2].text(10, 1E4, '$t^{%s}$' %mstr,
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
axarr[2].set_xlim(2E-1, 50)
axarr[2].set_xscale('log')
axarr[2].set_xlabel(r'Days since $t_0$', fontsize=16)
axarr[0].set_ylabel(r'$L_\mathrm{bol}$ (erg/s)', fontsize=16)
axarr[1].set_ylabel(r'$R_\mathrm{ph}$ (cm)', fontsize=16)
axarr[2].set_ylabel(r'$T_\mathrm{eff}$ (K)', fontsize=16)

plt.subplots_adjust(hspace=0)
plt.tight_layout()
#plt.show()
plt.savefig("bbfit_log.png")
