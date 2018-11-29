""" Comparison of bolometric light curves """

import numpy as np
import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import sys
from astropy.table import Table
from astropy.cosmology import Planck15
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/code")
from load_lum import load_lc

# Where the bolometric light curve compilation lives
ddir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/bol_lc"

# Initialize figure
fig,ax = plt.subplots(1,1, figsize=(8,6), sharex=True)

# Plot the bolometric light curve of ZTF18abukavn
dt, lum, llum, ulum = load_lc()
ax.errorbar(dt, lum, yerr=[llum,ulum], fmt='.', c='k')
ax.text(dt[0], lum[0], 'AT2018gep', fontsize=14)

# Next: SN2008D
dat = Table.read(
        ddir + "/sn2008d.dat", delimiter='&', format='ascii.fast_no_header')
dt = dat['col1']
TBB_raw = dat['col2']
TBB = np.array([val.split('+')[0] for val in TBB_raw])
uTBB = np.array([val.split('+')[1].split('_')[0] for val in TBB_raw])
lTBB = np.array([val.split('+')[1].split('_')[1] for val in TBB_raw])
lum_raw = dat['col4']
lum = 10**(np.array(
    [val.split('+')[0] for val in lum_raw]).astype(float))
ulum = 10**(np.array(
    [val.split('+')[1].split('_')[0] for val in lum_raw]).astype(float))
llum = 10**(np.array(
    [val.split('+')[1].split('_')[1] for val in lum_raw]).astype(float))
ax.errorbar(dt, lum, yerr=[llum,ulum], fmt='.', c='grey')
ax.plot(dt, lum, c='grey')
#ax.fill_between(dt, y1=lum-ulum, y2=lum+llum, color='grey')
ax.text(dt[0], lum[0], 'SN2008D', fontsize=14)

# Formatting
ax.tick_params(axis='both', labelsize=12)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylabel(r'$L_\mathrm{bol}$ (erg/s)', fontsize=16)
ax.set_xlabel(r'Days since $t_0$', fontsize=16)

plt.show()
