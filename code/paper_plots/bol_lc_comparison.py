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


# Plot the bolometric light curve of ZTF18abukavn
dt, lum, llum, ulum = load_lc()

fig,ax = plt.subplots(1,1, figsize=(8,6), sharex=True)
ax.errorbar(dt, lum, yerr=[llum,ulum], fmt='.', c='k')
ax.text(dt[0], lum[0], 'AT2018gep', fontsize=14)

ax.tick_params(axis='both', labelsize=12)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylabel(r'$L_\mathrm{bol}$ (erg/s)', fontsize=16)
ax.set_xlabel(r'Days since $t_0$', fontsize=16)

plt.show()
