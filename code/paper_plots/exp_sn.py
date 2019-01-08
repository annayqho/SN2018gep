""" Exponential function with a supernova """

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import numpy as np
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/code")
from astropy.table import Table
from load_lum import load_lc


# Initialize the figure
fig,ax = plt.subplots(1,1,figsize=(6,4))

# Import the bolometric luminosity
dt, lum, llum, ulum = load_lc()
ax.errorbar(
        dt, lum, yerr=[llum, ulum], fmt='s', c='k',
        label="ZTF18abukavn (SN2018gep)")

# Plot an exponential function with timescale = diffusion time
A = lum[1]
B = dt[1]
tau = 2
xexp = np.linspace(dt[0], dt[-1], 1000)
yexp = A*np.exp(-(xexp-B)/tau) 
ax.plot(xexp, yexp, c='k', ls='--', lw=1, label="Exp, $\\tau=%s$d" %tau)

# Plot SN2010bh x2
ddir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/bol_lc"
dat = Table.read(ddir + "/sn2010bh.dat", format='ascii.fast_no_header')
dt = dat['col1']
lum = dat['col2']
ysn = 2*lum
ax.plot(dt, ysn, c='k', ls=':', lw=1, label='SN2010bh, x2')

# Plot the sum of the two
tot = ysn + A*np.exp(-(dt-B)/tau)
ax.plot(dt, tot, c='k', ls='-', lw=1, label="Exp + SN")

ax.set_yscale('log')
ax.set_ylim(1E42, 1E45)
ax.legend(fontsize=14)
ax.set_xlabel("$\Delta t$ (days)", fontsize=16)
ax.set_ylabel("Bolometric Luminosity", fontsize=16)
ax.tick_params(axis='both', labelsize=14)

plt.tight_layout()
#plt.show()
plt.savefig("exp_sn.png")
