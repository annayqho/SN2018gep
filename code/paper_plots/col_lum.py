""" Color evolution, luminosity evolution phase space """

import numpy as np
import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)

# To start, just plot AT2018gep in this phase space.
DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/phot"

fig,ax = plt.subplots(1,1,figsize=(8,5))

f = DATA_DIR + "/ZTF18abukavn_opt_phot.dat"
dat = np.loadtxt(f, dtype=str, delimiter=' ')
instr = dat[:,0]
jd = dat[:,1].astype(float)
filt = dat[:,2]
mag = dat[:,3].astype(float)
emag = dat[:,4].astype(float)

det = np.logical_and(mag<99, ~np.isnan(mag))
nondet = np.logical_or(mag==99, np.isnan(mag))
zp = 2458370.6473
dt = jd-zp

band = filt=='g'
choose = np.logical_and(det, band)
order = np.argsort(dt[choose])
ax.plot(dt[choose][order], mag[choose][order], c='darkgreen')
ax.fill_between(dt[choose][order], mag[choose][order]-emag[choose][order],
        mag[choose][order]+emag[choose][order], color='lightgreen', alpha=0.3)

band = filt=='r'
choose = np.logical_and(det, band)
order = np.argsort(dt[choose])
ax.plot(dt[choose][order], mag[choose][order], c='red')
ax.fill_between(dt[choose][order], mag[choose][order]-emag[choose][order],
        mag[choose][order]+emag[choose][order], color='red', alpha=0.3)

ax.tick_params(axis='both', labelsize=14)
ax.invert_yaxis()

plt.tight_layout()

plt.show()
