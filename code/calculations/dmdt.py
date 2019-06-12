""" Calculate the dmag/dt from the first few min of the light curve """

import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/code")
from load_lc import get_forced_phot

jd, filt, mag, emag, limmag, code = get_forced_phot()
t0 = 2458370.6634
dt = (jd-t0)*24*60 # minutes
choose = np.logical_and.reduce(
        (code=='ZTF Camera', filt=='ztfg', ~np.isnan(mag), dt>0, dt<100))
plt.errorbar(
        dt[choose], mag[choose], yerr=emag[choose],
        c='k', ms=5, fmt='s', label="P48 $g$", zorder=2)
plt.gca().invert_yaxis()

# Do a Monte Carlo to get the error bar on the dmag/dt
nsim = 10000
x = dt[choose]
y = np.zeros((nsim,len(x)))

# Make 10000 new version of the mag array, repopulated
# from the Gaussian dist with STD of that uncertainty
for ii,mval in enumerate(mag[choose]):
    y[:,ii] = np.random.normal(loc=mval, scale=emag[choose][ii], size=nsim)

slopes = np.zeros(nsim)
bvals = np.zeros(nsim)

for ii in np.arange(nsim-1):
    m,b = np.polyfit(x, y[ii], deg=1)
    slopes[ii] = m
    bvals[ii] = b

plt.figure()
plt.errorbar(x, y[5], yerr=emag[choose], fmt='.')
xvals = np.linspace(min(x), max(x), 1000)
yvals = slopes[5]*xvals+bvals[5]
plt.plot(xvals, yvals)
print(slopes[5]*24)
plt.show()

mfit = np.mean(slopes)*60 # mag/hr
emfit = np.std(slopes)*60

print(mfit)
print(emfit)
