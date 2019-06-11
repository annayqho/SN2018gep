""" Calculate the dmag/dt from the first few min of the light curve """

import matplotlib.pyplot as plt
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

for ii in np.arange(nsim-1):
    m,b = np.polyfit(x, y[ii], deg=1)
    slopes[ii] = m

mfit = np.mean(slopes)*24*60 # mag/day
emfit = np.std(slopes)*24*60

