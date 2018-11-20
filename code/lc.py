""" Plot the light curve compared to other rapid transients """

import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.cosmology import Planck15
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF_Tools")
from query_lc import run

# ZTF18abukavn
t,mag,emag,fid = run('ZTF18abukavn')
distmod = 35.87616973979363

# g-band
plt.scatter(
        t[fid==1]-t[fid==1][0]-3, mag[fid==1]-distmod, c='green',
        label='ZTF18abukavn')

# iPTF16asu
# they estimate an explosion time of 2016 May 10.53 ± 0.17 days
distmod = 39.86226641498726
asu = np.loadtxt(
        "lc_16asu.txt", delimiter='&', dtype=str)
#t_asu_raw = asu[:,1].astype(float) # dt in days, phase
dt_raw = asu[:,1].astype(float) # dt in days, phase
#dt_raw = t_asu_raw - t_asu_raw[0]
filter_asu = np.array([val.strip() for val in asu[:,2]])
mag_asu_raw = np.array([str(val.strip()) for val in asu[:,3]]).astype(str)
islim = np.array(['>' in val for val in mag_asu_raw])

choose = np.logical_and(filter_asu == 'g', ~islim)
dt_asu = dt_raw[choose]
mag_asu = np.array(
        [float(val.split('$pm$')[0]) for val in mag_asu_raw[choose]])
emag_asu = np.array(
        [float(val.split('$pm$')[1]) for val in mag_asu_raw[choose]])
plt.scatter(dt_asu, mag_asu-distmod, c='grey', s=5)
plt.fill_between(
        dt_asu, mag_asu-distmod-emag_asu/2, mag_asu-distmod+emag_asu/2,
        color='grey', lw=2.0, alpha=0.5, label='iPTF16asu')

# plot the last upper limit
choose = np.logical_and(filter_asu == 'g', islim)
dt_asu = dt_raw[choose]
mag_asu = np.array(
        [float(val[2:7]) for val in mag_asu_raw[choose]])
plt.scatter(dt_asu, mag_asu-distmod, c='grey', marker='v')

plt.gca().invert_yaxis()
plt.xlim(-5,8)
plt.xlabel("$\Delta t$", fontsize=14)
plt.ylabel("Mag", fontsize=14)
plt.legend(loc='lower right')

plt.show()
#plt.savefig("lc_16asu_comparison.png")

