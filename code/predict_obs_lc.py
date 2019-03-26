""" 
From the model properties, predict what the observed light curve should be
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import blackbody
from load_lc import get_lc
from load_model import load
from astropy.cosmology import Planck15

d = Planck15.luminosity_distance(z=0.03154).cgs.value

dt, filt, mag, emag = get_lc()
det = np.logical_and(mag<99, ~np.isnan(mag))
use_f = 'g'
choose = np.logical_and(det, filt == use_f)
order = np.argsort(dt[choose])

mod_dt, mod_lum, mod_rad, mod_temp = load()

wl = 4669.53

flux = mod_lum / (4 * np.pi * d**2)

# Solve for the constant in front of the blackbody
# that gives you the correct bolometric luminosity

# then you'll have erg/cm2/sec/AA

fnu = 3.34E4 * wl**2 * flamb
mab = -2.5 * np.log10(fnu/3631)

plt.errorbar(
        dt[choose][order], mag[choose][order], emag[choose][order], 
        fmt='.', c='k', label="Obs. g-band LC")
plt.plot(mod_dt, mab, c='k', lw=0.5, ls='--', label="Predicted g-band LC")
plt.tick_params(axis='both', labelsize=14)
plt.xlabel("Days", fontsize=14)
plt.ylabel("Apparent Mag (AB)", fontsize=14)
plt.gca().invert_yaxis()
plt.legend(loc='lower center', fontsize=14)
plt.tick_params(axis='both', labelsize=14)
plt.xlim(-1,10)

plt.show()
