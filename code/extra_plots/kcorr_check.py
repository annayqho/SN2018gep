""" Check whether you need to do a K-correction """

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import trapz
from astropy.modeling.blackbody import blackbody_lambda

wl = np.linspace(1000,10000,1000)
f = blackbody_lambda(wl, 50000)
plt.plot(wl, f, c='k', label="BB at T=50,000K")
plt.axvline(x=2588,c='r',label="At z=0.03")
plt.axvline(x=1928+585, label="UVW2 at z=0")
plt.legend(loc='upper right', fontsize=14)
plt.xlabel("Wavelength (AA)", fontsize=14)
plt.ylabel("Flux (erg/AA/cm2/s/sr)", fontsize=14)
plt.tick_params(fontsize=12)
plt.tick_params(labelsize=12)
plt.tight_layout()
plt.show()

# Do the integral with trapz
# trapz(f,x=[1383,2588])
# trapz(f,x=[1343,2513])
# -2.5*np.log10(5.51E15/5.35E15)
