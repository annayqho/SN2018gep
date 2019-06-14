""" Check whether you need to do a K-correction """

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import trapz
from astropy.modeling.blackbody import blackbody_lambda
from astropy.io import ascii

z = 0.03154
dat = ascii.read("Swift_UVOT.UVW2.dat") # not really necessary

# UVW2: 
leff = 2085.7
weff = 667.7

wl = np.linspace(1000,10000,1000)
f = blackbody_lambda(wl, 50000)
plt.plot(wl, f, c='k', label="BB at T=50,000K")

# rest-frame
xi = leff-weff/2
xf = leff+weff/2
plt.axvline(x=xi, c='red', label="UVW2 at z=0")
plt.axvline(x=xf, c='red')

# redshifted
plt.axvline(x=xf*(1+z), c='r', ls='--', label="UVW2 at z=0.03")
plt.axvline(x=xi*(1+z), c='r', ls='--')

plt.legend(loc='upper right', fontsize=14)
plt.xlim(1000,5000)
plt.xlabel("Wavelength (AA)", fontsize=14)
plt.ylabel("Flux (erg/AA/cm2/s/sr)", fontsize=14)
plt.gca().tick_params(labelsize=12)
plt.tight_layout()
plt.show()

# Do the integral with trapz
precorr = trapz(f,x=[xi,xf])/1E15
postcor = trapz(f,x=[xi*(1+z),xf*(1+z)])/1E15
dmag = -2.5*np.log10(precorr/postcor)
print(dmag)
