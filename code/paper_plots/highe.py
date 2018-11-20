""" High-energy limit vs GRB fluence for various other things """

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import numpy as np
from astropy.cosmology import Planck15


# I think that what you want to do is convert to an isotropic peak luminosity?
# isotropic energy?

# iPTF16asu
# the limit is 1E-7 erg/cm2/s
# 10 keV - 10 MeV energy range
# associated isotropic peak luminosity limit is Liso <= 1E49 erg/s
# associated total energy is Eiso <= 1E50 erg 

# GRB 980425: Liso ~ 5E46 erg/s, Eiso ~ 1E48 erg... Galama 1998b

# AT2018gep
d = Planck15.luminosity_distance(z=0.033).cgs.value
fluence = 6E-9 # erg/cm2
lum = fluence * 4 * np.pi* d**2
print(lum)

# the associated 
