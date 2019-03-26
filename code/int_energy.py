""" Calculate the integrated energy output """

import numpy as np
from load_lum import load_lc

dt,lum,llum,ulum = load_lc()
print(np.trapz(dt*86400, lum))
