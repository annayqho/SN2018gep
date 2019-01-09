""" Calculate the integrated energy output """

import numpy as np
from load_lum import load_lc

dt,lum,llum,ulum = load_lc()
choose = dt <= 8
print(np.trapz(dt[choose]*86400, lum[choose]))
