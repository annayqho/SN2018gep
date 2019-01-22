""" Calculate the host redshift """

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sys
sys.path.append("/Users/annaho/Github/Fit_Redshift")
from fitlines import *

specfile = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/spec/ZTF18abukavn/ZTF18abukavn_20181109_Keck1_v1.ascii"
# Initial guess
z0 = 0.0322
# Window size, in angstroms
window = 20

# Fit for centroids of as many lines as you can tolerate
balmer = np.array([6564.61, 4862.68, 4341.68, 4102.89, 3970.072])
oiii = np.array([4363, 4932.6, 4960.295, 5008.24]) # O III
# Strong lines
lines = np.hstack((balmer[0], balmer[1], oiii[-1]))
zall = []
ezall = []

# Solve
for line in lines:
    z, ez = fit_redshift(specfile, line, z0, window)
    zall.append(z)
    ezall.append(ez)
zall = np.array(zall)
ezall = np.array(ezall)

# Use the STD of the best fits as the uncertainty
w = 1/ezall**2
zmean = np.average(zall, weights=w)
ezmean = np.std(zall)

# Print the best-fit redshift, and uncertainty
print("%s +/- %s" %(np.round(zmean,7), np.round(ezmean, 7)))
