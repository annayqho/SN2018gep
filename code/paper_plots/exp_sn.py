""" Exponential function with a supernova """

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import numpy as np
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/code")
from load_lum import load_lc

# Import the bolometric luminosity
dt, lum, llum, ulum = load_lc()

# Plot an exponential function with timescale = diffusion time

# Plot SN2010bh x2
