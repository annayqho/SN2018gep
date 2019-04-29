""" Reproduce this X-ray luminosity vs. isotropic gamma-ray energy plot from
Margutti (2014) """

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15

fig = plt.figure(figsize=(7,5))

datadir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data"
inputf = datadir + "/margutti2014.txt"
dat = np.loadtxt(inputf, dtype=str, delimiter=',')
names = dat[:,0]
x = dat[:,1].astype(float)
y = dat[:,2].astype(float)

# 98bw, 06aj, 03lw
use = ['98bw', '03lw', '06aj']
for name in use:
    choose = names == name
    if name=='98bw':
        plt.scatter(
                x[choose], y[choose], marker='o', c='k', s=50)
                #label="LLGRB-SNe")
    else:
        plt.scatter(
                x[choose], y[choose], marker='o', c='k', s=50)
    plt.text(
            x[choose]*1.01, y[choose]*1.01, name, fontsize=12,
            horizontalalignment='left', verticalalignment='bottom')

# 100316D, 120422A
use = ['12bz', '10bh']
for name in use:
    choose = names == name
    plt.scatter(x[choose], y[choose], marker='o', c='k', s=50)
    plt.arrow(x[choose][0], y[choose][0], 0, 3E49, length_includes_head=True,
            head_width=x[choose][0]/7, head_length=y[choose][0]/5, 
            fc='k', ec='k')
    plt.text(x[choose], y[choose]/1.4, name, fontsize=12,
            horizontalalignment='center', verticalalignment='top')

# 09bb
use = ['09bb']
for name in use:
    choose = names == name
    plt.scatter(
            x[choose], y[choose], marker='o', 
            edgecolor='k', facecolor='black', s=50)#, label="Ic-BL SN")
    plt.arrow(x[choose][0], y[choose][0], 0, -y[choose][0]/2, 
            length_includes_head=True,
            head_width=x[choose][0]/10, head_length=y[choose][0]/6, 
            fc='k', ec='k')
    plt.text(x[choose], y[choose]*1.1, name, fontsize=12,
            horizontalalignment='center', verticalalignment='bottom')

# 2012ap
use = ['12ap']
for name in use:
    choose = names == name
    plt.scatter(
            x[choose], y[choose], marker='o', 
            edgecolor='k', facecolor='k', s=50)
    plt.arrow(x[choose][0], y[choose][0], 0, -y[choose][0]/2, 
            length_includes_head=True,
            head_width=x[choose][0]/5, head_length=y[choose][0]/6, 
            fc='k', ec='k')
    plt.arrow(x[choose][0], y[choose][0], -x[choose][0]/2, 0, 
            length_includes_head=True,
            head_width=y[choose][0]/5, head_length=x[choose][0]/10, 
            fc='k', ec='k')
    plt.text(x[choose], y[choose]*1.1, name, fontsize=12,
            horizontalalignment='center', verticalalignment='bottom')

# AT2018gep
x = 5.391E39
y = 9.337E46
name = '18gep'
plt.scatter(
        x, y, marker='*', c='k', facecolor='white', s=100)
plt.arrow(x, y, 0, -y/2, 
        length_includes_head=True,
        head_width=x/5, head_length=y/6, 
        fc='k', ec='k')
plt.arrow(x, y, -x/4, 0, 
        length_includes_head=True,
        head_width=y/5, head_length=x/20, 
        fc='k', ec='k')
plt.text(x*1.01, y*1.01, name, fontsize=12,
        horizontalalignment='left', verticalalignment='bottom')

# iPTF16asu
x = 1.1E43
y = 8.2E50
name = '16asu'
plt.scatter(
        x, y, marker='o', edgecolor='k', facecolor='k', s=50)
plt.arrow(x, y, 0, -y/2, 
        length_includes_head=True,
        head_width=x/8, head_length=y/6, 
        fc='k', ec='k')
plt.arrow(x, y, -x/2, 0, 
        length_includes_head=True,
        head_width=y/5, head_length=x/8, 
        fc='k', ec='k')
plt.text(x*1.01, y*1.01, name, fontsize=12,
        horizontalalignment='left', verticalalignment='bottom')


plt.xscale('log')
plt.yscale('log')
plt.tick_params(labelsize=14)
plt.xlim(1E39, 1E44)
plt.ylim(3E46, 3E51)
plt.xlabel(r"$L_X$ [0.3-30 keV] (erg s$^{-1}$)", fontsize=16)
plt.ylabel(r"$E_\mathrm{\gamma,iso}$ [1-$10^4$ keV] (erg)", fontsize=16)
plt.tight_layout()
#plt.legend(fontsize=14, loc='upper left')

#plt.show()
plt.savefig("margutti2014.png")
