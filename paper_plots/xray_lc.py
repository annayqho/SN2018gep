import glob
import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15

def plot_lc(f, name=None):
    dt = []
    lum = []
    with open(f, "r") as inputf:
        for line in inputf.readlines():
            if len(line) > 40:
                dt_s = float(line.split('\t')[0])
                dt.append(dt_s/86400)
                flux = float(line.split('\t')[3])
                lum.append(flux * 4 * np.pi * d**2)
    if name:
        plt.plot(dt, lum, c='grey', alpha=1.0, lw=2, zorder=1)
    else:
        # just thin grey
        plt.plot(dt, lum, c='grey', alpha=0.1, lw=1)
    return dt, lum


# ZTF18abukavn
d = Planck15.luminosity_distance(z=0.033).cgs.value
ratio = 3.27E-11 # count-to-flux rate, erg/cm2/ct

xrtlc = Table.read(
    "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/from_brad.dat",
    format='ascii')
# Brad's zero-point is: MJD=58370.588137 = 2018 Sep 09 at 14:06:55.064 UT.
# Our zero-point is 2458370.6473
# So the difference is 58370.6473-58370.588137 = 0.059163
dt_sec = xrtlc['col1']+0.059163
ct = xrtlc['col4']
lum = ct * ratio * 4 * np.pi * d**2 
dt_day = dt_sec/86400
plt.scatter(dt_day, lum, marker='v', c='k', zorder=5)
plt.text(dt_day[0], lum[0]*1.2, 'AT2018gep', fontsize=14,
        horizontalalignment='center', verticalalignment='bottom')
# Chandra
plt.scatter(16, 3E-15*4*np.pi*d**2, marker='v', c='k', zorder=5)


# All the GRBs
data_dir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/xrtlc"
# flist = glob.glob(data_dir + "/*_xrt_bin.txt")
# for f in flist:
#     plot_lc(f)

# Individual LLGRBs
f = data_dir + "/060218_xrt_bin.txt" 
dt, lum = plot_lc(f, name="060218")
plt.text(
    0.06, 1.26E46, "060218", 
    horizontalalignment='center',
    verticalalignment='bottom', fontsize=12)

f = data_dir + "/100316d_xrt_bin.txt" 
dt, lum = plot_lc(f, name="100316D")
plt.text(
    dt[0], lum[0]/2, '100316D', 
    horizontalalignment='left',
    verticalalignment='top', fontsize=12)

f = data_dir + "/030329_xray.dat"
# this one has a different file
dat = np.loadtxt(f)
dt = dat[:,0] 
flux = dat[:,1] * 1E-12
lum = flux * 4 * np.pi * d**2
plt.plot(dt, lum, c='grey', lw=2, zorder=1)
plt.text(
    dt[0], lum[0], '030329',
    horizontalalignment='center',
    verticalalignment='bottom', fontsize=12)

# Need 2008D
f = data_dir + "/2008D.txt"
dat = np.loadtxt(f, dtype=str, delimiter='&', skiprows=1)
dt_sec = np.array(
        [string.replace('\xa0', '') for string in dat[:,0]]).astype(float)
dt_day = dt_sec/86400
lum = 1E42*np.array(
        [string.replace('\xa0', '') for string in dat[:,2]]).astype(float)
elum = 1E42*np.array(
        [string.replace('\xa0', '') for string in dat[:,3]]).astype(float)
plt.plot(dt_day, lum, c='grey', lw=2, zorder=1)
plt.text(
    3E-3, 2E43, '2008D',
    horizontalalignment='left',
    verticalalignment='center', fontsize=12)

# 2009bb
# Soderberg 2009
# detection from Chandra
d_09bb = 40*3.086E24
plt.scatter(31, 4.4E39, marker='o', c='grey', s=20)
plt.text(31, 4.4E39*2, "2009bb", fontsize=12,
        horizontalalignment='center', verticalalignment='center')
# upper limits from Swift/XRT
# dt = np.array([5, 19, 23, 31])
# flim = np.array([1.3E-13, 1.7E-13, 2.5E-13, 3E-13])
# llim = 4 * np.pi * d_09bb**2 * flim
# plt.scatter(dt, llim, marker='v', c='grey', s=20)
# plt.plot(dt, llim, c='grey', lw=2)

# 1998bw
# The original four points of 980425 reported by Pian et al. (2000)
# Fig 7b
d_98bw = 38*3.086E24
dt = np.array([1, 2, 7.5, 200])
f = np.array([4.3E-13, 4.2E-13, 2.8E-13, 1.7E-13])
lum = 4 * np.pi * d_98bw**2 * f
# Kouvelioutou 2004, Chandra: Day 1281, 1.2E39
# Sixth point is XMM measurement reported by Pian et al. (2004) from 2002 March 28,
# but XMM can't resolve two sources, so this luminosity value includes two sources
# They re-analyzed this data...
# dt = 1281
# lum = 1.2E39
plt.plot(dt, lum, c='grey', lw=2, zorder=1)
plt.text(
    dt[0], lum[0], '1998bw',
    horizontalalignment='right',
    verticalalignment='center', fontsize=12)

# 16asu limits
dt = np.array([7.4, 13.4, 19.2])
llim = np.array([2.5E43, 1.1E43, 1.5E43])
plt.scatter(dt, llim, marker='v', c='grey', s=40, facecolor='white')
plt.text(dt[0]*1.2, llim[0], '16asu', horizontalalignment='left', fontsize=12)

# One point for iPTF17cw
# Chandra: 8 Feb 2017
plt.scatter(41.8, 1E41, marker='o', c='grey', s=20)
plt.text(41.8, 1E41, '17cw', fontsize=12, horizontalalignment='right')

# Limits for 2012ap
plt.scatter(24, 2.4E39, marker='v', c='grey', s=40, facecolor='white')
plt.text(
        24/1.1, 2.4E39, '2012ap', verticalalignment='center', 
        horizontalalignment='right', fontsize=12)

# Formatting
plt.xlim(1E-3,1E2)
plt.xlabel(r"$\Delta t$ [days]", fontsize=14)
plt.ylabel(
r"X-ray luminosity [0.3-10 keV, erg\,s${}^{-1}$]", fontsize=14)
plt.yscale('log')
plt.xscale('log')
plt.ylim(1E39, 1E47)
plt.xlim(1E-3, 6E1)
plt.tick_params(axis='both', labelsize=14)


plt.tight_layout()
#plt.show()
plt.savefig("xray_lc.png")
