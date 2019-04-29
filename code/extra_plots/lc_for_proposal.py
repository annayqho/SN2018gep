import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import sys
from astropy.cosmology import Planck15
from astropy.table import Table
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF_Tools")
#from query_lc import run

distmod = 35.87616973979363

fig,ax = plt.subplots(1, 1, figsize=(10,3))

#t,mag,emag,fid = run('ZTF18abukavn')
dat = Table.read("../data/lc_marshal.txt", format='ascii')
jd = dat['jdobs']
filt = dat['filter']
M = dat['absmag']
inst = dat['instrument']
t0 = 2458370.6473

# Plot the P48 g and r band
choose = np.logical_and(filt=='g', inst=='P48+ZTF')
plt.scatter(jd[choose]-t0, M[choose], c='green', label='P48+ZTF g')
plt.plot(jd[choose]-t0, M[choose], c='green', label='_nolegend_')

choose = np.logical_and(filt=='r', inst=='P48+ZTF')
plt.scatter(jd[choose]-t0, M[choose], c='red', label='P48+ZTF r', marker='s')
plt.plot(jd[choose]-t0, M[choose], c='red', ls='--', label='_nolegend_')

# Plot Swift B-band
# choose = np.logical_and(filt=='B', inst=='Swift+UVOT')
# plt.scatter(
#         jd[choose]-t0, M[choose], facecolor='deepskyblue', 
#         edgecolor='k', label='Swift+UVOT B', marker='o')

# Plot Swift Plot Swift V
# choose = np.logical_and(filt=='V', inst=='Swift+UVOT')
# plt.scatter(
#         jd[choose]-t0, M[choose], facecolor='greenyellow', 
#         edgecolor='k', label='Swift+UVOT V', marker='s')

# Plot Swift UVM2
# choose = np.logical_and(filt=='UVM2', inst=='Swift+UVOT')
# plt.scatter(
#         jd[choose]-t0, M[choose], facecolor='darkblue', 
#         edgecolor='k', label='Swift+UVOT UVM2', marker='s')

# Plot Swift UVW1
# choose = np.logical_and(filt=='UVW1', inst=='Swift+UVOT')
# plt.scatter(
#         jd[choose]-t0, M[choose], facecolor='mediumpurple', 
#         edgecolor='k', label='Swift+UVOT UVW1', marker='o')

# Plot Swift UVW2
choose = np.logical_and(filt=='UVW2', inst=='Swift+UVOT')
plt.scatter(
        jd[choose][:-2]-t0, M[choose][:-2], facecolor='black', 
        edgecolor='k', label='Swift+UVOT UVW2', marker='v', s=60)
plt.plot(
        jd[choose][:-2]-t0, M[choose][:-2], c='k', label='_nolegend_')


# plt.scatter(
#         t[fid==2]-t[fid==1][0], mag[fid==2]-distmod, 
#         c='red', label='ZTF r')
plt.xlabel("Days since discovery", fontsize=16)
plt.ylabel("Abs Mag", fontsize=16)
plt.tick_params(axis='both', labelsize=14)
plt.xlim(-1,12)
#plt.xscale('log')
plt.ylim(-21, -15)
plt.gca().invert_yaxis()

plt.legend(ncol=2, loc='lower center')
plt.tight_layout()
plt.show()
#plt.savefig("lc.png")
