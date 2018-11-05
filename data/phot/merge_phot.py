""" Merge the photometry.

Use: P48 photometry as is, P60 photometry as is,
LT photometry from Christoffer's individual files.

Columns should be:
    JD, filter, mag, mag_error (statistical)
"""

import numpy as np
from astropy.io import ascii

# get the P48 and P60 photometry
dat = np.loadtxt("full_marshal_lc.dat", dtype=str, delimiter=',')

instr_all = dat[:,7].astype(str)
choose = np.logical_or(instr_all=='"P48+ZTF"', instr_all=='"P60+SEDM"')

instr = instr_all[choose]
mjd = dat[:,1][choose]
filt = dat[:,2][choose]
m = dat[:,4][choose]
e_m = dat[:,5][choose]
limmag = dat[:,6][choose]

# get the LT photometry
filts = ['u', 'g', 'r', 'i', 'z']
for f in filts:
    print(f)
    dat = np.loadtxt(
            "ZTF18abukavn_%s_band_final.ascii" %f, 
            dtype=str, delimiter=',')
    nvals = len(dat[:,1])
    instr = np.append(instr, ['LT']*nvals)
    mjd = np.append(mjd, dat[:,0])
    filt = np.append(filt, ['%s' %f]*nvals)
    m = np.append(m, dat[:,1])
    e_m = np.append(e_m, dat[:,2])

# formatting
instr = np.char.strip(instr, '"')
mjd = mjd.astype(float)
filt = np.char.strip(filt, '"')
m = m.astype(float)
e_m = e_m.astype(float)

ascii.write(
        [instr, mjd, filt, m, e_m], 
        'ZTF18abukavn_opt_phot.dat',
        names=['Instrument', 'MJD', 'Filter', 'Mag', 'eMag'],
        overwrite=True)
