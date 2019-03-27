""" Load the individual light curve at a given band """

import numpy as np
from astropy.cosmology import Planck15
from astropy.io import ascii

zp = 2458370.6473
d = Planck15.luminosity_distance(z=0.03154).cgs.value


def get_uv_lc():
    DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/uv"

    # Dictionary mapping the filter to the wavelength
    bands = {}
    bands['V'] = 5468
    bands['B'] = 4392
    bands['U'] = 3465
    bands['UVW1'] = 2600
    bands['UVM2'] = 2246
    bands['UVW2'] = 1928

    f = DATA_DIR + "/UVOT_lightcurve_maghist_full_hostsub.ascii"
    dat = np.loadtxt(f, dtype=str, delimiter=' ')
    jd = dat[:,0].astype(float)
    filt = dat[:,1]
    wl = np.array([bands[val] for val in filt]) # ang
    nu = 3E18 / wl # Hz
    fnu_mjy = dat[:,2].astype(float)
    efnu_mjy = dat[:,3].astype(float)
    fact = 1E-3 * 1E-23 * 4 * np.pi * d**2 
    lnu = fnu_mjy * fact 
    elnu = efnu_mjy * fact 

    zp = 2458370.6473
    dt = jd-zp

    return dt, filt, fnu_mjy, efnu_mjy


def get_forced_phot():
    DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/phot"

    # Full light curve from Danny
    f = DATA_DIR + "/ZTF18abukavn.csv"
    dat = ascii.read(f)
    mjd = dat['mjd']
    jd = mjd + 2400000.5
    filt = dat['filter']
    mag = dat['mag']
    emag = dat['magerr']
    limmag = dat['lim_mag']
    code = dat['instrument']
    # sn_det = np.logical_and.reduce(
    #         (code=='ZTF Camera', dt > -1, ~np.isnan(mag)))
    # prog_det = np.logical_and.reduce(
    #         (code=='final photometry', dt < 0, ~np.isnan(mag)))
    # prog_nondet = np.logical_and(code=='final photometry', np.isnan(mag))
    return jd, filt, mag, emag, limmag, code


def get_lc():
    # get optical light curves 
    DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/phot"
    f = DATA_DIR + "/ZTF18abukavn_opt_phot.dat"
    dat = np.loadtxt(f, dtype=str, delimiter=' ')
    instr = dat[:,0]
    jd = dat[:,1].astype(float)
    filt = dat[:,2]
    mag = dat[:,3].astype(float)
    emag = dat[:,4].astype(float)
    dt = jd-zp
    
    # add the UV light curves
    add_dt, add_filt, fnu_mjy, efnu_mjy = get_uv_lc()
    # convert to AB mag
    add_mag = -2.5 * np.log10(fnu_mjy*1E-3) + 8.90
    add_emag = (efnu_mjy/fnu_mjy) # I think it's just the ratio
    choose = add_emag < 50
    dt = np.append(dt, add_dt[choose])
    filt = np.append(filt, add_filt[choose])
    mag = np.append(mag, add_mag[choose])
    emag = np.append(emag, add_emag[choose])

    return dt, filt, mag, emag
