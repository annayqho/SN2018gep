import numpy as np

def load_lc():
    DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data"
    dat = np.loadtxt("%s/physevol.dat" %DATA_DIR, dtype=str)
    mjd = dat[:,0].astype(float)
    mjd0 = 58370.1473
    dt = mjd-mjd0
    lum = dat[:,8].astype(float) * 3.839E33 # original units L_sun
    llum = dat[:,9].astype(float) / dat[:,8].astype(float)
    ulum = dat[:,10].astype(float)  / dat[:,8].astype(float)
    elum = np.max((np.abs(llum), np.abs(ulum)), axis=0) # fractional unc
    return dt, lum, elum*lum # abs uncertainty
