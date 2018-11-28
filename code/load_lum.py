import numpy as np

def load_lc():
    DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data"
    lsun = 3.839E33

    dat = np.loadtxt("%s/physevol.dat" %DATA_DIR, dtype=str)
    mjd = dat[:,0].astype(float)
    mjd0 = 58370.1473
    dt = mjd-mjd0
    lum = dat[:,8].astype(float) * lsun # original units L_sun
    llum = np.abs(dat[:,9].astype(float)*lsun)
    ulum = dat[:,10].astype(float)*lsun
    return dt, lum, llum, ulum 
