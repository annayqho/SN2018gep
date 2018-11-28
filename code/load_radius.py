import numpy as np

def load_radius():
    DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data"
    dat = np.loadtxt("%s/physevol.dat" %DATA_DIR, dtype=str)
    mjd = dat[:,0].astype(float)
    mjd0 = 58370.1473
    dt = mjd-mjd0
    Rsun = 6.955E10
    AU = 1.496e+13
    rad = dat[:,2].astype(float) * AU # original units AU
    lrad = dat[:,3].astype(float)*AU
    urad = dat[:,4].astype(float)*AU
    erad = np.max((np.abs(lrad), np.abs(urad)), axis=0) 
    return dt, rad, erad # abs uncertainty
