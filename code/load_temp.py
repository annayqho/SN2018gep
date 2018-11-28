import numpy as np

def load_temp():
    DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data"
    dat = np.loadtxt("%s/physevol.dat" %DATA_DIR, dtype=str)
    dt = dat[:,1].astype(float)
    temp = dat[:,5].astype(float) * 1E3 # original units AU
    ltemp = np.abs(dat[:,6].astype(float)* 1E3)
    utemp = dat[:,7].astype(float)* 1E3
    return dt, temp, ltemp, utemp
