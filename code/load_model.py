""" Load David Khatami's model """
import numpy as np

def load():
    modeldir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data"
    model = np.loadtxt(modeldir + "/2018gep_csm_model.dat")
    mod_dt = model[:,0]
    mod_lum = model[:,1]
    mod_rad = model[:,2]
    mod_temp = model[:,3]

    return mod_dt, mod_lum, mod_rad, mod_temp
