""" See whether the LC is consistent with being powered by Ni decay at
late times """

import numpy as np
from scipy.integrate import cumtrapz
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/code")
import matplotlib.pyplot as plt
from load_lum import load_lc


def qgamma(t, mni):
    """
    t: time in days
    mni: nickel mass in solar masses
    """
    return (mni)*(6.45*np.exp(-t/8.76) + 1.38*np.exp(-t/111.4))*1E43


def qpos(t, mni):
    """
    t: time in days
    mni: nickel mass in solar masses
    """
    return 4.64 * (mni) * (-np.exp(-t/8.76) + np.exp(-t/111.4)) * 1E41


def get_qdep(t, t0, mni):
    """ The total energy deposition rate
    """
    return qgamma(t, mni) * (1-np.exp(-(t0/t)**2)) + qpos(t, mni)


if __name__=="__main__":
    # load the bolometric light curve
    dt, lum, llum, ulum = load_lc()

    # integral of bol LC over time
    IntLtdt=cumtrapz((lum*86400)*dt,dt,initial=0)

    # lum ratio over time
    lum_ratio = (lum*86400)/IntLtdt
    plt.scatter(dt, dt**(2.5) * lum_ratio, c='k')
    
    # now, calculate the energy deposition integral
    mni = 1 # doesn't matter
    for t0 in np.arange(10,20,2):
        dt_qdep = np.linspace(0,100,10000)
        qdep = get_qdep(dt_qdep, t0, mni)
        qdep_ratio = qdep / cumtrapz(qdep*dt_qdep, dt_qdep, initial=0)
        plt.plot(dt_qdep, dt_qdep**(2.5) * qdep_ratio, label="t0=%s" %t0) 

    plt.ylim(0,10)
    plt.xlim(0,50)
    plt.legend()
    plt.show()
