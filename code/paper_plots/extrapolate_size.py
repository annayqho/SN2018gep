""" Use the radius evolution to infer the size at the time of explosion """

import numpy as np
import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from scipy.integrate import quad
from scipy.optimize import curve_fit
from astropy.table import Table


def load_radius():
    DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data"
    dat = np.loadtxt("%s/physevol.dat" %DATA_DIR, dtype=str)
    dt = dat[:,1].astype(float)
    Rsun = 6.955E10
    AU = 1.496e+13
    rad = dat[:,2].astype(float) * AU # original units AU
    lrad = dat[:,3].astype(float)*AU
    urad = dat[:,4].astype(float)*AU
    erad = np.max((np.abs(lrad), np.abs(urad)), axis=0) 
    return dt, rad, erad # abs uncertainty


if __name__=="__main__":
    dt, rad, erad = load_radius()

    fig,ax = plt.subplots(1,1,figsize=(7,5))
    ax.errorbar(
            dt, rad, yerr=erad, fmt='o', 
            mec='k', mfc='k', ms=5, c='k')

    # fit a quadratic to the first few points
    choose = dt <2
    out = np.polyfit(dt[choose], rad[choose], deg=2)
    xfit = np.linspace(0, 1.3, 1000)
    yfit = out[0]*xfit**2  +out[1]*xfit + out[2]
    ax.plot(xfit, yfit, c='k', ls='--')

    # Functional form:
    # -2.63E14 * x**2 + 6.12E14 * x + 3.21E14
    ax.axhline(y=3.21E14, ls='--', c='k')
    ax.text(6, 3.5E14, r"$R_0 = 3.2 \times 10^{14}\,\mathrm{cm}$", fontsize=12)

    ax.set_ylabel("Photospheric Radius (cm)", fontsize=16)
    ax.set_xlabel(
        r"Days since $t_0=$JD 2458370.6473", fontsize=16)
    ax.yaxis.set_tick_params(labelsize=14)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.set_yscale('log')
    ax.set_xscale('linear')
    ax.set_ylim(1E14, 1E16)
    ax.set_xlim(-1, 10)
    #plt.show()
    plt.savefig("radfit.png")

