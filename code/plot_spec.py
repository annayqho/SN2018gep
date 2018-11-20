import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)



def early():
    dat = np.loadtxt("spec_sep11.txt")
    z = 0.033
    wl = dat[:,0]/(1+z)
    flux = dat[:,1]*10**15
    plt.plot(wl, flux, c='blue', label="-1 d")

    dat = np.loadtxt("spec_sep13.txt")
    z = 0.033
    wl = dat[:,0]/(1+z)
    flux = dat[:,1]*10**15
    plt.plot(wl, flux, c='red', label="+1 d")

    plt.xlabel("Rest Wavelength (Ang)", fontsize=14)
    plt.ylabel(r"Flux [$\times10^{-15}$] erg\,s$^{-1}$\,cm$^{-2}$", fontsize=14)
    plt.tick_params(axis='both', labelsize=12)
    plt.xlim(1000, 10000)

    plt.legend(loc='upper right', fontsize=14)

    plt.show()
    #plt.savefig("early_spec.png")


def late():
    dat = np.loadtxt("spec_sep17.txt")
    z = 0.033
    wl = dat[:,0]/(1+z)
    flux = dat[:,1]*10**15
    plt.plot(wl, flux, c='blue', label="-1 d")
    plt.ylim(0, 3)

    plt.xlabel("Rest Wavelength (Ang)", fontsize=14)
    plt.ylabel(r"Flux [$\times10^{-15}$] erg\,s$^{-1}$\,cm$^{-2}$", fontsize=14)
    plt.tick_params(axis='both', labelsize=12)

    plt.show()
    #plt.savefig("late_spec.png")


if __name__=="__main__":
    early()
