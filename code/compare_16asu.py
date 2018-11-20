""" compare light curves to PTF16asu """

import numpy as np
import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from astropy.table import Table
from astropy.cosmology import Planck15


DIST = Planck15.luminosity_distance(z=0.033).cgs.value


def optical():
    # PTF16asu
    dat = np.loadtxt("lc_16asu.txt", delimiter='&', dtype='str')
    t = dat[:,0].astype(float)
    dt = t-t[0]
    band = np.array([val.strip() for val in dat[:,2]])
    mag_raw = dat[:,3]
    mag = np.zeros(len(mag_raw))
    emag = np.zeros(len(mag))
    for ii,val in enumerate(mag_raw):
        if '>' not in val:
            mag[ii] = float(val.split('$pm$')[0])
            emag[ii] = float(val.split('$pm$')[1])

    choose = np.logical_and(mag > 0, band == 'g')
    Mag = mag[choose]-39.862 # distance modulus

    plt.errorbar(dt[choose], Mag, yerr=emag[choose], c='g', fmt='.')

    plt.gca().invert_yaxis()
    plt.tick_params(labelsize=14)
    plt.xlabel("dt [days]", fontsize=16)
    plt.ylabel("M", fontsize=16)

    plt.tight_layout()
    plt.show()



def radio():
    # radio
    dt = np.array([8, 8, 10, 10, 10, 10, 12])
    nu = np.array([4.8, 7.4, 6.1, 14, 20, 95, 6.1])
    f = np.array([79, 82, 33, 42, 90, 3600, 33]) # uJy
    ef = np.array([16, 15, -99, -99, -99, -99, -99])
    # -99 ef means upper limit

    f_cgs = f * 1E-6 * 1E-23 * 4 * np.pi * DIST**2

    # detections
    choose = ef > 0
    plt.errorbar(dt[choose], f_cgs[choose], yerr=ef[choose], fmt='.', c='k')

    # non-detections
    choose = ef < 0
    plt.scatter(dt[choose], f_cgs[choose], marker='v', c='k')

    # text saying the frequency
    for ii,val in enumerate(nu):
        if ii%2==0:
            plt.text(
                    dt[ii]*1.01, f_cgs[ii], "%s GHz" %val, fontsize=14,
                horizontalalignment='left', verticalalignment='top')
        else:
            plt.text(
                    dt[ii]*1.01, f_cgs[ii], "%s GHz" %val, fontsize=14,
                horizontalalignment='left', verticalalignment='bottom')

    
    # and now limits for this new source...
    # dist was Sept 9 I think
    # AMI: Sept 13
    d = Planck15.luminosity_distance(z=0.033).cgs.value
    plt.scatter(4, 35*1E-6*1E-23*4*np.pi*d**2, s=50, c='lightblue', marker='v')
    plt.text(4*1.01, 35*1E-6*1E-23*4*np.pi*d**2, "15.5 GHz", 
            horizontalalignment='left', fontsize=14)

    # SMA: 
    plt.scatter(4, 3.5*1E-3*1E-23*4*np.pi*d**2, s=50, c='lightblue', marker='v')
    plt.text(4*1.01, 3.5*1E-3*1E-23*4*np.pi*d**2, "231.5 GHz", 
            horizontalalignment='left', fontsize=14)

    dt = 6
    L = 0.59*1E-3*1E-23*4*np.pi*d**2
    plt.scatter(dt, L, s=50, c='lightblue', marker='v')
    plt.text(dt*1.01, L, "230 GHz", 
            horizontalalignment='left', fontsize=14)

    # VLA:
    plt.scatter(5, 91*1E-6*1E-23*4*np.pi*d**2, s=50, c='lightblue', marker='v')
    plt.text(5*1.01, 27*1E-6*1E-23*4*np.pi*d**2, "10 GHz", 
            horizontalalignment='left', fontsize=14)


    plt.yscale('log')
    plt.xlim(3,14)
    plt.tick_params(labelsize=14)
    plt.xlabel("dt [day]", fontsize=14)
    plt.ylabel(r"$L_\nu $[erg/s/Hz]", fontsize=14)
    plt.show()

    #plt.savefig("radio_comparison.png")


if __name__=="__main__":
    #optical()
    radio()
