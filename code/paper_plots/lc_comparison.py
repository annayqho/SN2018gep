""" Do a light curve comparison """

from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
import numpy as np

DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/phot"

if __name__=="__main__":
    dat = np.loadtxt("../data/lc_16asu.txt", delimiter='&', dtype='str')
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

    plt.errorbar(dt[choose]-11, Mag, yerr=emag[choose], c='grey', fmt='o')
    plt.plot(dt[choose]-11, Mag, c='grey', label='iPTF16asu')

    f = DATA_DIR + "/ZTF18abukavn_g_band_final.ascii"
    band = f.split("_")[1]
    dat = np.loadtxt(f, dtype=str, delimiter=',')
    mjd = dat[:,0].astype(float)
    mag = dat[:,1].astype(float)
    emag = dat[:,2].astype(float)
    distmod = 35.87616973979363
    plt.errorbar(
            mjd-mjd[0], mag-distmod, emag, fmt='o', ms=5, c='k',
            mec='black', mfc='white', label='18gep, %s-band' %band)

    plt.gca().invert_yaxis()
    plt.tick_params(labelsize=14)
    plt.xlim(-5,20)
    plt.ylim(-17, -21)
    plt.xlabel("dt [days]", fontsize=16)
    plt.ylabel("Absolute Magnitude", fontsize=16)
    plt.legend(loc='upper right', fontsize=12)
    plt.tight_layout()

    #plt.show()
    plt.savefig("lc_comparison.png")
