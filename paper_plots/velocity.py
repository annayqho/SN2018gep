""" Plot of velocity vs. time """

import glob
import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import numpy as np
from astropy.io import fits as pyfits
from astropy.table import Table
from astropy.cosmology import Planck15

DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data"

def plot_18gep():
    """ Temporary values from Ragnhild """
    dt = np.array([4, 9, 11.5, 16.5, 22])
    vel = np.array([30000, 23000, 24000, 22000, 22000])/1E3
    evel = np.array([1000, 4000, 2000, 1000, 2000])/1E3
    el = np.array(['O', 'Fe', 'Fe', 'Fe', 'Fe'])

    #plt.errorbar(
    #        dt[el=='O'], vel[el=='O'], yerr=evel[el=='O'], fmt='o', c='k')
    plt.errorbar(
            dt[el=='Fe'], vel[el=='Fe'], yerr=evel[el=='Fe'], 
            fmt='s', c='k', label="AT2018gep", zorder=10)


def plot_16asu():
    """ values from Table 4 of Whitesides et al. 2017
    
    for the explosion date, they get 57518.53, 2016 May 10.53 (Section 3.1)
    let's use the Fe II velocity since that's the closest to what we have
    """
    # dt between May 10.53 and May 31.99 (inclusive) = 21.47
    dt = np.array(
            [24.97-10.53, 27.36-10.53, 21.47+4.39, 21.47+7.36, 21.47+10.42])
    print(dt)
    vel = np.array(
            [28.3, 29.5, 25.7, 21.6, 22.0])*1000/1E3
    print(vel)
    evel = np.array(
            [1.3, 1.4, 0.3, 0.4, 1.3])*1000/1E3
    #plt.errorbar(dt, vel, yerr=evel, fmt='s', c='grey', label="iPTF16asu")
    plt.plot(dt, vel, c='grey', label='Rel SN, no GRB')
    plt.text(dt[0], vel[0], 'iPTF16asu',
            horizontalalignment='right', fontsize=12,
            verticalalignment='center')


def plot_population():
    """ Modjaz et al. 2016 """
    inputf = pyfits.open(DATA_DIR + "/asu.fit")
    dat = inputf[1].data
    inputf.close()
    names = dat['SN']
    phase = dat['Phase']
    vel = -1*dat['vabs']/1E3

    # Rel SN, LLGRB
    name = 'sn2006aj'
    # epoch of max: 2453794.7 (Modjaz 2014)
    # epoch of explosion: 2453784.649 (Campana 2006)
    offset = 2453794.7-2453784.649
    choose = names == name
    dt = phase[choose]+offset
    v = vel[choose]
    plt.plot(dt, v, c='magenta', label="LLGRB-SN")
    plt.text(dt[-1], v[-1], 'SN2006aj',
            horizontalalignment='left', fontsize=12,
            verticalalignment='center')

    name = 'sn1998bw'
    # explosion epoch: 2450929.41
    # epoch of max light: May 12.2 for V-band (Galama 1998)
    offset = 2450945.7-2450929.41
    choose = names == name
    dt = phase[choose] + offset
    v = vel[choose]
    plt.plot(dt, v, c='magenta')
    plt.text(dt[0], v[0], 'SN1998bw',
            horizontalalignment='right', fontsize=12,
            verticalalignment='center')


    # Rel SN, no GRB
    name = 'sn2009bb'
    # explosion epoch: March 19
    # maximum light: Apr 1 (Soderberg 2009)
    offset = 13
    choose = names == name
    dt = phase[choose]+offset
    v = vel[choose]
    plt.plot(dt, v, c='grey')
    plt.text(dt[0], v[0], 'SN2009bb',
            horizontalalignment='right', fontsize=12,
            verticalalignment='center')

    # regular Ic
    name = 'sn2007gr'
    # time of max: 2454338.5
    # time of discovery: 2007 Aug 15.51, expl must be < 5 days before disc
    # arbitrary rise time of 13.5 is used in Mazzali 2010,
    # do we use that here
    offset = 13.5
    choose = names == name
    dt = phase[choose]+offset
    v = vel[choose]
    plt.plot(dt, v, c='lightblue', label='SN Ic')
    plt.text(dt[0], v[0], 'SN2007gr',
            horizontalalignment='right', fontsize=12,
            verticalalignment='center')


    # regular Ic
    name = 'sn2005az'
    # time of max: 2453473
    # arbitrary rise time of 13.5
    offset = 13.5
    choose = names == name
    dt = phase[choose]+offset
    v = vel[choose]
    plt.plot(dt, v, c='lightblue')
    plt.text(dt[0], v[0], 'SN2005az',
            horizontalalignment='right', fontsize=12,
            verticalalignment='center')



if __name__=="__main__":
    plot_18gep()
    plot_16asu()
    plot_population()

    # Formatting
    plt.legend(fontsize=14, loc='lower left', ncol=2)
    plt.xlabel(r"$\Delta t$ [days]", fontsize=14)
    plt.ylabel(
            r"Fe II Velocity ($10^3$ km/s)", fontsize=14)
    #plt.yscale('log')
    #plt.xscale('log')
    plt.xlim(0, 30)
    plt.ylim(0, 35)
    plt.tick_params(axis='both', labelsize=14)
    plt.tight_layout()

    #plt.show()
    plt.savefig("vel.png")
