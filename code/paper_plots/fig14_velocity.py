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
    """ V-band max is 3 days after t_0 """
    # From early spectra
    dt_early = np.array([1, 2, 4.2])
    vel_early = np.array([45000, 33000, 30000]) / 1E3
    plt.scatter(
            dt_early, vel_early, marker='s', facecolor='white', 
            edgecolor='k', s=100, lw=2)

    # From Ragnhild
    dt = np.array([4, 9, 11.5, 16.5, 22])
    vel = np.array([30000, 23000, 24000, 22000, 22000])/1E3
    evel = np.array([1000, 4000, 2000, 1000, 2000])/1E3
    el = np.array(['O', 'Fe', 'Fe', 'Fe', 'Fe'])

    plt.errorbar(
            dt[el=='Fe'], vel[el=='Fe'], yerr=evel[el=='Fe'], 
            fmt='s', c='k', label="SN2018gep", 
            zorder=10, ms=10, lw=2)

    plt.plot(
        np.hstack((dt_early, dt)), np.hstack((vel_early, vel)), 
        c='k', ls='-', lw=3, zorder=5)


def plot_16asu():
    """ values from Table 4 of Whitesides et al. 2017
    
    for the explosion date, they get 57518.53, 2016 May 10.53 (Section 3.1)
    let's use the Fe II velocity since that's the closest to what we have
    the uncertainty on their explosion date is 0.17d
    """
    # dt between May 10.53 and May 31.99 (inclusive) = 21.47
    dt = np.array(
            [24.97-10.53, 27.36-10.53, 21.47+4.39, 21.47+7.36, 21.47+10.42])
    vel = np.array(
            [28.3, 29.5, 25.7, 21.6, 22.0])*1000/1E3
    evel = np.array(
            [1.3, 1.4, 0.3, 0.4, 1.3])*1000/1E3
    plt.errorbar(
            dt, vel, xerr=0.17, yerr=evel, ms=10,
            marker='o', c='#84206b', fmt='-', lw=3, zorder=10)
    plt.text(dt[1]*1.1, vel[1], 'iPTF16asu',
            horizontalalignment='left', fontsize=12,
            verticalalignment='bottom')


def plot_1998bw():
    """ Modjaz et al. 2016
    offset from Galama 1998 """
    dat = np.loadtxt(DATA_DIR + "/modjaz_vel.txt", dtype=str, delimiter=';')
    names = np.array([val.strip() for val in dat[:,0]])

    name = 'sn1998bw'
    choose = names==name
    phase = dat[:,1][choose].astype(float)
    vel = dat[:,2][choose].astype(float)*-1
    evel = dat[:,3][choose].astype(float)
    offset = 2450945.7-2450929.41
    dt = phase+offset
    #plt.errorbar(dt, vel/1E3, yerr=evel/1E3, marker='v', c='#f6d746')
    plt.plot(dt, vel/1E3, c='#f6d746', lw=3, alpha=0.5)
    plt.text(
            dt[8], vel[8]/1E3, 'SN1998bw', fontsize=12,
            verticalalignment='bottom', horizontalalignment='left')


def plot_2006aj():
    """ Modjaz et al. 2016
    offset from Campana 2006 """
    dat = np.loadtxt(DATA_DIR + "/modjaz_vel.txt", dtype=str, delimiter=';')
    names = np.array([val.strip() for val in dat[:,0]])

    name = 'sn2006aj'
    choose = names==name
    phase = dat[:,1][choose].astype(float)
    vel = dat[:,2][choose].astype(float)*-1
    evel = dat[:,3][choose].astype(float)
    offset = 2453794.7-2453784.649
    dt = phase+offset
    plt.plot(
            dt, vel/1E3, c='#f6d746', lw=3, alpha=0.5)
    plt.text(
            dt[0], vel[0]/1E3/1.05, 'SN2006aj', fontsize=12,
            verticalalignment='top', horizontalalignment='center')


def plot_2010bh():
    """ Modjaz et al. 2016
    offset from Bufano 2012 """
    dat = np.loadtxt(DATA_DIR + "/modjaz_vel.txt", dtype=str, delimiter=';')
    names = np.array([val.strip() for val in dat[:,0]])

    name = 'sn2010bh'
    choose = names==name
    phase = dat[:,1][choose].astype(float)
    vel = dat[:,2][choose].astype(float)*-1
    evel = dat[:,3][choose].astype(float)
    offset = 8 
    dt = phase+offset
    #plt.errorbar(
    #        dt, vel/1E3, yerr=evel/1E3, marker='D', c='#f6d746', alpha=0.3)
    plt.errorbar(
            dt, vel/1E3, c='#f6d746', alpha=0.5, lw=3)
    plt.text(
            dt[3], vel[3]/1E3, 'SN2010bh', fontsize=12,
            verticalalignment='bottom', horizontalalignment='left')


def plot_2003lw():
    """ Modjaz et al. 2016
    offset from Malesani 2014"""
    dat = np.loadtxt(DATA_DIR + "/modjaz_vel.txt", dtype=str, delimiter=';')
    names = np.array([val.strip() for val in dat[:,0]])

    name = 'sn2003lw'
    choose = names==name
    phase = dat[:,1][choose].astype(float)
    vel = dat[:,2][choose].astype(float)*-1
    evel = dat[:,3][choose].astype(float)
    offset = 8 
    dt = phase+offset
    plt.plot(
            dt, vel/1E3, c='#f6d746', lw=3, alpha=0.5)
    #plt.errorbar(
    #        dt, vel/1E3, yerr=evel/1E3, marker='D', c='#f6d746', alpha=0.3)
    plt.text(
            dt[0], vel[0]/1E3, 'SN2003lw', fontsize=12,
            verticalalignment='bottom', horizontalalignment='left',
            label="LLGRB-SNe")


def grb171205a():
    # For GRB171205A, the g-band max is about 12 days after the GRB
    dat = np.loadtxt(
            DATA_DIR + "/grb171205a_vel.txt", dtype=float, delimiter=',')
    dt = dat[:,0]
    vel = dat[:,1]/1E3
    #plt.scatter(dt, vel, marker='o', c='grey', alpha=0.5)
    plt.plot(dt, vel, ls='-', c='#f6d746', lw=3, alpha=0.5, label="LLGRB-SN")
    plt.text(dt[0], vel[0], 'SN2017iuk', fontsize=12)


def plot_12gzk():
    """ Modjaz et al. 2016 """
    dat = np.loadtxt(
            DATA_DIR + "/modjaz_vel.txt", dtype=str, delimiter=';')
    n = 'PTF12gzk'

    names = np.array([val.strip() for val in dat[:,0]])
    choose = np.where(names==n)[0]

    name = names[choose]
    phase = dat[:,1][choose].astype(float)
    vel = dat[:,2][choose].astype(float)*-1
    evel = dat[:,4][choose].astype(float)

    # for the dt, use the fact that Ic usually reach V-band max in 12-20 days
    # to add a sort of uncertainty to the dt, so offset would be 16 +/- 4 days
    dt = phase+16
    plt.errorbar(
            dt, vel/1E3, yerr=evel/1E3, fmt='.', ls='--', c='#f6d746')
    plt.text(dt[1], vel[1]/1E3, 'PTF12gzk',
            horizontalalignment='right', fontsize=12,
            verticalalignment='bottom')



def plot_population():
    """ Modjaz et al. 2016 """
    inputf = pyfits.open(DATA_DIR + "/asu.fit")
    dat = inputf[1].data
    inputf.close()
    names = dat['SN']
    phase = dat['Phase']
    vel = -1*dat['vabs']/1E3

    # Rel SN, no GRB
    name = 'sn2009bb'
    # explosion epoch: March 19
    # maximum light: Apr 1 (Soderberg 2009)
    offset = 13
    choose = names == name
    dt = phase[choose]+offset
    v = vel[choose]
    plt.plot(dt, v, c='grey')
    plt.scatter(dt, v, c='grey')
    plt.text(dt[2], v[2]/1.02, 'SN2009bb',
            horizontalalignment='right', fontsize=12,
            verticalalignment='top', zorder=5)

    name = 'sn2012ap'
    # explosion date: 2012 Feb 5
    # maximum light: let's use the B-band max from Milisavljevic 2014,
    # which is on Feb 18.2, so that's 13 days
    offset = 13
    choose = names == name
    dt = phase[choose]+offset
    v = vel[choose]
    plt.plot(dt, v, c='grey')
    plt.scatter(dt, v, c='grey')
    plt.text(dt[0], v[0], 'SN2012ap',
            horizontalalignment='right', fontsize=12,
            verticalalignment='center', zorder=5)

    # iPTF17cw
    # max light: 57760.666 is 8 days before maximum light
    # explosion time set to be 57750.552
    # so it's 18 days to maximum
    dt = np.array([17, 43]) + 18 # since maximum light
    v = np.array([17300, 17500]) / 1E3
    plt.plot(dt, v, c='grey')
    plt.scatter(dt, v, c='grey')
    plt.text(dt[0], v[0]/1.03, 'iPTF17cw',
            horizontalalignment='left', fontsize=12,
            verticalalignment='top', zorder=5)

    # 2008D: measured using He I 5876 lines
    dt = np.array(
            [5.75, 6.77, 9.00, 11.66,  
             22.78, 23.88, 30.56, 32.61, 
             33.87, 37.79, 49.81, 60.67])
    v = np.array(
            [14200, 13300, 12300, 11500, 
             10600, 10600, 10500, 10400, 
             10300, 10000, 9300, 8600]) /1E3
    plt.scatter(dt, v, c='grey')
    plt.plot(dt, v, c='grey')
    plt.text(dt[0], v[0], '2008D(HeI)',
            horizontalalignment='right', fontsize=12,
            verticalalignment='top', zorder=5)



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
    plt.plot(dt, v, c='#57106e', label='SN Ic')
    plt.scatter(dt, v, c='#57106e')
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
    plt.plot(dt, v, c='#57106e')
    plt.scatter(dt, v, c='#57106e')
    plt.text(dt[0], v[0], 'SN2005az',
            horizontalalignment='right', fontsize=12,
            verticalalignment='center')



if __name__=="__main__":
    fig,ax = plt.subplots(1, 1, figsize=(6,5))

    plot_18gep()
    plot_16asu()
    grb171205a()
    plot_1998bw()
    plot_2006aj()
    plot_2003lw()
    plot_2010bh()
    #plot_12gzk()

    # Formatting
    plt.legend(fontsize=14, loc='upper right', ncol=1)
    plt.xlabel(r"Time since explosion (days)", fontsize=16)
    plt.ylabel(
            r"Fe II Velocity ($10^3$ km/s)", fontsize=16)
    plt.yscale('log')
    #plt.xscale('log')
    plt.xlim(0, 25)
    plt.ylim(10, 100)
    plt.tick_params(axis='both', labelsize=16)
    plt.tight_layout()

    #plt.show()
    plt.savefig("vel.eps", format='eps', dpi=1000)
