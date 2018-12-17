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
    # From Avishay
    dt = np.array([2.07, 4.25])
    vel = np.array([35000, 25000]) / 1E3
    plt.scatter(
            dt, vel, marker='s', facecolor='white', 
            edgecolor='k', s=100, lw=2)

    # From Ragnhild
    dt = np.array([4, 9, 11.5, 16.5, 22])
    vel = np.array([30000, 23000, 24000, 22000, 22000])/1E3
    evel = np.array([1000, 4000, 2000, 1000, 2000])/1E3
    el = np.array(['O', 'Fe', 'Fe', 'Fe', 'Fe'])

    plt.errorbar(
            dt[el=='Fe'], vel[el=='Fe'], yerr=evel[el=='Fe'], 
            fmt='s', c='k', label="AT2018gep", zorder=10, ms=10, lw=2)


def plot_16asu():
    """ values from Table 4 of Whitesides et al. 2017
    
    for the explosion date, they get 57518.53, 2016 May 10.53 (Section 3.1)
    let's use the Fe II velocity since that's the closest to what we have
    """
    # dt between May 10.53 and May 31.99 (inclusive) = 21.47
    dt = np.array(
            [24.97-10.53, 27.36-10.53, 21.47+4.39, 21.47+7.36, 21.47+10.42])
    vel = np.array(
            [28.3, 29.5, 25.7, 21.6, 22.0])*1000/1E3
    evel = np.array(
            [1.3, 1.4, 0.3, 0.4, 1.3])*1000/1E3
    plt.errorbar(dt, vel, yerr=evel, marker='.', c='purple', fmt='--')
    plt.text(dt[0]/1.03, vel[0], 'iPTF16asu',
            horizontalalignment='right', fontsize=12,
            verticalalignment='center')



def plot_grbsne():
    """ Modjaz et al. 2016 """
    dat = np.loadtxt(DATA_DIR + "/modjaz_vel.txt", dtype=str, delimiter=';')
    n = np.array(['sn1998bw', 'sn2006aj', 'sn2003lw', 'sn2010bh'])

    names = np.array([val.strip() for val in dat[:,0]])
    is_grbsn = np.array([val in n for val in names])

    name = names[is_grbsn]
    phase = dat[:,1][is_grbsn].astype(float)
    dt = np.zeros(len(phase))

    # correct phase to dt from GRB
    choose = name == 'sn1998bw'
    # explosion epoch: 2450929.41
    # epoch of max light: May 12.2 for V-band (Galama 1998)
    offset = 2450945.7-2450929.41
    dt[choose] = phase[choose] + offset

    choose = name == 'sn2006aj'
    # epoch of max: 2453794.7 (Modjaz 2014)
    # epoch of explosion: 2453784.649 (Campana 2006)
    offset = 2453794.7-2453784.649
    dt[choose] = phase[choose] + offset

    choose = name == 'sn2010bh' # 100316D
    # epoch of burst: March 16 2010
    # epoch of max: 8 days post-explosion (Bufano 2012)
    offset = 8
    dt[choose] = phase[choose] + offset

    choose = name == 'sn2003lw' # 030329
    # epoch of burst: 2452977.41769 (Malesani 2004)
    # epoch of V-band max:  16 days (not sure about this number)
    offset = 16
    dt[choose] = phase[choose] + offset

    vel = dat[:,2][is_grbsn].astype(float)*-1
    evel = dat[:,3][is_grbsn].astype(float)

    # To calculate the rolling mean, leave out 100316D
    tofit = name != 'sn2010bh'
    tfit = dt[tofit]
    vfit = vel[tofit]
    evfit = evel[tofit]

    dt_mean = np.arange(0,30,1)
    vel_mean = np.zeros(len(dt_mean))
    evel_mean = np.zeros(len(dt_mean))
    for ii,t in enumerate(dt_mean):
        choose = np.abs(tfit-t) <= 5
        if sum(choose) >= 3:
            w = 1/(evfit[choose])**2
            mean,wsum = np.average(
                    vfit[choose], weights=w, returned=True)
            emean = np.sqrt(1/np.sum(w))
            vel_mean[ii] = mean
            evel_mean[ii] = emean
        else:
            vel_mean[ii] = -100
            evel_mean[ii] = 1

    plt.fill_between(
            dt_mean, (vel_mean-evel_mean)/1E3, (vel_mean+evel_mean)/1E3,
            color='grey', alpha=0.5, label="LLGRB-SNe")


def plot_icbl():
    """ Modjaz et al. 2016 """
    dat = np.loadtxt(
            DATA_DIR + "/modjaz_vel.txt", dtype=str, delimiter=';')
    n = np.array(
            ['sn1997ef', 'sn2002ap', 'sn2003jd', 'sn2007D', 'sn2007bg', 
             'sn2007ru', 'sn2010ay', 'PTF10bzf', 'PTF10qts', 'PTF10vgv'])

    names = np.array([val.strip() for val in dat[:,0]])
    is_icbl = np.array([val in n for val in names])

    name = names[is_icbl]
    phase = dat[:,1][is_icbl].astype(float)
    vel = dat[:,2][is_icbl].astype(float)*-1
    evel = dat[:,3][is_icbl].astype(float)

    # for the dt, use the fact that Ic usually reach V-band max in 12-20 days
    # to add a sort of uncertainty to the dt, so offset would be 16 +/- 4 days
    dt = phase+16
    edt = 4

    dt_mean = np.arange(0,30,1)
    vel_mean = np.zeros(len(dt_mean))
    evel_mean = np.zeros(len(dt_mean))
    for ii,t in enumerate(dt_mean):
        choose = np.abs(dt-t) <= 5
        if sum(choose) >= 3:
            w = 1/(evel[choose])**2
            mean,wsum = np.average(
                    vel[choose], weights=w, returned=True)
            emean = np.sqrt(1/np.sum(w))
            vel_mean[ii] = mean
            evel_mean[ii] = emean
        else:
            vel_mean[ii] = -100
            evel_mean[ii] = 1

    plt.fill_between(
            dt_mean, (vel_mean-evel_mean)/1E3, (vel_mean+evel_mean)/1E3,
            color='purple', alpha=0.5, label="Ic-BL SNe")


def plot_ic():
    """ Modjaz et al. 2016 """
    dat = np.loadtxt(
            DATA_DIR + "/modjaz_vel.txt", dtype=str, delimiter=';')
    n = np.array(
            ['sn1983V', 'sn1990B', 'sn1992ar', 'sn1994I', 'sn2004dn',
             'sn2004fe', 'sn2004ge', 'sn2004gt', 'sn2005az', 'sn2004aw',
             'sn2005kl', 'sn2005mf', 'sn2007cl', 'sn2007gr', 'sn2011bm',
             'sn2013dk'])

    names = np.array([val.strip() for val in dat[:,0]])
    is_ic = np.array([val in n for val in names])

    name = names[is_ic]
    phase = dat[:,1][is_ic].astype(float)
    vel = dat[:,2][is_ic].astype(float)*-1
    evel = dat[:,4][is_ic].astype(float)

    # for the dt, use the fact that Ic usually reach V-band max in 12-20 days
    # to add a sort of uncertainty to the dt, so offset would be 16 +/- 4 days
    dt = phase+16
    edt = 4

    dt_mean = np.arange(0,30,1)
    vel_mean = np.zeros(len(dt_mean))
    evel_mean = np.zeros(len(dt_mean))
    for ii,t in enumerate(dt_mean):
        choose = np.abs(dt-t) <= 5
        if sum(choose) >= 3:
            w = 1/(evel[choose])**2
            mean,wsum = np.average(
                    vel[choose], weights=w, returned=True)
            emean = np.sqrt(1/np.sum(w))
            vel_mean[ii] = mean
            evel_mean[ii] = emean
        else:
            vel_mean[ii] = -100
            evel_mean[ii] = 1

    plt.fill_between(
            dt_mean[vel_mean>0], 
            (vel_mean[vel_mean>0]-evel_mean[vel_mean>0])/1E3, 
            (vel_mean[vel_mean>0]+evel_mean[vel_mean>0])/1E3,
            color='orange', alpha=0.5, label="Ic SNe")


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
            dt, vel/1E3, yerr=evel/1E3, fmt='.', ls='--', c='orange')
    plt.text(dt[0]/1.03, vel[0]/1E3, 'PTF12gzk',
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
    fig,ax = plt.subplots(1, 1, figsize=(8,5))

    plot_18gep()
    plot_grbsne()
    plot_16asu()
    plot_icbl()
    plot_ic()
    plot_12gzk()

    # Formatting
    plt.legend(fontsize=14, loc='lower left', ncol=2)
    plt.xlabel(r"$\Delta t$ [days]", fontsize=16)
    plt.ylabel(
            r"Fe II Velocity ($10^3$ km/s)", fontsize=16)
    #plt.yscale('log')
    #plt.xscale('log')
    plt.xlim(0, 30)
    plt.ylim(0, 40)
    plt.tick_params(axis='both', labelsize=16)
    plt.tight_layout()

    #plt.show()
    plt.savefig("vel.png")
