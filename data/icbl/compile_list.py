""" Compile a list of Ic-BL SNe with z < 0.1 """

import numpy as np
from astropy.table import Table


def opensn():
    """ The link to the catalog is here: https://sne.space/
    I searched for "BL" and got rid of anything that was not Ic-BL
    (e.g. "blue" or "variable"). I also deleted anything that did not
    have a known position. The problem is that there is no reported redshift
    for these sources. """
    dat = Table.read("open_sn_catalog.csv", format="csv", delimiter=',')
    nsn = len(dat)
    ra = []
    dec = []

    for ii,val in enumerate(dat['R.A.']):
        ra.append(str(val).split(',')[0])
        dec.append(str(dat['Dec.'][ii]).split(',')[0])

def ptf():
    """ the PTF/iPTF sample of 34 Ic-BL SNe
    I copied the table directly from the .tex file downloaded from the arXiv,
    then ran the following two commands
    %s/\\//g
    %s/ //g
    %s/\*//g
    %s/xx//g
    I also removed the commented-out lines
    """
    dat = Table.read(
            "taddia2018.dat", delimiter='&', format='ascii.fast_no_header')
    name = dat['col1']
    ra = dat['col2']
    dec = dat['col3']
    z = dat['col5']
    choose = z <= 0.1
    return name[choose], ra[choose], dec[choose], z[choose]

def sdssII():
    """ the sample from SDSS-II (Taddia et al. 2015)
    I copied the table from the .tex file, and only kept the Ic-BL ones
    also got rid of the IIb that was commented out
    """
    dat = Table.read(
            "taddia2015.dat", delimiter='&', format='ascii.fast_no_header')
    name = dat['col1']
    ra = dat['col2']
    dec = dat['col3']
    z = dat['col5']
    choose = z <= 0.1 # only one of them!
    return name[choose], ra[choose], dec[choose], z[choose]

def cano2013():
    """ The table of GRB/XRF-less Ic-BL SNe,
    from Cano et al. 2013 (since the ones with GRBs """
    dat = Table.read(
            "cano2013.dat", delimiter='&', format='ascii.fast_no_header')
    name = dat['col1']
    z = dat['col3']
    choose = z <= 0.1 
    return name[choose], z[choose]

def cano2016():
    """ The table of GRB-SNe from Cano et al. 2016 """
    dat = Table.read(
            "cano2016.dat", delimiter='&', format='ascii.fast_no_header')
    name = dat['col1'].filled()
    z = dat['col4'].filled()
    choose = z <= 0.1 
    return name[choose], z[choose]

def lyman2016():
    """ The list of Ic-BL SNe from Lyman et al. 2016 """
    dat = Table.read(
            "lyman2016.dat", delimiter='&', format='ascii.fast_no_header')
    name = dat['col1']
    z = dat['col3']
    choose = z <= 0.1 
    return name[choose], z[choose]

def prentice2016():
    """ The list of Ic-BL SNe from Prentice et al. 2016 """
    dat = Table.read(
            "prentice2016.dat", delimiter='&', format='ascii.fast_no_header')
    name = dat['col1']
    cl = dat['col2']
    z = dat['col4']
    is_icbl = np.logical_or(cl=='Ic-BL', cl=='GRB-SN')
    choose = np.logical_and(is_icbl, z <= 0.1)
    return name[choose], z[choose]

def modjaz2016():
    """ The list of Ic-BL SNe from Modjaz et al. 2016 
    Removed the one with a lower limit on redshift """
    dat = Table.read(
            "modjaz.dat", delimiter='&', format='ascii.fast_no_header')
    name = dat['col1']
    z = dat['col2']
    choose = z <= 0.1
    return name[choose], z[choose]

def ztf():
    """ The list of Ic-BL discovered in ZTF """
    dat = Table.read(
            "ztf.dat", delimiter='&', format='ascii.fast_no_header')
    name = dat['col1']
    ra = dat['col2']
    dec = dat['col3']
    z = dat['col4']
    choose = z <= 0.1
    return name[choose], ra[choose], dec[choose], z[choose]

if __name__=="__main__":
    name, z = modjaz2016()
