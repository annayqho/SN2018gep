""" Print table of ground-based opt measurements """

import numpy as np
from math import floor, log10
from astropy.time import Time
from astropy.cosmology import Planck15

def round_sig(x, sig=2):
    if x < 0:
        return -round(-x, sig-int(floor(log10(-x)))-1)
    return round(x, sig-int(floor(log10(x)))-1)

d = Planck15.luminosity_distance(z=0.03154).cgs.value

headings = np.array(
        ['Date (JD)', '$\Delta t$', 'Instrument', 'Filter', 
         'AB Mag', 'Error in AB Mag'])
label = "opt-phot"
caption = "Ground-based optical photometry for SN2018gep"

# Print the table headers
ncol = len(headings)
colstr = ""
colstr += 'l'
for col in np.arange(ncol-1): colstr+="r"
print(colstr)

colheadstr = ""
for col in np.arange(ncol-1):
    colheadstr += "\colhead{%s} & " %headings[col]
colheadstr += "\colhead{%s}" %headings[-1]

rowstr = ""
for col in np.arange(ncol-1):
    rowstr += "%s & "
rowstr += "%s \\\ \n"

outputf = open("table_%s.txt" %label, "w")
outputf.write("\\startlongtable \n")
outputf.write("\\begin{deluxetable}{%s} \n" %colstr)
outputf.write("\\tablecaption{%s\label{tab:%s}} \n" %(caption,label))
outputf.write("\\tablewidth{0pt} \n")
outputf.write("\\tablehead{ %s } \n" %colheadstr)
#outputf.write("\\rotate \n")
outputf.write("\\tabletypesize{\scriptsize} \n")
outputf.write("\startdata \n")

datadir = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/phot"
dat = np.loadtxt(datadir + "/ZTF18abukavn_opt_phot.dat",
delimiter=' ', dtype=str)
tel = dat[:,0]
mjd = Time(dat[:,1].astype(float), format='mjd')
filt = dat[:,2]
mag = dat[:,3].astype(float)
emag = dat[:,4].astype(float)
nrows = dat.shape[0]

for ii in np.arange(nrows):
    # Print the date as is
    # Convert the date into a dt
    t0 = Time(2458370.6473, format='mjd')
    dt = round_sig((mjd[ii]-t0).value, 1)
    # Convert the flux into a fluxstr
    if mag[ii] < 99.0:
        # If not an upper limit, print row
        row = rowstr %(mjd[ii], dt, tel[ii], filt[ii], mag[ii], emag[ii])
        outputf.write(row)

outputf.write("\enddata \n")
outputf.write("\end{deluxetable} \n")
