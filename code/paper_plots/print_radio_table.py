""" Print table of radio flux measurements """

import numpy as np
from astropy.time import Time
from astropy.cosmology import Planck15

d = Planck15.luminosity_distance(z=0.03154).cgs.value

headings = np.array(
        ['Date (UTC)', '$\Delta t$', 'Instrument', r'$\nu$ (GHz)', 
         'Flux $\mu$Jy', 'Luminosity ($10^27\,\erg\,\psec$)', 
         r'$\theta_\mathrm{FWHM}$',
         'Int. time (h)', 'Notes'])
label = "radio-flux"
caption = "Radio flux density measurements for SN2018gep. \
           An upper limit refers to the image RMS."

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

outputf = open("radio_table_%s.txt" %label, "w")
outputf.write("\\startlongtable \n")
outputf.write("\\begin{deluxetable}{%s} \n" %colstr)
outputf.write("\\tablecaption{%s\label{tab:%s}} \n" %(caption,label))
outputf.write("\\tablewidth{0pt} \n")
outputf.write("\\tablehead{ %s } \n" %colheadstr)
#outputf.write("\\rotate \n")
outputf.write("\\tabletypesize{\scriptsize} \n")
outputf.write("\startdata \n")

dat = np.loadtxt(
"/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/radio.txt",
delimiter=',', dtype=str)
date = Time(dat[:,0])
tel = dat[:,1]
nu = dat[:,2]
f = dat[:,3]
ef = dat[:,4]
nrows = dat.shape[0]

for ii in np.arange(nrows):
    # Print the date as is
    # Convert the date into a dt
    t0 = Time('2018-09-09')
    dt = int((date[ii]-t0).value)
    # Convert the flux into a fluxstr
    if '<' in f[ii]:
        # if upper limit, print as such
        fval = float(f[ii][1:])
        fstr = '$%s$' %f[ii]
    else:
        # if not an upper limit, include the uncertainty
        fstr = '$%s \pm %s$' %(f[ii], ef[ii])
        fval = float(f[ii])
    # Convert the flux into a luminosity
    lum = fval * 1E-6 * 1E-23 * 4 * np.pi * d**2 * float(nu[ii]) / 1E27
    lumstr = np.round(lum, 2)
    # Print row
    row = rowstr %(date[ii], dt, tel[ii], nu[ii], fstr, lumstr, "", "", "")
    outputf.write(row)

outputf.write("\enddata \n")
outputf.write("\end{deluxetable} \n")
