""" Print table of radio flux measurements """

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
        ['Start Time', '$\Delta t$', 'Instrument', r'$\nu$', 
         r'$f_\nu$', r'$L_\nu$', r'$\theta_\mathrm{FWHM}$',
         'Int. time (h)'])
subheadings = np.array(
        ['(UTC)', '(days)', '', '(GHz)', '($\mu$Jy)', '(erg\,\psec\,\phz)',
         '$^{\prime\prime}$', '(hr)']
label = "radio-flux"
caption = "Radio flux density measurements for SN2018gep."

# Print the table headers
ncol = len(headings)
colstr = ""
colstr += 'l'
for col in np.arange(ncol-1): colstr+="r"

colheadstr = ""
for col in np.arange(ncol-1):
    colheadstr += "\colhead{%s} & " %headings[col]
colheadstr += "\colhead{%s}" %headings[-1]

colsubheadstr = ""
for col in np.arange(ncol-1):
    colsubheadstr += "\colhead{%s} & " %subheadings[col]
colsubheadstr += "\colhead{%s}" %subheadings[-1]

rowstr = ""
for col in np.arange(ncol-1):
    rowstr += "%s & "
rowstr += "%s \\\ \n"

outputf = open("%s.txt" %label, "w")
outputf.write("\\startlongtable \n")
outputf.write("\\begin{deluxetable}{%s} \n" %colstr)
outputf.write("\\tablecaption{%s\label{tab:%s}} \n" %(caption,label))
outputf.write("\\tablewidth{0pt} \n")
outputf.write("\\tablehead{ %s \\ %s} \n" %(colheadstr,colsubheadstr))
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
majbeam = dat[:,5]
minbeam = dat[:,6]
inttime = dat[:,7]
nrows = dat.shape[0]

for ii in np.arange(nrows):
    # Print the date as is
    # Convert the date into a dt
    t0 = Time('2018-09-09T03:32')
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
    # Convert the flux into a luminosity density
    lum = fval * 1E-6 * 1E-23 * 4 * np.pi * d**2 
    if '<' in f[ii]:
        lumstr = "$<%s$" %round_sig(lum,2)
    else:
        lumstr = round_sig(lum, 2)
    # Print row
    row = rowstr %(
            date[ii], dt, tel[ii], nu[ii], fstr, lumstr, 
            '%s x %s' %(majbeam[ii],minbeam[ii]), inttime[ii])
    outputf.write(row)

outputf.write("\enddata \n")
outputf.write("\tablecomments{For VLA measurements: The quoted errors are calculated as the quadrature sums of the image rms, plus a 5\% nominal absolute flux calibration uncertainty. When the peak flux density within the circular region is less than three times the RMS, we report an upper limit equal to three times the RMS of the image. For AMI measurements: non-detections are reported as 3-$\sigma$ upper limits. For SMA measurements: non-detections are reported as a 1-$\sigma$ upper limit.")
outputf.write("\end{deluxetable} \n")
