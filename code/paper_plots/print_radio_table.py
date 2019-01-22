""" Print table of radio flux measurements """

import numpy as np

headings = np.array(
        ['Date (UTC)', '$\Delta t$', 'Instrument', r'$\nu$ (GHz)', 
         'Flux $\mu$Jy', 'Luminosity (erg\,\psec)', r'$\theta_\mathrm{FWHM}$',
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

fstr = '$%s \pm %s$' %(f[ii], ef[ii])
row = rowstr %(np.round(t,2), 'ALMA', nu[ii], fstr)
outputf.write(row)

outputf.write("\enddata \n")
outputf.write("\end{deluxetable} \n")
