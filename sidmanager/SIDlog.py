"""
03-29-2013
@author: nmorenoc
         nicolas.morenochaparro@kasut.edu.sa

It createsa latex file with the simulations details I've run so far.
This script is based on: http://users.obs.carnegiescience.edu/~cburns/site/?p=22, and uses the Table script from this source.

"""

#import Table
import sys, os
import numpy as np


def header(fout, flag):
    if flag =='1':   ### if flag is zero, then no header I think is not being used, because the header is not too big and the table is always generated full from SID.sid
        o = open(fout, 'w')
        print >>o, "\documentclass[8pt]{article}"
        print >>o,"\\usepackage[english]{babel}"
        print >>o,"\\usepackage{graphicx}"
        print >>o,"\\usepackage{url}"
        print >>o,"\\usepackage[landscape, left= 0.5in, legalpaper]{geometry}"
        print >>o,"\\usepackage{lastpage}"
        print >>o,"\\usepackage{deluxetable}"
        print >>o,"\\usepackage[section]{placeins}  %% This is a package I put to avoid error due to the use of too many figures"
        print >>o,"\\usepackage{caption}"
        print >>o,"\\usepackage{subcaption} %this is to use the subfigure options"
        print >>o,"\\usepackage{amsmath}"
        print >>o,"\\usepackage{amssymb}"

        print >>o,"%%\\usepackage{anysize}"
        print >>o,"%%\marginsize{0.5cm}{0.5cm}{1cm}{1cm}  %%%\marginsize{left}{right}{top}{bottom}"

        print >>o,"\n \n \n"
        print >>o,"\\begin{document}"
        print >>o,"\\title{Summary of SID carried out for DBC project}"
        print >>o,"\\author{Moreno-Chaparro Nicolas, Nunes P. Suzana, Calo Victor M.}"
        print >>o,"\date{April 26th 2013}"
        print >>o,"\maketitle"
        print >>o,"\n \n \n"        


#def createTable(fout, c, justs, header, variables):
#
#    ###http://en.wikibooks.org/wiki/LaTeX/Tables  Info about tables in latex
#    t = Table.Table(variables, justs=justs, caption='Simulations ID details', label="tab:SID")
#    t.add_header_row(header)
#    t.add_data(c, sigfigs=0)
#    t.print_table(fout)
#    print >>fout, "\end{document}"
#    fout.close()

def read(file, spacer):
    data = []
    f = open(file, 'r')
    for i in f:
        data.append(i.split(spacer))
    f.close()
    return np.array(data)


def doLog(sidFile, Texfile, flag, spacer='\t'):
    all = read(sidFile, spacer)
    header(Texfile, flag)
    headerRow = all[0]
    variables = len(headerRow)
    

    justs = np.chararray(variables)
    justs[0] = 'c'
    justs[1:] = 'r'
    justs = justs.tostring()

    c = all[1:].T
    c = c.tolist()
    #fout = open(Texfile,'a')
    #createTable(fout,c, justs, headerRow, variables)



if sys.argv[0] == "SIDlog.py":
    file = sys.argv[2]
    sidFile = sys.argv[1]
    flag = sys.argv[3]
    try:
        spacer = sys.argv[4]
    except: pass
    else: spacer = sys.argv[4]

    doLog(sidFile, file, flag, spacer)
else: print( "sidLog: It suppose to be running from another script?" )

