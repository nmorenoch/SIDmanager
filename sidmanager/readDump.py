"""
Generic script to read the data at the dump file. This data can be from a position, velocity or force dump file.

parameters: dump file, shifted. Shifted 0 or 1, depending of the module to used readDump, the particles position falling out of the box will be shifted in order to fit during binning process. By default not shifted.

output: returns the dumped variables stored as a numpy arrays.

"""
def readDumpFile(fileName,skiprows=9,cols=12):
    import pandas as pd
    filed = fileName
    d = pd.read_csv(filed, skiprows=skiprows,sep=" ",header=None,usecols=range(cols))
    d =d.sort_values(0)
    return d.values


def readOLD(f, shifted=0):

    import numpy as np
    import sys, os, re

    #o   = sys.argv[2]


    ##Just the input dump file, and eventually an output file can be created
    dump   = open(f, 'r')
    #output = open(o, 'w')

    lbox = []


    ##Flags for the dump file to read the data
    flagAtoms = False
    flagBoxB  = False

    all=[]   ##all the data

    for line in dump:
        l = line
        if (flagAtoms == True):    
            y = map(float, l.split())
            all.append(y)
            # output.write('%g %g %g %g\n' % (y[3],y[4],y[5], y[2]))
        if re.match('ITEM: ATOMS', l) != None:   
            flagAtoms = True
            flagBoxB  = False   #stop reading box bounds 

        if (flagBoxB == True):    
            y = map(float, l.split())
            len = np.abs(y[0]) + y[1] 
            lbox.append(len)
        
        if re.match('ITEM: BOX BOUNDS', l) != None:   
            flagBoxB = True

    dump.close()
    #print lbox, 'lbox'
    L = np.array(lbox)

    allAr=np.array(all)


    pIndex = allAr[:,0]  ##particle index
    mIndex = allAr[:,1]  ##molecule index
    type   = allAr[:,2]    ## type of particle

    #Particles coordinates
    x = allAr[:,3]
    y = allAr[:,4]
    z = allAr[:,5]

    #shifted = 1  ## This variable is used when boundary isses arise in old lammps saved files
    if shifted ==1:
        l=L[0]/2
        lz=L[2]/2
        x[x<=-l] = -l+0.01
        y[y<=-l] = -l+0.01
        z[z<=-lz] =  -lz+0.01
    
        x[x>=l] =  l-0.01
        y[y>=l] =  l-0.01
        z[z>=lz] =  lz-0.01

    return pIndex, mIndex, type, x, y, z, L

## get methods. In case I made this a class, I wont need to keep this methods, since I can call the variables directly

def getpIndex():
    return pIndex
def getmIndex():
    return mIndex
def getType():
    return type
def getX():
    return X
def getY():
    return Y
def getZ():
    return Z
    
    



