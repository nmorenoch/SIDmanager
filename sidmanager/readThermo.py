"""
05-24-2013
@author: Nicolas Moreno Ch.
e-mail: nicolas.morenochaparro@kaust.edu.sa

It reads log or out files with the thermo information saved in lammps simulations. It mainly reads the thermo information and the data related with the efficiency of that run: Processor, ttime, steps and atoms. 
"""

import numpy as np
import sys, re


def readThermo(thermoFile):
    all=[]   ##all the data
    names = []
    effic = []  ##Data to be used in efficiency analisys Total time, in N processor, for M steps with Np particles

    flagData = False
    flagmin  = False
    f = open(thermoFile, 'r')
    cluster = 'WS'
    for line in f:
        l = line
        if re.match('Setting up minimization',l) != None:
            flagmin = True
        if re.match('TACC',l) != None:
            cluster = 'TACC'
        if re.match('WARNING: ',l) != None:
            flagData = False 
        if re.match('Loop ',l )  !=None:
            flagData = False
            y = l.split(' ')
            if flagmin ==False:
                effic.append([y[3], y[5], y[8], y[11]]) ##Time, processors, steps, atoms
            flagmin = False
        if (flagData == True):    
            y = map(float, l.split())
            if len(y)>6:
                all.append(y) #If is not the minimization info
        if re.match('Step ', l) != None:   
            flagData = True
            names = l.split(' ')
        
    f.close()

    e     = np.array(effic)
    allAr = np.array(all)
    
    return allAr, names, e, cluster

def readVoronoiVol(voronoiVolFile):
    flagData = False
    f = open(voronoiVolFile, 'r')
    all = []
    for line in f:
        l = line
        if (flagData == True):    
            y = map(float, l.split())
            all.append(y)
        if re.match('# TimeStep ', l) != None:   
            flagData = True

    f.close()

    allAr = np.array(all)
    
    return allAr

