"""
06-09-2013
@author: nmorenoc
         nicolas.morenochaparro@kasut.edu.sa

It reads a bunch of output files generated by lammps and extract the processors efficiency data to contruct efficiency plots and so. The calculations are based on a wall time of 24 hours. 

Input: initial anf final SIDs, fitting function: lineal, powerLaw or exp, file name of the plot (without extension)
It returns an array with processors density,and the corresponding number of steps reached in a 24h walltime. The las 4 returned values correspond to the fitting parameters, slope (B), intercept (A), correlation coefficient and standar error.

The regression can be based on a function lineal, powerlaw or logarithmic.
lineal:   NParticles/Nprocessors = A+TimeSteps*B
powerlaw: NParticles/Nprocessors = A*TimeSteps^B
exp:      NParticles/Nprocessors = A*exp(TimeSteps*B)
"""

import numpy as np
import itertools as it
import pandas as pd
from scipy import stats
import readThermo as rt
s

typeAnalysis = 'effic'

class analysis():

    def __init__(self, sid, fit='lineal'):
        self.sidi = sid[0]
        self.sidf = sid[1]
        self.fit - fit
    
    def run(init, fin, fit):

        EALL = np.array([], dtype=float)
        for j in range(init,fin):
            try:
                allThermo, name, e, cluster = rt.readThermo('%s/micelle.out' %j)  ##e contains: Time, processors, steps, atoms
            except: "no SID %s" %j
            else:
                if cluster == 'TACC':
                    EALL = np.append(EALL,e)


        sortedEff = EALL.reshape(-1,4).astype(float)
        secpday = 60*60*24

        ppp = sortedEff[1:,3]/sortedEff[1:,1]
        steps = sortedEff[1:,2]/sortedEff[1:,0]*secpday/.

        ppplog = np.log(ppp)
        stepslog = np.log(steps)

        if fit =='lineal':
            slope, intercept, r_value, p_value, std_err = stats.linregress(steps,ppp)
        elif fit == 'powerLaw':
            slope, intercept, r_value, p_value, std_err = stats.linregress(stepslog,ppplog)
            intercept = np.exp(intercept)
        elif fit == 'exp':
            slope, intercept, r_value, p_value, std_err = stats.linregress(steps,ppplog)
            intercept = np.exp(intercept)

        return ppp, steps, slope, intercept, r_value, std_err

    def test(init, fin, fit,fileName):
        ppp, steps, r_value, p_value, std_err= run(init, fin,fit)
        pa.plotTrajectories(ppp, steps, '${t_{steps}\\times 10^{-3}}/{day}$','$\\rho_{processor}$',fileName, [0.0,0], [0,0],True, '.',' ', None)
        #plt.plot(stepslog, interceptlog + slopelog*stepslog, '--')
        #plt.plot(stepslog, ppplog, '.')
        #plt.show()
