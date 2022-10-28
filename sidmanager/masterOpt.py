"""
08-01-2013 (Optimized)
@author: Nicolas Moreno Ch.
e-mail: nicolas.morenochaparro@kaust.edu.sa

It executes all the postprocessing implemented methods for a set of simulations ID.
Input: Initial simulation ID, Final simulation ID, type of run, and folder where the simulations results are located. Where type of run can be a particular analysis or can be set as 'all' in order to run all the scripts
Output: Latex and dvi report file. Denpending on the type of run this script generates diferent outputs as histograms, trajectories and plain text. All of this information is saved in the related simulation ID folder.


"""

import os, sys, glob
import numpy as np
import pandas as pd


import readDump as rd
import gofr as gofr
import msd as msd
import vacf as vacf
import angacf as angacf
import rOfG as rog
import thermoVars as tv

#import structureFactor as sf
import dbScanRead as dbsR


import micelleScaling as mSca
import plotAll as pa


class master():

    def __init__(self, sids, types, analysis, step = 0):
        self.lista = {'angacf':angacf, 'gofr':gofr, 'msd':msd, 'vacf':vacf, 'rog':rog, 'dbscan':dbsR, 'thermo':tv}  ##Any additional analysis test should be updated in this list.
        self.sids = sids
        self.type = types   
        self.pwd  = os.getcwd()
        self.analysis = analysis ##list of values with the setup for each analysis. First entry is the analysis name, second is a shift value for the time steps, and the third entry is a list with the labels for the plots
        ##Creating list to store as many snapshot as needed.
        self.a = []
        self.m = []
        self.t = []
        self.x, self.y, self.z = [], [], []
        self.lbox = []
        self.snapshots = 0  ##Counter for the number of snapshots or frames dumped
        self.step = step
        #self.info = pd.read_table('SIDlog/SID.sid', sep='\t', header=0, index_col=0)  ##Loading data from SID log file
        self.stepsCounter = []
    
    """
    This function load the data of each snapshot in global variables, sorted by particle ID
    """
    def readData(self, dumpFile):
        step = dumpFile.split('.')[1]  ##Assuming the current standard naming of the files
        #print step, 'step in readData'
#    for dumpFile in sorted(glob.glob(os.path.join(self.folderID+'/dumps/', '*.dump')), key=os.path.getmtime):
        #   i+=1.
        #    r = i/bin
        #    if (int(r)-r) == 0 and points<max:
        #        points+=1
        #        print dumpFile
        #        a,m, t, x, y,z, lbox = rd.read(dumpFile)
        #        index = np.argsort(a)   ##sorted index according to the ID name fro each particle
        #        x = x[index]
        #        y = y[index]
        #        z = z[index]
        #        xcom, ycom, zcom = self.comi(x), self.comi(y), self.comi(z)
        #        X.append([x,y,z,xcom, ycom, zcom])
        #M=points
        try:
            data = rd.read(dumpFile)
        except:
            print ("Dump file %s corrupted or empty" % step)
        else:
            index = np.argsort(data[0])
            self.a.append(data[0][index])
            self.m.append(data[1][index])
            self.t.append(data[2][index])
            self.x.append(data[3][index])
            self.y.append(data[4][index])
            self.z.append(data[5][index])
            self.lbox.append(data[6])
            self.snapshots+=1
            self.stepsCounter.append(step)
            #print self.snapshots, 'snpashot in readData'

    def iterateFolder(self, functions, dumps):
        for dumpFile in dumps:          ### For loop over the dump files in a folder 
            functions(dumpFile)

    def restartValues(self):
        self.a = []
        self.m = []
        self.t = []
        self.x, self.y, self.z = [], [], []
        self.lbox = []
        self.snapshots = 0  ##Counter for the number of snapshots or frames dumped
        self.step = 0
        self.stepsCounter = []

    listDumps = lambda self, SID, folder: sorted(glob.glob(os.path.join("%s/%s/" % (SID, folder), '*')), key=os.path.getmtime)
    #listDumps = lambda self, SID, folder: sorted(glob.glob(os.path.join("%s/%s/" % (SID, "allDumps"), '*')), key=os.path.getmtime)   ##Modified on 09/08/2014 comment this line in order to make it worl properly

    ##This function verify if the dump file information needs to be loaded.
    def loadDumps(self, sid):
        flagLoad = False
        for i in self.analysis:
            if self.lista[i[0]].typeAnalysis == 'traj' or self.lista[i[0]].typeAnalysis == 'single': 
                flagLoad = True
        if flagLoad:
            try:
                dumps = self.listDumps(sid,'dumps')
            except: 
                print ('SID %s or dump does not exist' % sid)
            else: 
                print ("Processing SID: %s" % sid)
                self.loaded = True
                #print dumps, 'here dumps in loadDumps'
                
            if self.step == 0:  #If there is not an specific time step to process
                self.iterateFolder(self.readData, dumps)  ##This line load all the snapshots of sid
            else:
                dumpFile = ("%s" % (dumps[self.step]))
                self.readData(dumpFile) #only one snapshot will be analysed
                
    
    def run(self):
        lista = self.lista
        for sid in self.sids:

            self.loadDumps(sid)

            for t in self.analysis:
                aType = lista[t[0]].typeAnalysis
            
               
                print ("Doing %s" % (t[0]))
                
                if lista[t[0]].typeAnalysis == 'traj':
                    test = lista[t[0]].analysis([self.a, self.m, self.t, self.x, self.y, self.z, self.lbox])
                    out = test.run()
                    self.testLabels(t,test)

                    self.saveAnalysis(out, sid, t[0], test.label, test.name, test.indepVar) # This function plots and saves as plain text data the result of each analysis

                elif lista[t[0]].typeAnalysis == 'single':
                    out = []
                    for step in range(self.snapshots):
                        # print step, 'step'
                        test = lista[t[0]].analysis([self.a[step], self.m[step], self.t[step], self.x[step], self.y[step], self.z[step], self.lbox[step]], t[1])
                        
                        out.append(test.run())
                        #print test, 'this is test'
                        self.testLabels(t,test)
                    self.saveAnalysis(np.array(out, dtype=object), sid, t[0], test.label, test.name, test.indepVar, test.limsup) # This function plots and saves as plain text data the result of each analysis

                else: # lista[t[0]].typeAnalysis == 'indep' or lista[t[0]].typeAnalysis == 'thermo':
                    try:
                        test = lista[t[0]].analysis(sid)
                        out = test.run()
                    except Exception, e:
                        print str(e)
                        print ('There is a problem running analysis %s for sid %s, maybe the required input file does not exist.' % (sid, t[0]))
                    else:
                        self.testLabels(t,test)
                        self.saveAnalysis(out, sid, t[0], test.label, test.name, test.indepVar) # This function plots and saves as plain text data the result of each analysis

            self.restartValues() #It is necesary to restart the particles arrays and global counters

    def testLabels(self, t, test):
        if len(t)>2:
            test.label, test.name =  t[1], t[2] 
            #return self.lista[t[0]].label, self.lista[t[0]].name

    
    def runComparison(self):
        for t in self.analysis:
            readAnalysis(self, t, self.SID)
            #plot()
            #save()

        v=[]
        x=[]
        label =[]
        vmean=[]
        stds=[]

        for j in compSID:
            allThermo, name, e = rt.readThermo('%s/micelle.out' %j)   
            compVar = 'Temp'
            compVar = name.index(compVar) ##Index of this varsr
            v.append(allThermo.T[compVar])#/(nch[j]*pl[j]))
            x.append(allThermo.T[0])
            vmean.append(allThermo.T[compVar].mean())#/(nch[j]*pl[j]))
            stds.append(allThermo.T[compVar].std())#/(nch[j]*pl[j]))
            label.append(j)
            v=np.array(v)
            x=np.array(x)
            stds=np.array(stds)
            pa.plotTrajectories(v.T, x.T, 'Time', name[compVar], '%s%s_%s'% (name[compVar], compSID[0],compSID[-1]), 0,True, 'o', label)

            er = vmean+stds
            er2 = vmean-stds
            vmean = np.array(vmean).T
            label = np.array(label).T

        pa.plotTrajectories(vmean, label, 'SID', name[compVar], '%s%s_%s_mean'% (name[compVar], compSID[0],compSID[-1]), 0,True, 'o', ['mean'], stds)
    
        print 'Finish comparison'
    
    def saveAnalysis(self, data, sid, name, labels, captions, indepVar, max=0):
        
        self.saveBinary(sid,name,[data, labels, captions, indepVar, max])

        types = self.lista[name].typeAnalysis
        print types
        if types =='traj':
            pa.plotTrajectories(data[1], data[0], labels[0], labels[1], ("%s/data/%s"%(sid,captions[1])), [0,max], [0,0], True,'o','-','')
        elif types == 'single':
            temp = data.T
            
            for i in range(len(temp)):
                if temp[i][0] != None: #So I define this None flag when combining arrays and single values.
                    try:
                        pa.plotTrajectories(temp[i], self.stepsCounter, 't', labels[i+1], ("%s/data/%s"%(sid,captions[i+1])), [0,0], [0,0], True, '.', ' ', '') ##In this case I dont define the limit in the source scirpt.
                    except:
                        d = len(temp[i][0])
                        spaceTraj = np.zeros(d)
                        totalSteps = len(self.stepsCounter)
                        for j in range(totalSteps):
                            spaceTraj=temp[i][j]
                            #spaceTraj /= (1.*totalSteps)
                        pa.plotTrajectories(spaceTraj, indepVar, labels[0], labels[i+1], ("%s/data/%s"%(sid,captions[i+1])), [0,max], [0,0], True, '.', ' ', '')

        elif types == 'indep':
            #print len(data), data[0], len(data[0])
            for j in data:
                    for i in range(len(j)):
                        try:
                            pa.plotTrajectories(j[i], indepVar[i], labels[0], labels[i+1], ("%s/data/%s"%(sid,captions[i+1])), [0,0], [0,0], True, '.', ' ', '') ##In this case I dont define the limit in the source scirpt.
                        except:
                            print 'error'
        elif types =='thermo':
            #print len(data)
            for i in range(len(data)):
                pa.plotTrajectories(data[i], indepVar, labels[0], labels[i+1], ("%s/data/%s"%(sid,captions[i+1])), [0,0], [0,0], True, '.', ' ', '') ##In this case I dont define the limit in the source scirpt
                
        elif types =='clust':
            shift = 2  # shift value to start plotting. I should put this value as an input from the constructor.
                       #print len(data)
            for i in range(len(data)):
                try:
                    pa.plotTrajectories(data[i][shift:], indepVar[i][shift:], labels[0], labels[i+1], ("%s/data/%s"%(sid,captions[i+1])), [0,0], [0,0], True, '.', ' ', '') ##In this case I dont define the limit in the source scirpt
                except:
                    pa.plotHistogram(data[i], 50,  labels[i+1], "Number of micelles", ("%s/data/%s"%(sid,captions[i+1])), False, indepVar[i])
             
        print ("Saving results of %s for SID: %s" % (name, sid))

    
    def plot(self, vaf, t, fft, freqs, max):
        self.id = "msd"+self.id
        fftID = "FFTmsd"+self.id
        #                  array, xdata, xlabel, ylabel, name, maxX, maxY, grid, marker, linestyle,legend,stds='a',exten='.eps'
        pa.plotTrajectories(vaf, t, 't',"msd(t)", self.id, [0,max], [0,0], True, 'o', '-', '')
        pa.plotTrajectories(fft, freqs, 'freq',"fft(C)", fftID, [0, 0], [0, 0], True, 'o', '-', '')

        pa.plotTrajectories(allThermo.T[i,:], allThermo.T[0,:], 'Time', name[i], '%s/%s'% (j,name[i]), [0,0], [0,0],True, ' ','-', [j])


    def efficiency(self):
            EALL = np.append(EALL,e)
            
            sortedEff = EALL.reshape(-1,4).astype(float)
            secpday = 60*60*24

            rho = sortedEff[1:,3]/sortedEff[1:,1]
            ti = sortedEff[1:,2]/sortedEff[1:,0]*secpday/1000.

        
    def saveBinary(self, sid, tag,data):
        datadir = ("%s/data"%sid)
        if not os.path.exists(datadir):
            os.makedirs(datadir)
        
        oFile = ("%s/%s.npy"% (datadir,tag))
        o = open(oFile, 'w')
        np.save(o,data)
        #times.tofile(o, sep=" ")
        #o.write('\n')
        #nofparticles = len(data[0])   ####This should be the total number of particlese
        #for i in range(0,nofparticles):
        #    data[:,i].tofile(o, sep=" ", format="%f")
        #    o.write('\n')
        o.close()
        
    def test(self):
        return 0
