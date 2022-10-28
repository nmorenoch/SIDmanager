"""
@author: Nicolas Moreno Chaparro

This script generates system set up for rigid core brush structure. To model Si_2-PMMA membranes

r_c =rwerew
"""

import sys, os, shutil
import numpy as np
from random import sample
#import pylab as plt 

##Constants
K_Boltz = 1.3806503e-23 ## Boltzmman constant m2 kg s-2 K-1
NAvog = 6.0221415e23  ## Avogadro Number


###Default values for a micelle simulation - DPD
rc = 1.0
rc_map = 1.0
m = 1.0
T = 1.0
rho = 3.0
pFrac = 0.15        ## polymer fraction
blockFracc = 0.25   ## PVP fraction
beadChain = 16  ## Number of beads
lx,ly,lz = 75,75,75   ## box size
delta_t = 0.04

#####Spring and bending harmonics parameters
ks = 50.0 ## spring constant
kb = 0.0  ## bending constant (angle restriction)
ro = 0.75
angle = 180

#####Physical data fof the particular system

mwBlockA = 65000.0 ##molecular weight block A g/mol
mwBlockB = 175000.0 ##molecular weight block B g/mol

mwA = 105.0   #4VP
mwB = 104.0   #Styrene

mVolA = 92.44  #molar volume cm^3/mol polymer  From Paper of P4VP-PS modeling
mVolB = 96.98  #molar volume cm^3/mol polymer
densPoly = np.array([1.15,1.06])  #density p4vp and PS g/cm^3


mul = 1.   ###Defining a multiplicative to the dispersion contribution, when the plot of Rs is donce
R = 8.6  ##Solubiloty radius for PMMA

deltaA = np.array([mul*18.6, 10.5, 7.5, 18.58])   ###Solubility parameter for the component PMMA: van der Walls, polar, hydrogen bonds, kind of radius Average. \delta= 22.64 calculated by Hansen contribution
deltaB = np.array([mul*18.6, 1.0, 4.1, 19.1])
#deltaB = np.array([mul*21.3, 5.8, 4.3, 19.1])

            #Solvents
mwTHF = 72.11   
mwDMF = 73.09
mwDOX = 88.11
mwBut = 86.10
mwMpy = 99.13
mwAce = 58.08


mwSolv   = np.array([72.11, 73.09, 88.11]) ###molecular weight of the 3 solvents g/mol
densSolv = np.array([0.8892,0.948,1.033, 1.13, 1.03, 0.79]) ###densities of the 3 solvents g/ml
                                          

deltaTHF = np.array([mul*16.8, 5.7, 8.0, 19.5])
deltaDMF = np.array([mul*17.4, 13.7, 11.3, 24.9])
deltaDOX = np.array([mul*19.0, 1.8, 7.4, 20.5])
deltaBut = np.array([mul*19.0, 16.6, 7.4, 26.2])
deltaMPy = np.array([mul*18.0, 12.3, 7.2, 22.9])

deltaAce = np.array([mul*15.5, 10.4, 7.0, 20.1])
deltaWater = np.array([mul*15.5, 16.0, 42.3, 47.8])



wFracTHF = 0.333
wFracDMF = 0.333
wFracDOX = 1.0 - wFracTHF - wFracDMF


class dataFile():

    def __init__(self, file, nPart, nBonds, nAngles, aType, bType, anType, xlo, xhi, ylo, yhi, zlo, zhi, m, rho, nSpheres, fgraf, polyL, ro, bondScal, sType, sBeads):

        self.f = oFile
      
        self.nBonds = nBonds
        self.nAngles = nAngles
        self.aType = aType
        self.bType = bType
        self.anType = anType
        self.xlo = xlo
        self.ylo = ylo
        self.zlo = zlo
        self.xhi = xhi
        self.yhi = yhi
        self.zhi = zhi
        self.m = m
        self.rho = rho
        self.sphereBeads = 0
        self.fgraf = fgraf
        self.polyL = polyL
        self.ro = ro
        self.bondScal = bondScal
        #self.totalBeads = int((self.xhi-self.xlo)*(self.xhi-self.xlo)*(self.xhi-self.xlo)*self.rho)
        self.nPart = int((self.xhi-self.xlo)*(self.xhi-self.xlo)*(self.xhi-self.xlo)*self.rho)

        self.nSpheres = nSpheres
        self.fracNPbeads = nSpheres*sBeads/float(self.nPart)
        self.polyFrac    = fgraf*self.nSpheres*sBeads*polyL/float(self.nPart)
        
        self.setup = "# Number of spheres: %d  -  Grafting Fraction: %f -  Polymer length: %d  -   Bond distance: %f  -  Bond scaling: %f  -  Particle N. Density: %f  -  Type of sphere %d  -  Fract. NP beads: %f  -  Fract. polymer: %f" % (nSpheres, fgraf, polyL, ro, bondScal, rho, sType, self.fracNPbeads, self.polyFrac)

    def header(self, inFile, inEdge):

          f = open(self.f, 'w')
          print >>f,"# Spherical shape generator"
          print >>f,"# @Author: Nicolas Moreno"
          print >>f,self.setup

          print >>f, " "
          print >>f, self.nPart.__str__()+" atoms"
          print >>f, self.nBonds.__str__()+" bonds"
          print >>f, self.nAngles.__str__()+" angles"
          print >>f, "0 dihedrals"
          print >>f, "0 impropers"

          print >>f, " "
          print >>f, self.aType.__str__()+" atom types"
          print >>f, self.bType.__str__()+" bond types"
          print >>f, self.anType.__str__()+" angle types"
          print >>f, "0 dihedral types"
          print >>f, "0 improper types"

          #Box limits
          print >>f, " "
          print >>f, self.xlo.__str__()+" "+self.xhi.__str__()+" xlo xhi"
          print >>f, self.ylo.__str__()+" "+self.yhi.__str__()+" ylo yhi"
          print >>f, self.zlo.__str__()+" "+self.zhi.__str__()+" zlo zhi"

          #Masses
          print >>f, " "
          print >>f,"Masses \n"
          for i in range(len(self.m)):
              print >>f, (i+1).__str__()+" "+self.m[i].__str__()

          print >>f,"\n Atoms"
          print >>f, " "

          self.readPointDataFromBlender(inFile, inEdge, f)
          f.close()


    def readPointDataFromBlender(self, vertFile, edgeFile, f):
        iVer = open(vertFile, 'r')
        iEd = open(edgeFile, 'r')
        radius = 3
        a=[]
        array = np.array
        arrayBond = np.array

        for line in iVer:
            l = line
            y = map(float, l.split())
            array = np.append(array,y)

        for line in iEd:
            l = line
            y = map(float, l.split())
            arrayBond = np.append(arrayBond,y)


        array = np.delete(array,0)
        arrayBond = np.delete(arrayBond, 0)
        array = array.reshape(-1,3)
        arrayBond = arrayBond.reshape(-1,2)

        arrayBondChains = []
        arrayGraftPoint = []

        polyChains = len(array)*self.fgraf   ##number of chains per nano particle
        polyBeads = polyChains*self.polyL

        polyFract = self.nSpheres*float(polyBeads)/float(self.nPart)
        print 'polymer fraction', polyFract
        
        k=0
            #print array
        x, y, z = self.randomXYZarray(self.nSpheres)
        for ns in range(self.nSpheres):
            grp = self.fgraf*len(array)
            rp = np.zeros(len(array))
            index = self.randomID(len(array), int(grp))
            rp[index] = 1
            
            for j in range(len(array)):
                k+=1
                px = array[j][0]*radius+x[ns]
                py = array[j][1]*radius+y[ns]
                pz = array[j][2]*radius+z[ns]
                print >>f, "%d %d 1 %f %f %f" % ((k),ns+1,px,py,pz)
                #ax,ay,az = self.randomXYZarray(1)
                #polyx = ax #px
                #polyy = ay #py
                #polyz = az #pz
                grp-=1
                flag = np.random.randint(2)
                #flag  = np.mod(j,1./self.fgraf)  ##flag is 0 when the current j is a multiple of the current grafting fraction

                if rp[j]==1:   ##If grafting polymer chain in this particle vertice
                    arrayGraftPoint.append([k, px, py, pz, ns+1, x[ns], y[ns], z[ns]]) ##x,y,z are the center coordinates for the current sphere
                    
        for gp in arrayGraftPoint:
            #Printing the polymer chains coordinates
             polyx = gp[1]
             polyy = gp[2]
             polyz = gp[3]
             sign = [1,1,1]

             
             if polyx<gp[5]:
                 sign[0] = -1
             if polyy<gp[6]:
                 sign[1] = -1
             if polyz<gp[7]:
                 sign[2] = -1
             
             molec = gp[4]
             prev  = gp[0]  ##The ID of the point where the polymer chain will be grafted
             type  = 1  #Bond type equals to the beads interparticle.
             for p in range(self.polyL):
                 temp = self.randomPositionAtRo(self.ro/self.bondScal) ##getting coord for a random sphere of ro over bondScaling value in order to avoid overlaping and box crossing of the polyme chains
                 polyx = polyx+sign[0]*temp[0]
                 polyy = polyy+sign[1]*temp[1]
                 polyz = polyz+sign[2]*temp[2]
                 k+=1
                 print >>f, "%d %d 2 %f %f %f" % ((k),molec,polyx, polyy, polyz)

                 type = 2   ## polymer - np beads use a particluar spring type 2 

                 if p>0:
                     prev = k-1
                     type = 3  ## The polymer spring type will be 3 for polymer-polymer beads, and 2 for polymer-particle beads.
                 arrayBondChains.append([prev,k, type])  ##loading an array with the conectivity info for the polymer chains


        arrayBondChains = np.array(arrayBondChains)
        
        self.sphereBeads = k
        
        solventBead = int((self.xhi-self.xlo)*(self.xhi-self.xlo)*(self.xhi-self.xlo)*self.rho - self.sphereBeads)
        print solventBead, k, self.nPart, len(arrayBond), 'solventbeads  spherebeads  totalBeads   spherebonds'
        x, y, z = self.randomXYZarray(solventBead)

        for p in range(0,solventBead):
            print >>f, "%d %d 3 %f %f %f" % (self.sphereBeads+p+1, self.nSpheres+p+1,x[p],y[p],z[p])

        k=0
        print >>f,"\n Bonds \n"
        for ns in range(self.nSpheres):
            part = len(array)*ns    ##This part value takes into account that there are more than one sphere, and shift by this value the global index of the particle in the bond topology construction
            for j in range(len(arrayBond)):
                k+=1
                print >>f, "%d 1 %d %d" %((k), (arrayBond[j][0]+1)+part, (arrayBond[j][1]+1)+part)
        # Adding the information for bonds of polymer chains                
        for bond in arrayBondChains:  
            k+=1
            print >>f, "%d %d %d %d" %((k), bond[2],bond[0], bond[1])
        
    def randomID(self, N, k):
        r = sample(xrange(N), k)  ##N is the total number of active points, and k is the number of grafted points
        return r

    def randomPositionAtRo(self, ro):
        x = np.random.rand(1)*(ro-ro/2)
        y = np.sqrt(np.random.rand(1)*(ro**2-x**2))
        z = np.sqrt(ro**2-x**2-y**2)
        return x,y,z

    def randomXYZarray(self, randomPoints):
        print randomPoints
        x = np.random.rand(randomPoints)*(self.xhi-self.xlo)+self.xlo
        y = np.random.rand(randomPoints)*(self.yhi-self.ylo)+self.ylo
        z = np.random.rand(randomPoints)*(self.zhi-self.zlo)+self.zlo

        return x, y, z
        
    #oFile = sys.argv[1]  #outputfile
    #sType = (sys.argv[2])  # Sphere type
    #particles = int(sys.argv[3]) # number of nanoparticles
    #grafFrac = float(sys.argv[4]) #grafting fraction
    #plen = int(sys.argv[5])  #polymer length

lx = 25.
ly = 25.
lz = 25.

#if sType=='1':
#    nbonds = 30
#    sBeads = 12
#if sType=='2':
#    nbonds = 120
#    sBeads = 42
#if sType=='3':
#    nbonds = 480
#    sBeads = 162
#if sType=='4': 
#    nbonds = 1920
#    sBeads = 642


#polymerBonds = plen*grafFrac*sBeads   #plymerLenght*Grafting FRaction*sphere beads

#totalB = int((polymerBonds+nbonds)*particles)

SID = np.array([31,41], dtype=int)

nofs = SID[1]-SID[0]

sType = np.ones(nofs)*2    # Sphere type
#particles = np.array([20, 40, 60, 80, 100, 100, 100, 100, 100, 100]) # number of nanoparticles
particles = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]) # number of nanoparticles
grafFrac = np.ones(nofs) #grafting fraction
plen = np.ones(nofs)*45  #polymer length

for f in range(nofs):
    os.system("mkdir %s"% (SID[0]+f))
    oFile = "%s/data_%s.gn" %(SID[0]+f, SID[0]+f)   #outputfile
    lx = 20.
    ly = 20.
    lz = 20.

    if sType[f]==1:
        nbonds = 30
        sBeads = 12
    if sType[f]==2:
        nbonds = 120
        sBeads = 42
    if sType[f]==3:
        nbonds = 480
        sBeads = 162
    if sType[f]==4: 
        nbonds = 1920
        sBeads = 642
    
    polymerBonds = plen[f]*grafFrac[f]*sBeads   #plymerLenght*Grafting FRaction*sphere beads
    totalB = int((polymerBonds+nbonds)*particles[f])

    data = dataFile(oFile, 197, totalB, 0, 3, 3, 1, -lx, lx, -ly, ly, -lz, lz, [1,1,1], 3, particles[f], grafFrac[f], int(plen[f]), 0.75, 3., int(sType[f]), sBeads)
    data.header('nanoP'+int(sType[f]).__str__(), 'nanoB'+int(sType[f]).__str__())




