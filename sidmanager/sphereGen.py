"""
@author: Nicolas Moreno Chaparro

This script generates a lammps data file with the coordinates of a sphere, based on a vtk file with the strucutre and the conectivity.


"""



import sys, os, re
import numpy as np



#f = sys.open(inFile, 'r')
#o = sys.open(oFile, 'w')


class dataFile():

    def __init__(self, file, nPart, nBonds, nAngles, aType, bType, anType, xlo, xhi, ylo, yhi, zlo, zhi, m, rho, nSpheres):

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
        #self.totalBeads = int((self.xhi-self.xlo)*(self.xhi-self.xlo)*(self.xhi-self.xlo)*self.rho)
        self.nPart = int((self.xhi-self.xlo)*(self.xhi-self.xlo)*(self.xhi-self.xlo)*self.rho)

        self.nSpheres = nSpheres

    def header(self, inFile, inEdge):

          f = open(self.f, 'w')
          print >>f,"# Spherical shape generator"
          print >>f,"# @Author: Nicolas Moreno"

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

          #self.readPointDataFromVTK(inFile,f)
          self.readPointDataFromBlender(inFile, inEdge, f)
          f.close()


    def readPointDataFromVTK(self, vtkFile,f):
        i = open(vtkFile, 'r')
        flag = False
        flagBond = False
        a=[]
        array = np.array
        arrayBond = np.array
        
        for line in i:
            l = line
           
            if re.match('POINT_DATA', l) != None:
                flagBond = False
                
            if flagBond == True:
                y = map(float, l.split())
                arrayBond = np.append(arrayBond,y)

            if re.match('POLYGONS', l) != None:
                flag = False
                flagBond = True
              
            if flag == True:
                y = map(float, l.split())
                array = np.append(array,y)

            if re.match('POINTS', l) != None:
                flag = True
                self.sphereBeads = int(l.split()[1])
            
                #array = np.array(a)
        
        array = np.delete(array,0)
        arrayBond = np.delete(arrayBond, 0)
        array = array.reshape(-1,3)
        arrayBond = arrayBond.reshape(-1,4)

        
       #print array
        for j in range(len(array)):
            print >>f, "%d 1 1 %f %f %f" % ((j+1),array[j][0],array[j][1],array[j][2])

        
        solventBead = int((self.xhi-self.xlo)*(self.xhi-self.xlo)*(self.xhi-self.xlo)*self.rho - self.sphereBeads)
        x = np.random.rand(solventBead)*(self.xhi-self.xlo)+self.xlo
        y = np.random.rand(solventBead)*(self.yhi-self.ylo)+self.ylo
        z = np.random.rand(solventBead)*(self.zhi-self.zlo)+self.zlo
          
        for p in range(0,solventBead):
            print >>f, "%d 2 2 %f %f %f" % (self.sphereBeads+p+1,x[p],y[p],z[p])
        
        print >>f,"\n Bonds \n"
        for j in range(len(arrayBond)):
            if (arrayBond[j][3] == 0):
                print >>f, "%d 1 %d %d" %((2*j+1), (arrayBond[j][1]+1), (arrayBond[j][2]+1))
                print >>f, "%d 1 %d %d" %((2*j+2), (arrayBond[j][2]+1), (arrayBond[j][3]+1))
                k=j+1
            else:
                print >>f, "%d 1 %d %d" %((j+k+1), (arrayBond[j][1]+1), (arrayBond[j][2]+1))


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

        k=0
       #print array
        x, y, z = self.randomXYZarray(self.nSpheres)
        for ns in range(self.nSpheres):
            for j in range(len(array)):
                k+=1
                print >>f, "%d %d 1 %f %f %f" % ((k),ns+1,array[j][0]*radius+x[ns],array[j][1]*radius+y[ns],array[j][2]*radius+z[ns])
                


        self.sphereBeads = k
        solventBead = int((self.xhi-self.xlo)*(self.xhi-self.xlo)*(self.xhi-self.xlo)*self.rho - self.sphereBeads)
        x, y, z = self.randomXYZarray(solventBead)
          
        for p in range(0,solventBead):
            print >>f, "%d %d 2 %f %f %f" % (self.sphereBeads+p+1, self.nSpheres+p,x[p],y[p],z[p]) ##XXX this can not be 2 since this is the number of the molecule.
        k=0
        print >>f,"\n Bonds \n"
        for ns in range(self.nSpheres):
            part = len(array)*ns    ##This part value takes into account that there are more than one sphere, an shift by this value the global index of the particle in the bond topology construction
            for j in range(len(arrayBond)):
                k+=1
                print >>f, "%d 1 %d %d" %((k), (arrayBond[j][0]+1)+part, (arrayBond[j][1]+1)+part)
                

    def randomXYZarray(self, randomPoints):
        x = np.random.rand(randomPoints)*(self.xhi-self.xlo)+self.xlo
        y = np.random.rand(randomPoints)*(self.yhi-self.ylo)+self.ylo
        z = np.random.rand(randomPoints)*(self.zhi-self.zlo)+self.zlo

        return x, y, z
        
oFile = sys.argv[1]
data = dataFile(oFile, 197, 1920, 0, 2, 1, 1, -15.0, 15.0, -15.0, 15.0, -15.0, 15.0, [1,1], 3, 30)
data.header('vertices', 'faces')

