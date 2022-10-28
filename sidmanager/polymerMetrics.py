"""

@author: nmorenoc
         nicolas.morenochaparro@kasut.edu.sa

Generates      with different metrics for the polymer chains, such as Rf, contour, etc


"""


import numpy as np
import sys, os, time 
import readDump as rd
import pylab as plt
import plotAll as PA


class polymerMetrics():

   def __init__(self,data,bin,lbox,plen):
      self.bin = bin
      self.e = np.array([bin,bin,bin])
      self.lbox = lbox
      self.a, self.m, self.t, self.x, self.y, self.z = data
      self.vol = self.lbox[0]*self.lbox[1]*self.lbox[2]
      self.O = np.array(0)
      self.midBox = np.array(self.lbox)/2.0
      deltaMol = self.m[1:] - self.m[:-1]
      #print(deltaMol)
      self.plen = plen
      self.molchange = np.append(0,1+np.where(deltaMol)[0])     ##appending 0 at the begining in order to star the index counting from 0
      

   def oldinit(self, file, bin,cols,lbox):
        
      self.file = file
      self.bin = bin
      self.e = np.array([bin,bin,bin])
      self.lbox = lbox
      self.a, self.m, self.t, self.x, self.y, self.z = rd.readDumpFile(self.file,cols=cols)
      self.vol = self.lbox[0]*self.lbox[1]*self.lbox[2]
      self.O = np.array(0)
      self.midBox = self.lbox[0]/2.0
      index = np.argsort(self.a)
      self.a = self.a[index]
      self.m = self.m[index]
      self.t = self.t[index]
      self.x = self.x[index]
      self.y = self.y[index]
      self.z = self.z[index]


      deltaMol = (self.m[1:] - self.m[:-1])
      self.molchange = np.append(0,np.where(deltaMol)[0])     ##appending 0 at the begining in order to star the index counting from 0
      
      
   def distance(self, x1,x2,y1,y2,z1,z2):
      dx = np.abs(x1-x2)
      dy = np.abs(y1-y2)
      dz = np.abs(z1-z2)
 
      if dx>self.midBox[0]:
          dx=np.abs(dx-self.lbox[0])
      if dy>self.midBox[1]:
          dy=np.abs(dy-self.lbox[1])
      if dz>self.midBox[2]:
          dz=np.abs(dz-self.lbox[2])
         
          
         #lDist = dist[dist>self.midBox]
         #lDist = np.abs(lDist-self.lbox[0])
         #dist[np.where(dist>self.midBox)[0]] = lDist
      dist = np.sqrt(dx**2+dy**2+dz**2)
      return dist

   def distanceBlockData(self, x1,x2,y1,y2,z1,z2):
      dx = np.abs(x1-x2)
      dy = np.abs(y1-y2)
      dz = np.abs(z1-z2)
 
      if dx>0:
          
         ldx = dx[dx>self.midBox[0]]
         ldx = np.abs(ldx-self.lbox[0])
         dx[np.where(dx>self.midBox[0])[0]] = ldx
         
         ldy = dy[dy>self.midBox[1]]
         ldy = np.abs(ldy-self.lbox[1])
         dy[np.where(dy>self.midBox[1])[0]] =  ldy
         
         ldz = dz[dz>self.midBox[2]]
         ldz = np.abs(ldz-self.lbox[2])
         dz[np.where(dz>self.midBox[2])[0]] =  ldz
         
          
         #lDist = dist[dist>self.midBox]
         #lDist = np.abs(lDist-self.lbox[0])
         #dist[np.where(dist>self.midBox)[0]] = lDist
      dist = np.sqrt(dx**2+dy**2+dz**2)
      return dist


   def angleCalc(self, p1,p2,p3):
      a = p1-p2
      b = p1-p3
      A = np.linalg.norm(a)
      B = np.linalg.norm(b)
      theta = np.arccos(a.dot(b)/(A*B))
      return theta
      
   
   def rf(self):
      rf = [] 
      theta = []
      cm = []  # center of mass of polymerchain
      #dimMol = self.molchange[1]   ## this number indicates the first position where there is a change in the type of molecule which is related with the molecule size, this will be used in the rf calculation assuming the same size for all the molecs.
     # print(self.molchange)

      for i in range (1,len(self.molchange)):

         dimMol = self.molchange[i]-self.molchange[i-1]   ## this number indicates the first position where there is a change in the type of molecule which is related with the molecule size, this will be used in the rf calculation assuming the same size for all the molecs.
         if self.plen==dimMol:
            x1 = self.x[self.molchange[i]-1]
            x2 = self.x[self.molchange[i]-dimMol]
            xcm = np.mean(self.x[self.molchange[i]-dimMol:self.molchange[i]-1])

            y1 = self.y[self.molchange[i]-1]
            y2 = self.y[self.molchange[i]-dimMol]
            ycm = np.mean(self.y[self.molchange[i]-dimMol:self.molchange[i]-1])
            
            z1 = self.z[self.molchange[i]-1]
            z2 = self.z[self.molchange[i]-dimMol]
            zcm = np.mean(self.z[self.molchange[i]-dimMol:self.molchange[i]-1])
            #print (x1, x2)
            #print (self.distance(x1,x2,y1,y2,z1,z2))
            rf.append(self.distance(x1,x2,y1,y2,z1,z2))
            cm.append([xcm,ycm,zcm])

            ##Taking advantage that the last value of the 'i' and therefore the molecchange index, teh rf for the last molec is calculated
         #   x1 = self.x[self.molchange[i]+1]
         #   x2 = self.x[self.molchange[i]+1+dimMol]

         #   y1 = self.y[self.molchange[i]+1]
         #   y2 = self.y[self.molchange[i]+1+dimMol]

          #  z1 = self.z[self.molchange[i]+1]
          #  z2 = self.z[self.molchange[i]+1+dimMol]
         
            p1 = np.array([x1,y1,z1]).transpose()  ##It is necessary to transpose in order to iterate in the next for loop over the realted p1,p2,p3 points
            p2 = np.array([x2,y2,z2]).transpose()
            p3 = np.array([x2,y1,z1]).transpose()
            
            theta.append(self.angleCalc(p1,p2,p3)) #With this for loop only the coordinates for each three points realted are sent together.  At the end of this loop a list with all the angles in a molecules is produced.
         #else:
         #       print(dimMol)
      return rf, theta, cm

   def lc(self):
      lc = []
      dimMol = self.molchange[1]
      for j in range (1,len(self.molchange)):
         x1 = self.x[self.molchange[j]-dimMol:self.molchange[j]]    #array[0....m-1]
         x2 = self.x[self.molchange[j]-dimMol+1:self.molchange[j]+1]    #array[1....m]

         y1 = self.y[self.molchange[j]-dimMol:self.molchange[j]]    #array[0....m-1]
         y2 = self.y[self.molchange[j]-dimMol+1:self.molchange[j]+1]    #array[1....m]

         z1 = self.z[self.molchange[j]-dimMol:self.molchange[j]]    #array[0....m-1]
         z2 = self.z[self.molchange[j]-dimMol+1:self.molchange[j]+1]    #array[1....m]

         l =self.distance(x1,x2,y1,y2,z1,z2) 
         lc.append(l.sum())  #Since the x,y,z are actually vectors of coordinates, l is a vector of calculated distance between the beads in a molecule, then lc should be stimated as the sumation of this values

         ##Taking advantage that the last value of the 'i' and therefore the molecchange index, the lc for the last molec is calculated

      x1 = self.x[self.molchange[j]+1:self.molchange[j]+dimMol+1]    #array[0....m-1]
      x2 = self.x[self.molchange[j]+2:self.molchange[j]+dimMol+2]    #array[1....m]

      y1 = self.y[self.molchange[j]+1:self.molchange[j]+dimMol+1]    #array[0....m-1]
      y2 = self.y[self.molchange[j]+2:self.molchange[j]+dimMol+2]    #array[1....m]

      z1 = self.z[self.molchange[j]+1:self.molchange[j]+dimMol+1]    #array[0....m-1]
      z2 = self.z[self.molchange[j]+2:self.molchange[j]+dimMol+2]    #array[1....m]

      l =self.distance(x1,x2,y1,y2,z1,z2) 
      lc.append(l.sum())   #Since the x,y,z are actually vectors of coordinates, l is a vector of calculated distance between the beads in a molecule, then lc should be stimated as the sumation of this values
      
      return lc

   

   def ang(self):
      ang = []
      theta = []
      dimMol = self.molchange[1]
      for j in range (1,len(self.molchange)):
         x1 = self.x[self.molchange[j]-dimMol:self.molchange[j]-1]    #array[0....m-2]
         x2 = self.x[self.molchange[j]-dimMol+1:self.molchange[j]]    #array[1....m-1]
         x3 = self.x[self.molchange[j]-dimMol+2:self.molchange[j]+1]    #array[2....m]

         y1 = self.y[self.molchange[j]-dimMol:self.molchange[j]-1]    #array[0....m-2]
         y2 = self.y[self.molchange[j]-dimMol+1:self.molchange[j]]    #array[1....m-1]
         y3 = self.y[self.molchange[j]-dimMol+2:self.molchange[j]+1]    #array[2....m]

         z1 = self.z[self.molchange[j]-dimMol:self.molchange[j]-1]    #array[0....m-2]
         z2 = self.z[self.molchange[j]-dimMol+1:self.molchange[j]]    #array[1....m-1]
         z3 = self.z[self.molchange[j]-dimMol+2:self.molchange[j]+1]    #array[2....m]

         p1 = np.array([x1,y1,z1]).transpose()        ##It is necessary to tranpose in order to iterate in the next for loop over the realted p1,p2,p3 points
         p2 = np.array([x2,y2,z2]).transpose()
         p3 = np.array([x3,y3,z3]).transpose()

         # print self.molchange
         #print self.x[0:18]
         for i in range(len(x1)):
            #print self.angleCalc(p1[i],p2[i],p3[i])
            theta.append(self.angleCalc(p1[i],p2[i],p3[i])) #With this for loop only the coordinates for each three points realted are sent together.  At the end of this loop a list with all the angles in a molecules is produced.

         ang.append(np.array(theta).sum()/len(theta))
         ##Taking advantage that the last value of the 'i' and therefore the molecchange index, the lc for the last molec is calculated

      x1 = self.x[self.molchange[j]+1:self.molchange[j]+dimMol]    #array[0....m-2]
      x2 = self.x[self.molchange[j]+2:self.molchange[j]+dimMol+1]    #array[1....m-1]
      x3 = self.x[self.molchange[j]+3:self.molchange[j]+dimMol+2]    #array[2....m]
      
      y1 = self.y[self.molchange[j]+1:self.molchange[j]+dimMol]    #array[0....m-2]
      y2 = self.y[self.molchange[j]+2:self.molchange[j]+dimMol+1]    #array[1....m-1]
      y3 = self.y[self.molchange[j]+3:self.molchange[j]+dimMol+2]    #array[2....m]
      
      z1 = self.z[self.molchange[j]+1:self.molchange[j]+dimMol]    #array[0....m-2]
      z2 = self.z[self.molchange[j]+2:self.molchange[j]+dimMol+1]    #array[1....m-1]
      z3 = self.z[self.molchange[j]+3:self.molchange[j]+dimMol+2]    #array[2....m]
      
      p1 = np.array([x1,y1,z1]).transpose()  ##It is necessary to transpose in order to iterate in the next for loop over the realted p1,p2,p3 points
      p2 = np.array([x2,y2,z2]).transpose()
      p3 = np.array([x3,y3,z3]).transpose()
      
      for i in range(len(x1)):
         theta.append(self.angleCalc(p1[i],p2[i],p3[i])) #With this for loop only the coordinates for each three points are sent together.
      ang.append(np.array(theta).sum()/len(theta))
      return ang
      
   def averages(self, rf, lc, ang):
      rfAv  = np.array(rf).sum()/len(rf)
      lcAv  = np.array(lc).sum()/len(lc)   
      angAv = np.array(ang).sum()/len(ang)

      return rfAv,lcAv, angAv


   def test(self):
      print("This is a test of polymerMetrics.py")
      #pm = self.polymerMetrics(self.file,self.bin)

      #p = pMet(a, m, t, x, y, z, lbox, folderID)

      rf = self.rf()
      lc = self.lc()
      ang = self.ang()
      self.averages(rf,lc,ang)

      ind = "201_noShift"
      binHist = 5
      units = "nm"
      ylabel = "distribution"
      PA.plotHistogram(rf, binHist, "End to end distance ($R_f$) "+units, ylabel, "rfHistogram"+ind, False)
      PA.plotHistogram(lc, binHist, "Contour length ($L_c$) "+units, ylabel, "lcHistogram"+ind, False)
      PA.plotHistogram(ang, binHist, "Angles ($R_f$) "+units, ylabel, "angleHistogram"+ind, False)

      print ('End test')

#file = sys.argv[1]
#p = polymerMetrics(file,1)
#p.test()
