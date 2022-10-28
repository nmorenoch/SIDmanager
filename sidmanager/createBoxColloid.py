"""
08-30-2013

@author: Nicolas Moreno Chaparro
         nicolas.morenochaparro@oist.jp

This script generates a lammps data file with the coordinates of colloid strcutures. Initially constructed for Mobius colloids.

"""
import sys, os, re
import numpy as np



class dataFile():

   def __init__(self, oFile, lbox, m, rho, colloidNx, colloidK, ro, nColloids, eFact, nSolv, solvFrac, nTwists=0.5,location=-1,ks=50,simtype='ring'):
      #Variable location is used to place initially all the blocks in one side of the box.
      self.simtype = simtype
      self.f = oFile
      self.lbox = lbox
      self.m = m.split(',')
      self.rho = rho
      self.colloidNx = int(colloidNx)
      self.colloidK  = int(colloidK)
      #self.colloidNy = self.colloidK+1
      self.ro = ro
      self.nColloids = nColloids
      self.nTwists =nTwists

      self.eFact = eFact
      self.nSolv = nSolv
      print (nSolv, "number of solvent types")
      self.solvFrac = np.array(solvFrac.split(','), dtype=float)
      self.location = location
     
      self.nPart         = int(self.lbox**3*self.eFact*self.rho)
      self.parPerColloid = self.colloidNx*(self.colloidK+1)
      self.colloidBeads  = self.parPerColloid*self.nColloids

      self.ks=ks #Spring constant of the bonds, if 0 bonds are not created
      if self.ks>0:
        self.nBonds        = (self.colloidNx*(self.colloidK+1)+2*self.colloidK*self.colloidNx)*self.nColloids
      else:
          self.nBonds = 0
          
      self.nAngles      = 0 #self.nChains*(self.polyL-2)#+self.addNAngles
      self.solventBeads = self.nPart - self.colloidBeads


      self.partPerSolvent = self.solventBeads*self.solvFrac

      self.xlo, self.xhi = -self.lbox/2., self.lbox/2.
      self.ylo, self.yhi = -self.lbox/2., self.lbox/2.
      self.zlo, self.zhi = -self.lbox/2.*self.eFact, self.lbox/2.*self.eFact

      self.aType  = self.nSolv + 1 ##Only one type of particle for the colloid
      self.bType  = 1 #+ len(self.blockFracs) - 1   ##Assuming that not consecutive block definition are not connected (e.g: for types [1,2,3] 1-1, 2-2, 3-3, 1-2, 2-3)
      self.anType = 1 # len(self.blockFracs) + len(self.addBlockFracs) #+ len(self.blockFracs) - 1   ##Same previous assumption

      if self.aType != len(self.m):
         print ("Stypes:%s, Colloid:%s" % (self.nSolv ,1))
         print ('ERROR ATOM TYPES AND MASS DO NOT MATCH: %d, %d' %(self.aType, len(self.m)))
         exit()

      ##Initialize global counters for atom and molecule ID
      self.currentA = 1
      self.currentM = 1
      self.currentBond = 1
      self.currentColloid = 0

      ## Creating Header
      self.header()

   def header(self):
      f = open(self.f, 'w')
      print >>f,"# Lammps data file generator"
      print >>f,"# @Author: Nicolas Moreno"
      #print >>f,self.setup
      
      print >>f, " "
      print >>f, self.nPart.__str__()+" atoms"
      print >>f, (self.nBonds).__str__()+" bonds"
      print >>f, (self.nAngles).__str__()+" angles"
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
      print >>f, "0.0 0.0 0.0 xy xz yz"
      
      #Masses
      print >>f, " "
      print >>f,"Masses \n"
      for i in range(len(self.m)):
         print >>f, (i+1).__str__()+" "+self.m[i].__str__()
         
      print >>f,"\n Atoms"
      print >>f, " "
      f.close()

   def createMobius(self, f):
       
        h = np.sqrt(3)/2.*self.ro ##Height of the triangles dividing the strip
        X = self.colloidNx*self.ro

        px0 = np.arange(-X/2.,X/2.,self.ro)
        py0 = np.zeros(int(self.colloidNx))
        pz0 = np.zeros(int(self.colloidNx))

        bandX,bandY, bandZ  = px0, py0, pz0

        for i in np.arange(1,self.colloidK+1):
            if i%2==0:
                shiftX = 0
            else: shiftX = 1
            bandX = np.hstack([bandX,px0+shiftX*self.ro/2.])
            bandY = np.hstack([bandY,py0+i*h])
            bandZ = np.hstack([bandZ,pz0])

        bandY = bandY-(self.colloidK*h/2.)  ###Shift the band to the origin on Y
        exponRate = 1.
        rotRate = self.nTwists*np.pi/((X-self.ro/2.)/2.)**exponRate ##roation rate for the cylinder around X.
        
        allPoints = np.vstack([bandX,bandY,bandZ]).T
        rotBands = np.zeros_like(allPoints)
        R = (X-0/2.)/(2*np.pi)

        if self.nTwists == 0: ##Number of rotations of the strip
            allPoints = np.round(allPoints*self.rotate(np.pi/2,0),4)

        for i in range(0,len(allPoints)):
            rotBands[i] = np.round(allPoints[i]*self.rotate((allPoints[i][0])**exponRate*rotRate,0),4)#*i.T

           

            if self.simtype=='ring': ##If the ntwist variable is selected with a negative value helicoidal ribbon can be created.
                
                rotBands[i][0] = rotBands[i][0]-allPoints[i][0] ##return to the origin in x axis
                rotBands[i][1] = rotBands[i][1]+R  #Translate by Radius in the Y axis
                rotBands[i][2] = rotBands[i][2]
                
                angle = allPoints[i][0]*np.pi/(X/2.)  ## As it is closes the  cylinder ok Equivalent angle given the x position
                rotBands[i] = np.round(rotBands[i]*self.rotate(angle,2),4) ###Rotating around Z asis to fill the whole circunference.

        an1,an2,an3 = np.random.rand(3)*np.pi*2
        p1,p2,p3 = 0.,0.,0.
        if self.location ==-1: #If other flag specified put the band in the fraction of the box specified along the z axis
            p1,p2,p3 = np.random.rand(3)*(self.zhi-self.zlo)+self.zlo
        ##location in this case is the fraction of of the box size in the z direction where the particle is located.
        else:
            p3 = self.location*(self.zhi-self.zlo)+self.zlo
            
        rotBands = np.round(rotBands*self.rotate(an1,2)*self.rotate(an2,0)*self.rotate(an3,1),4)
        rotBands = rotBands+[p1,p2,p3] ##Translate randomly
        print ("random angles:%s, %s, %s" % (an1,an2,an3))

        for i in rotBands:
            f.write("%d %d %d %f %f %f\n" % (self.currentA, self.currentM, 1, i[0],i[1], i[2]))
            self.currentA+=1
        self.currentM+=1

   def wBond(self, f, idbond,p1,p2):
        IDshift = self.currentColloid*self.parPerColloid ##shifting in the ID according to the current saved bond data for colloids
        #print IDshift, self.currentColloid
        f.write("%s 1 %s %s \n"% (idbond, p1+IDshift,p2+IDshift))

   def defineBonds(self,fil):
        k = self.colloidK
        nx = self.colloidNx
        idbond = self.currentBond
        idlocal = 0
        firstColumn = False ##Flag for first column
       # print nx,k,self.currentBond

        for i in range(1,k+2):
            for j in range(1,nx+1):
                idlocal+=1
                
                if (idlocal-1)%nx ==0:##First Column
                    firstColumn=True
                if i%(k+1)!=0: ##Do it for all but the last row
                    self.wBond(fil, idbond, idlocal, idlocal+nx)
                    idbond+=1
                    #Conditional to avoid do it for last row
                    if j%nx==0: ##Last column of the band
                        self.wBond(fil,idbond, idlocal, idlocal-nx+1)
                        idbond+=1
                        if i%2==0:
                            self.wBond(fil,idbond, idlocal, idlocal+1)
                            idbond+=1
                        else:
                            self.wBond(fil,idbond, idlocal, idlocal+nx-1)
                            idbond+=1        
                    elif j%nx!=0: ##all  (nx-1) columns
                        self.wBond(fil,idbond, idlocal, idlocal+1)
                        idbond+=1
    
                        ##Assuming left-shifted coordinates for the first row of particles in the band
                        if i%2==0:  
                            self. wBond(fil,idbond, idlocal, idlocal+nx+1)
                            idbond+=1                   
                        elif i%2!=0:
                            if firstColumn==False:
                                self.wBond(fil,idbond, idlocal, idlocal+nx-1)
                                idbond+=1
                            else:
                                self.wBond(fil,idbond, idlocal, idlocal+2*nx-1)
                                idbond+=1
                                firstColumn = False  ##Reset flag
                else:
                    if j%nx==0: ##Last column of the band
                        self.wBond(fil,idbond, idlocal, idlocal-nx+1)
                        idbond+=1
                    elif j%nx!=0: ##all  (nx-1) columns
                        self.wBond(fil,idbond, idlocal, idlocal+1)
                        idbond+=1
        self.currentBond  = idbond
        self.currentColloid +=1   ##Adding new colloid to counter

   def atoms(self):
       f = open(self.f, 'a')
      
       # WRITTING ON FILE POLYMER BEADS COORDS
       for nC in range(self.nColloids):
           self.createMobius(f)

          
       # WRITTING ON FILE SOLVENT BEADS COORDS
       if self.location>2:
#          xc, yc, zc = self.randomXYZarrayAtLocation(self.solventBeads)
          xc, yc, zc = self.randomXYZarrayAtLevels(self.solventBeads)
       else:
          xc, yc, zc = self.randomXYZarray(self.solventBeads)
       a = self.currentA  ##Takes the last indezx used in polymer bead definition
       t = 2
       m = self.currentM

       k=0 #starting counter
       flagS = 0 ## this is a flag to check when the type of solvent must be changed. I prefer to do this in this way rather that define an array with the type of all the solvents since this data is ussualy the larger, so its better not to stored.
       pps = self.partPerSolvent[flagS]
       for i in range(self.solventBeads):
          f.write("%d %d %d %f %f %f\n" % (a+k, m+k, t+flagS, xc[i],yc[i], zc[i]))
          if pps <= k:
             flagS+=1
             pps += self.partPerSolvent[flagS]
          k+=1

       if self.nBonds>0:
            # WRITTING ON FILE BONDS
            k=0
            print >>f,"\n Bonds"
            print >>f, " "
            for i in range(self.nColloids):
                self.defineBonds(f)

        # WRITTING ON FILE ANGLES
       if self.nAngles>0:
            k=0
            print >>f,"\n Angles"
            print >>f, " "
            for i in range(self.nChains):
                for j in range(1,self.polyL-1):
                    f.write("%d %d %d %d %d\n" % (allia[k],allta[k], j+(i*self.polyL), j+1+(i*self.polyL), j+2+(i*self.polyL)))
                    k+=1
        
       f.close()


   def rotate(self,theta,axis):
        c, s = np.cos(theta), np.sin(theta)
        if axis ==0:
            R = np.matrix([[1.,0.,0.], [0,c, -s], [0,s, c]])
        elif axis == 2:
            R = np.matrix([[c, -s,0], [s, c,0],[0.,0.,1.]])
        elif axis == 1:
            R = np.matrix([[c, 0.,s], [0., 1.,0.],[-s,0.,c]])
        return R.T
        

   def randomPositionAtRoBiased(self, ro):
        x = np.random.rand(1)*(ro-ro/2)
        y = np.sqrt(np.random.rand(1)*(ro**2-x**2))
        z = np.sqrt(ro**2-x**2-y**2)
        return x,y,z

   def randomPositionAtRo(self, ro):
       alpha = np.random.rand(1)*2*np.pi
       beta = np.random.rand(1)*2*np.pi
       x = ro*np.sin(alpha)*np.cos(beta)
       y = ro*np.sin(alpha)*np.sin(beta)
       z = ro*np.cos(alpha)
       return x,y,z

   def randomXYZarray(self, randomPoints):
        print (randomPoints, 'randomPoints')
        x = np.random.rand(randomPoints)*(self.xhi-self.xlo)+self.xlo
        y = np.random.rand(randomPoints)*(self.yhi-self.ylo)+self.ylo
        z = np.random.rand(randomPoints)*(self.zhi-self.zlo)+self.zlo

        return x, y, z

   def randomXYZarrayAtLocation(self, randomPoints):
        print (randomPoints, 'randomPoints')
        z = self.polyLimits +np.random.rand(randomPoints)*(self.zhi-self.polyLimits)#(1-self.pFrac+self.addFrac)-self.polyLimits
                                                                                    #        x = np.random.rand(randomPoints)*(self.lbox*(1-self.pFrac+self.addFrac))-self.lbox/2
        print (self.polyLimits)
        y = np.random.rand(randomPoints)*(self.yhi-self.ylo)+self.ylo
        x = np.random.rand(randomPoints)*(self.xhi-self.xlo)+self.xlo

        return x, y, z  


   def randomXYZarrayAtLevels(self, randomPoints):
        print (randomPoints, 'randomPoints')
        beadS1 = int(randomPoints*self.solvFrac[0])
        beadS2 = randomPoints-beadS1

        x = np.random.rand(randomPoints)*(self.xhi-self.xlo)+self.xlo
        y = np.random.rand(randomPoints)*(self.yhi-self.ylo)+self.ylo

        frOfBox = self.pFrac+self.addFrac+self.solvFrac[0]*(1-self.pFrac+self.addFrac)
        z1 = np.random.rand(beadS1)*(self.lbox*self.eFact*frOfBox)-self.zhi
        z2 = self.zhi-np.random.rand(beadS2)*(self.lbox*self.eFact*(1-frOfBox))#(1-self.pFrac+self.addFrac)-self.polyLimits
                                                                                    #        x = np.random.rand(randomPoints)*(self.lbox*(1-self.pFrac+self.addFrac))-self.lbox/2
        print (self.polyLimits)
        #print x.shape, z2.shape, z1.shape, self.solvFrac[0], beadS1, beadS2,randomPoints-beadS1,'shaoes'
        if len(z2)>0:
            z = np.append(z1,z2,1)
            print ("%s Particles for second solvent included" % len(z2))
        else:
            z = z1
        #print x.shape, z.shape, z2.shape, z1.shape, self.solvFrac[0], beadS1, beadS2,randomPoints-beadS1,'shaoes'
        return x, y, z 




if sys.argv[0] == 'createBoxColloid.py':
   oFile = sys.argv[1]
   ro = 0.5
   L = 63
   W = 10 
   lbox = 20
    
   data =   dataFile(oFile, lbox, '1, 1', 3., L, W, ro, 1, 1., 1, "1.0, 0.0, 0.0", nTwists=0.0,location=1)
#   def __init__(self, oFile, lbox, m, rho, colloidNx, colloidK, ro, nColloids, eFact, nSolv, solvFrac, nTwists=0.5,location=0):
   data.atoms()
else: print ("createBox: Running from another script" )

   
"""
What JAVA module reads is:
self.init, sim['$P_l$'], sim['$P_{f,i}$'], sim['$D_poly$'], sim['$r_o$'], sim['$K_s$'], sim['$K_b$'], sim['$l_{box}$'], sim['$H_{frac}$'], sim['$\\rho$'], sim
                        ['T'], sim['SID'], sim['dFromSurf'], sim['isSurface'], sim['isCaps'], sim['dOsurf'], sim['nofSegments'], sim['surfBondL'],
                         sim['bodyPatch'], sim['$e_fact$'], sim['$N_{s}$'], sim['$S_{frac}^M$']

"""
