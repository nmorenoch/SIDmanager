# %load sidmanager/createBox.py
"""
06-24-2020

@author: Nicolas Moreno Chaparro
         nicolas.morenochaparro@oist.jp

Adapted from createBox script from DBC projects. Most of the functionalities still remain.

This script generates a lammps data file with the coordinates of a polymer system. Initially this code is only for linear polymer (diblock) constructions

"""
import sys, os, re
import numpy as np



class dataFile():

   def __init__(self, oFile, dimension, lbox, rho, polyL, ro, pFrac, blockFracs, nSolv, solvFrac, additBF=None, addFrac=0.0, addLen=1, location=0,lcore=1):
      #Variable location is used to place initially all the blocks in one side of the box.
      self.dimension = dimension
      self.lcore = lcore ### size of the core for fs simulations
      self.f = oFile
      self.lbox = lbox
      self.eFact = 1. ##Elongation factor is an old parameter to scale z direction, ommited so far so set in 1
      #self.m = m.split(',')
      self.rho = rho
      self.polyL = polyL
      self.ro = ro
      self.pFrac = pFrac
      self.blockFracs = np.array(blockFracs.split(','), dtype=float)
      self.nSolv = nSolv
      print (nSolv)
      self.solvFrac = np.array(solvFrac.split(','), dtype=float)
      self.location = location
      ##If DBC addtive used
      try:
         self.addBlockFracs =  np.array(additBF.split(','), dtype=float)
         print('lenght of additive %s' % addLen)
      except:
         self.addBlockFracs =  np.array([], dtype=float)
         print ("No DBC additive")
      self.addFrac = addFrac
      self.addLen = addLen

      self.nPart        = int(self.lbox**self.dimension*self.rho)
      
      #Calculated Values from input if additive
      self.addBeads = int(np.round(self.nPart*self.addFrac))
      self.addNChains      = int(self.addBeads/self.addLen)
      self.addNBonds       = self.addNChains*(self.addLen-1)
      self.addNAngles      = self.addNChains*(self.addLen-2)
      if  len(self.addBlockFracs)>0:
        self.addBlockLen = np.round(self.addLen*self.addBlockFracs).astype(int) ## Calculating the numebr of beads per block in additive
        self.addBlockLen[-1] = self.addLen-self.addBlockLen[:-1].sum() # and normalizing.

      #Calculated Values from input
    
      self.polymerBeads = int(np.round(self.nPart*self.pFrac))#+self.addBeads
      self.nChains      = int(self.polymerBeads/self.polyL)#+self.addNChains
      self.nBonds       = self.nChains*(self.polyL-1)#+self.addNBonds
      self.nAngles      = self.nChains*(self.polyL-2)#+self.addNAngles
      if self.nAngles<0:
        self.nAngles = 0
      
      self.solventBeads = self.nPart - (self.nChains*self.polyL) - (self.addNChains*self.addLen)
      print (self.addNChains*self.addLen), "number of add beads"

      self.blockLen = np.round(self.polyL*self.blockFracs).astype(int) ## Calculating the numebr of beads per block 
      self.blockLen[-1] = self.polyL-self.blockLen[:-1].sum() # and normalizing.

      self.partPerSolvent = self.solventBeads*self.solvFrac

      if self.location==5: ###location of the whole box centered in 0
          self.xlo, self.xhi = -self.lbox/2., self.lbox/2.
          self.ylo, self.yhi = -self.lbox/2., self.lbox/2.
          self.zlo, self.zhi = -self.lbox/2.*self.eFact*(self.dimension-2), self.lbox/2.*self.eFact*(self.dimension-2)
      else:
          self.xlo, self.xhi = 0, self.lbox
          self.ylo, self.yhi = 0, self.lbox
          self.zlo, self.zhi = 0, self.lbox*self.eFact*(self.dimension-2)

      self.aType  = self.nSolv + len(self.blockFracs) + len(self.addBlockFracs)
      self.bType  = len(self.blockFracs) + len(self.addBlockFracs) #+ len(self.blockFracs) - 1   ##Assuming that not consecutive block definition are not connected (e.g: for types [1,2,3] 1-1, 2-2, 3-3, 1-2, 2-3)
      self.anType = len(self.blockFracs) + len(self.addBlockFracs) #+ len(self.blockFracs) - 1   ##Same previous assumption

      #if self.aType != len(self.m):
      #   print ("Stypes:%s, BlockTypes:%s, addBlockTypes:%s" % (self.nSolv , len(self.blockFracs) , len(self.addBlockFracs)))
      #   print ('ERROR ATOM TYPES AND MASS DO NOT MATCH: %d, %d' %(self.aType, len(self.m)))
      #   exit()

      self.header()

      ##If polymer
      if self.location==1:
         self.polyLimits = self.zlo+self.lbox*self.eFact*(self.solvFrac[0])#+self.pFrac+self.addFrac)
      else:
         self.polyLimits = self.zlo
      #self.setup = "# Polymer length: %d  -   Bond distance: %f  - Particle N. Density: %f   -  Fract. polymer: %f \n #solventFracts:  %f-  #nsolvents: %d" % (self.polyL, self.ro, self.rho, self.polyFrac, self.solvFrac[0], self.nSolv)


   def header(self):
      f = open(self.f, 'w')
      f.write("# Lammps data file generator \n")
      f.write("# @Author: Nicolas Moreno \n")
      #f.writeself.setup
      
      f.write("  \n")
      f.write(self.nPart.__str__()+" atoms \n")
      f.write((self.nBonds+self.addNBonds).__str__()+" bonds \n")
      f.write((self.nAngles+self.addNAngles).__str__()+" angles \n")
      f.write("0 dihedrals \n")
      f.write("0 impropers \n")
          
      f.write("  \n")
      f.write(self.aType.__str__()+" atom types \n")
      f.write(self.bType.__str__()+" bond types \n")
      f.write(self.anType.__str__()+" angle types \n")
      f.write("0 dihedral types \n")
      f.write("0 improper types \n")

      #Box limits
      f.write("  \n")
      f.write(self.xlo.__str__()+" "+self.xhi.__str__()+" xlo xhi \n")
      f.write(self.ylo.__str__()+" "+self.yhi.__str__()+" ylo yhi \n")
      
      if self.dimension==2:
        f.write("0.0 0.2 zlo zhi \n")
      else:
        f.write(self.zlo.__str__()+" "+self.zhi.__str__()+" zlo zhi \n")
      
      #Masses
      #f.write("  \n")
      #f.write("Masses \n \n")
      #for i in range(len(self.m)):
      #   f.write((i+1).__str__()+" "+self.m[i].__str__()+" \n")
         
      f.write("\n Atoms")
      f.write("  \n \n")
      f.close()
         

   def additive(self, f, preK):

       addNOfBl = len(self.addBlockLen) ## If additive > 0
       
       allp, allt, allm, allmb, allma, alltb, allta = [], [], [], [], [], [], []

       typesat = np.arange(3,addNOfBl+3)
       typesm = np.arange(self.nChains+1,self.nChains+self.addNChains+1)
       typesb = np.arange(3,addNOfBl+3)
       typesa = np.arange(3,addNOfBl+3)

       alli = np.arange(preK+1, preK+1+self.addNChains*self.addLen)
       allib = np.arange(preK+1,1+self.addNBonds)
       allia = np.arange(preK+1,1+self.addNAngles)

       if self.location==4:
          seedx  = np.random.rand(self.addNChains,1)*(self.lbox*(self.pFrac+self.addFrac))-self.lbox/2
          seedy = np.random.rand(self.addNChains,2)*self.lbox-self.lbox/2
          #seedz = np.random.rand(self.addNChains,1)*self.lbox-self.lbox/2

          seed = np.append(seedy,seedx,1) ##The first row correspond to the x coordinates
       elif self.location==6: ##Polymer localized at the center of the box in a cube equivalent to it fraction.
          seed = (np.random.rand(self.addNChains,3)-0.5)*(self.lcore)+self.lbox/2
       elif self.location==7: ##Polymer localized at the center of the box in a cube equivalent to it fraction.
          seed = (np.random.rand(self.addNChains,3)-0.5)*((self.lbox*self.eFact*(self.pFrac+self.addFrac))**(1./self.dimension))+self.lbox/2
       else:
          seed = np.random.rand(self.addNChains,3)*self.lbox-self.lbox/2
       
      
       for k in range(len(typesat)):
          allt += [typesat[k]]*self.addBlockLen[k]
       allt=allt*self.addNChains
       
       for k in range(len(typesm)):
          allm += [typesm[k]]*self.addLen
          allmb += [typesm[k]]*(self.addLen-1)
          allma += [typesm[k]]*(self.addLen-2)

       for k in range(len(typesb)-1):
          alltb += [typesb[k]]*(self.addBlockLen[k])  ##Initial blocks will take te id of bonds at the interface
       alltb += [typesb[k+1]]*(self.addBlockLen[k+1]-1) ##In the last block one one is substracted since the las beads is not connected further.

       for k in range(len(typesa)-1):
          allta += [typesa[k]]*(self.addBlockLen[k])  ##Initial blocks will take te id of bonds at the interface
       allta += [typesa[k+1]]*(self.addBlockLen[k+1]-2) ##In the last block one one is substracted since the las beads is not connected further.

       allta=allta*self.addNChains
       alltb=alltb*self.addNChains

       # WRITTING ON FILE POLYMER additive BEADS COORDS
       k=0
       for i in seed:
          xc, yc, zc =  i[0],i[1],i[2]*(self.dimension-2)  #Current seed position (initial bead per chain)
          for j in range(self.addLen):
             xr, yr, zr = self.randomPositionAtRo(self.ro)  
             xc += xr
             yc += yr
             zc += zr    
             f.write("%d %d %f %f %f %d 1.0 0.0 125.0\n" % (alli[k], allt[k], xc,yc, zc, allm[k]))  ###Swap molecule to the end for hybrid
             k+=1

      
   def atoms(self):
       f = open(self.f, 'a')
       #addNOfBl = len(self.addBlockLen) ## If additive > 0
       nOfBl = len(self.blockLen) 
       
       allp, allt, allm, allmb, allma, alltb, allta = [], [], [], [], [], [], []

       typesat = np.arange(1,nOfBl+1)
       typesm = np.arange(1,self.nChains+1)
       typesb = np.arange(1,nOfBl+1)
       typesa = np.arange(1,nOfBl+1)

       alli = np.arange(1, self.nChains*self.polyL+1)
       allib = np.arange(1,self.nBonds+1)
       allia = np.arange(1,self.nAngles+1)

       
       if self.location==1: ##Polymer localized with one fluid 
          seedx  = np.random.rand(self.nChains,1)*(self.lbox*self.eFact*(self.pFrac+self.addFrac))-self.zhi
          seedy = np.random.rand(self.nChains,2)*self.lbox-self.lbox/2
          #seedz = np.random.rand(self.nChains,1)*self.lbox-self.lbox/2

          seed = np.append(seedy,seedx,1) ##The first row correspond to the x coordinates

       elif self.location==2:   ##Polymer localized within the whole domain.
          #          seed = np.random.rand(self.nChains,3)*self.lbox-self.lbox/2
          print ("located in all the box")
          seedz  = np.random.rand(self.nChains,1)*(self.lbox*self.eFact)-self.zhi
          seedxy = np.random.rand(self.nChains,2)*self.lbox-self.lbox/2
          #seedz = np.random.rand(self.nChains,1)*self.lbox-self.lbox/2

          seed = np.append(seedxy,seedz,1) ##The first row correspond to the x coordinates

       elif self.location==3: ##Polymer within one of the solvent when modeling two inmicible fluids
          seedx  = np.random.rand(self.nChains,1)*(self.lbox*self.eFact*(self.pFrac+self.addFrac+self.solvFrac[0]*(1-self.pFrac+self.addFrac)))-self.zhi
          seedy = np.random.rand(self.nChains,2)*self.lbox-self.lbox/2
          #seedz = np.random.rand(self.nChains,1)*self.lbox-self.lbox/2

          seed = np.append(seedy,seedx,1) ##The first row correspond to the x coordinates

       elif self.location==4: ##Polymer localized at the 0,0,0 in a cube equivalent to it fraction.
          seed = np.random.rand(self.nChains,3)*(self.lbox*self.eFact*(self.pFrac+self.addFrac)**(1./3.))-(self.lbox*self.eFact*(self.pFrac+self.addFrac)**(1./3.))/2
       elif self.location==6: ##Polymer localized at the center of the box in a cube equivalent to it fraction.
          seed = (np.random.rand(self.nChains,3)-0.5)*(self.lcore)+self.lbox/2
       elif self.location==7: ##Polymer localized at the center of the box in a cube equivalent to it fraction.
          seed = (np.random.rand(self.nChains,3)-0.5)*((self.lbox*self.eFact*(self.pFrac+self.addFrac))**(1./self.dimension))+self.lbox/2

       elif self.location>4: ##Polymer localized at the center of the box in a flat slabequivalent to it fraction.
          print ("located slab")
          lcube = self.lbox*self.eFact*(self.pFrac+self.addFrac)**(1./self.dimension)
          seedz  = lcube*(np.random.rand(self.nChains,1)-1/2.)*1./self.location
          seedxy = (lcube*(self.location)**0.5)*(np.random.rand(self.nChains,2)-0.5)#*self.lbox-self.lbox/2
          #seedz = np.random.rand(self.nChains,1)*self.lbox-self.lbox/2

          seed = np.append(seedxy,seedz,1) ##The first row correspond to the x coordinates
       elif self.location==5: ### Ccenter on the box at 0
          seed = np.random.rand(self.nChains,3)*self.lbox-self.lbox/2
       else: ###from 0 to lbox thi sis standar for latest sdpd however should be updated
          seed = np.random.rand(self.nChains,3)*self.lbox
      
       for k in range(len(typesat)):
          allt += [typesat[k]]*self.blockLen[k]
       allt=allt*self.nChains
       
       for k in range(len(typesm)):
          allm += [typesm[k]]*self.polyL
          allmb += [typesm[k]]*(self.polyL-1)
          allma += [typesm[k]]*(self.polyL-2)

       if self.nBonds>0: 
           for k in range(len(typesb)-1):
               alltb += [typesb[k]]*(self.blockLen[k])  ##Initial blocks will take te id of bonds at the interface
           alltb += [typesb[k+1]]*(self.blockLen[k+1]-1) ##In the last block one one is substracted since the las beads is not connected further.
       if self.nAngles:
           for k in range(len(typesa)-1):
               allta += [typesa[k]]*(self.blockLen[k])  ##Initial blocks will take te id of bonds at the interface
           allta += [typesa[k+1]]*(self.blockLen[k+1]-2) ##In the last block one one is substracted since the las beads is not connected further.

       allta=allta*self.nChains
       alltb=alltb*self.nChains

       # WRITTING ON FILE POLYMER BEADS COORDS
       k=0
       #print alli
       print (len(seed), seed.shape)
       for i in seed:
          xc, yc, zc =  i[0],i[1],i[2]*(self.dimension-2)  #Current seed position (initial bead per chain)
          for j in range(self.polyL):
             xr, yr, zr = self.randomPositionAtRo(self.ro)  
             xc += xr
             yc += yr
             zc += zr
             #print k
             f.write("%d %d %f %f %f %d 0.0 0.0 125.0\n" % (alli[k],allt[k], xc,yc, zc, allm[k]))
             k+=1

       if len(self.addBlockFracs) > 0:
          print ("writing additive")
          self.additive(f, k)
          
       # WRITTING ON FILE SOLVENT BEADS COORDS
       if self.location>0 :# or self.location<4:
#          xc, yc, zc = self.randomXYZarrayAtLocation(self.solventBeads)
          xc, yc, zc = self.randomXYZarrayAtLevels(self.solventBeads)
       else:
          xc, yc, zc = self.randomXYZarray(self.solventBeads)
       a = alli[k-1]+1+self.addNChains*self.addLen  ##Takes the last indezx used in polymer bead definition
       t = allt[k-1]+1+len(self.addBlockFracs)
       m = allm[k-1]+1+self.addNChains

       k=0 #starting counter
       flagS = 0 ## this is a flag to check when the type of solvent must be changed. I prefer to do this in this way rather that define an array with the type of all the solvents since this data is ussualy the larger, so its better not to stored.
       pps = self.partPerSolvent[flagS]
       for i in range(self.solventBeads):
          f.write("%d %d %f %f %f %d 0.0 0.0 125.0\n" % (a+k, t+flagS, xc[i],yc[i], zc[i], m+k))
          if pps <= k:
             flagS+=1
             pps += self.partPerSolvent[flagS]
          k+=1
          
        # WRITTING ON FILE BONDS
       k=0
       if self.nBonds>0:
           f.write(" \n")
           f.write("Bonds \n")
           f.write(" \n")
       for i in range(self.nChains):
          for j in range(1,self.polyL): ##Since counter starts in 1, there is no need to substract the value of pl
             f.write("%d %d %d %d\n" % (allib[k],alltb[k], j+(i*self.polyL), j+1+(i*self.polyL)))
             k+=1

       addK=1
       for i in range(self.addNChains):
          for j in range(1,self.addLen): ##Since counter starts in 1, there is no need to substract the value of pl
             f.write("%d %d %d %d\n" % (allib[k-1]+addK,alltb[k-1]+1, j+(i*self.addLen)+self.nChains*self.polyL, j+1+(i*self.addLen)+self.nChains*self.polyL))
             addK+=1

       if self.nAngles>0: 
           # WRITTING ON FILE ANGLES
           k=0
           f.write("\n Angles \n")
           f.write("\n")
           for i in range(self.nChains):
               for j in range(1,self.polyL-1):
                   f.write("%d %d %d %d %d\n" % (allia[k],allta[k], j+(i*self.polyL), j+1+(i*self.polyL), j+2+(i*self.polyL)))
                   k+=1

               addK=1
           for i in range(self.addNChains):
               for j in range(1,self.addLen-1):
                   f.write("%d %d %d %d %d\n" % (allia[k-1]+addK,allta[k-1]+1, j+(i*self.addLen)+self.nChains*self.polyL, j+1+(i*self.addLen)+self.nChains*self.polyL, j+2+(i*self.addLen)+self.nChains*self.polyL))
               addK+=1



       print(self.nChains)
       f.close()

   def randomPositionAtRoBiased(self, ro):
        x = np.random.rand(1)*(ro-ro/2)
        y = np.sqrt(np.random.rand(1)*(ro**2-x**2))
        z = np.sqrt(ro**2-x**2-y**2)
        return x,y,z

   def randomPositionAtRo(self, ro):
        if self.dimension==3:
           alpha = np.random.rand(1)*2*np.pi
        else:
            alpha =np.pi/2  ###zeroing z coord
        beta = np.random.rand(1)*2*np.pi
        x = ro*np.sin(alpha)*np.cos(beta)
        y = ro*np.sin(alpha)*np.sin(beta)
        z = ro*np.cos(alpha)
  
        return x,y,z

   def randomXYZarray(self, randomPoints):
        print (randomPoints, 'randomPoints')
        x = np.random.rand(randomPoints)*(self.xhi-self.xlo)+self.xlo
        y = np.random.rand(randomPoints)*(self.yhi-self.ylo)+self.ylo
        z = (self.dimension-2)*(np.random.rand(randomPoints)*(self.zhi-self.zlo)+self.zlo)

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


        lcube = self.lbox*self.eFact*(self.pFrac+self.addFrac)**(1./3.)
        gapz  = lcube*1./self.location
        
        z1 = np.random.rand(beadS1)*(self.lbox*self.eFact*frOfBox)-(self.zhi+gapz)
        z2 = gapz+self.zhi-np.random.rand(beadS2)*(self.lbox*self.eFact*(1-frOfBox))#(1-self.pFrac+self.addFrac)-self.polyLimits
                                                                                    #        x = np.random.rand(randomPoints)*(self.lbox*(1-self.pFrac+self.addFrac))-self.lbox/2
        print (self.polyLimits)
        print (x.shape, z2.shape, z1.shape, self.solvFrac[0], beadS1, beadS2,randomPoints-beadS1,'shaoes')
        if len(z2)>0:
            z = np.append(z1,z2,0)*(self.dimension-2)
            print ("%s Particles for second solvent included" % len(z2))
        else:
            z = z1*(self.dimension-2)
        
        #print x.shape, z.shape, z2.shape, z1.shape, self.solvFrac[0], beadS1, beadS2,randomPoints-beadS1,'shaoes'
        return x, y, z 

if sys.argv[0] == 'createBox.py':
    oFile = sys.argv[1]
    dim = 2.0
    rho = np.round(1/0.2**dim,0)
    mass = 0.2**dim
    Lbox = 12+1.6*2+4*2
    #(oFile, dimension, lbox, m, rho, polyL, ro, pFrac, blockFracs, nSolv, solvFrac, additBF="0", addFrac=0.0, addLen=1.0, location=0)
    data = dataFile(oFile,dim, Lbox, '%s, %s, %s'%(mass,mass,mass), rho, 2, 0.2, 1.0, '0.5, 0.5', 1, '1.0, 0.0, 0.0', location=6)

    data.atoms()
else: print("createBox: Running from another script") 
