"""
04-15-2013
@author: Nicolas Moreno Ch.
e-mail: nicolas.morenochaparro@kaust.edu.sa



CalculateZ routine is based on the code presented in: http://www.timteatro.net/2010/09/29/velocity-autocorrelation-and-vibrational-spectrum-calculation/
It computes the velocity autocorrelation funtion of a system, defined as c(t)=<v(0) * v(t)>.  c(0)=<v^2>=d k_b T /m

input: input: Atom index, type, molecule and positions (dumpfile)
output: msd



"""
import numpy as np
import sys



typeAnalysis = 'traj'

class analysis():

    def __init__(self, data):
        self.name = ["time", "mean square displacement"]
        self.label = ['t','msd']
        self.outType = ['x','y']
        self.a, self.m, self.t, self.x, self.y, self.z, self.lbox = data
        self.vol = self.lbox[0][0]*self.lbox[0][1]*self.lbox[0][2]  #constant volume
        self.nofp = float(len(self.a[0]))  ##Assuming constant number of particles, however this can be easily modified for varying n of p.
        #self.id = self.file.split("/")[2]
        self.sigma = 1./2.5
        self.indepVar = 0;



    def run(self):
        coms = []#np.array([self.x, self.y, self.z]).T
        vaf = []
        time = []
        bin = 1.
        max = 2000
        M=len(self.a)
        velvel = 0
        i=0
        
        for m in range(M):
             xcom, ycom, zcom = self.comi(self.x[m]), self.comi(self.y[m]), self.comi(self.z[m])
             coms.append([xcom, ycom, zcom])
        
        for m in range(M): ## the index of this loop basically correspond to the dt, this is m = 3 means all the set of displacement for snapshots with difference of 3.
            C = 0.
            s = 0.
            for n in range(0,M-m): ## since the value at t0 in each loop correspond to the n time step the interval is reduced, so I cannot look for correlation 1 step forward for the last snapshot for example.
                                   
                deltax = self.filt(self.x[n+m], self.x[n])
                deltay = self.filt(self.y[n+m], self.y[n])
                deltaz = self.filt(self.z[n+m], self.z[n])
                
                ##Computing center of mass drifting
                comDrift = 0.#(X[n+m][3] - X[n][3])**2 + (X[n+m][4] - X[n][4])**2 + (X[n+m][5] - X[n][5])**2
                #print comDrift, 'com'
                ## looping over all the time steps looking forward only for correlation.
                #s += (deltax-(X[n+m][3] - X[n][3]))**2+(deltay-(X[n+m][4] - X[n][4]))**2+(deltaz-(X[n+m][5] - X[n][5]))**2
                s += (deltax)**2+(deltay)**2+(deltaz)**2 - comDrift     ## square distance --I NEED TO VERIFY THIS COMPUTATION. ARE THEY CORRECTLY ''SQUARED''?
            C = 1./self.nofp * s.sum()  ##It seems this was wrongly idented.
            C/=(M-m)
            vaf.append(C)


        pad = np.zeros(max*5)
        padTime = (np.array(time)+max)

        
        #vaf = np.append(vaf, pad)
        #time = np.append(time, padTime)

        time = np.arange(len(vaf))
        
        #vaf = self.smooth(vaf)
        
        #fft =ftb.fft3DBox(vaf)
        #print (freqs)

        #F = np.fft.fftn(vaf)
        #F = np.fft.fftshift( F )
        #F= (np.abs(F))**2
        #F=np.log(F)
        
        vv = np.argsort(np.array(vaf))
        tt= np.array(time)[vv]
        
        #freqs = (np.fft.fftfreq(fft.size, 1.))**2
        #plt.plot(time, vaf)
        #plt.plot(freqs, F)
        #plt.show()
        #self.plot(vaf, time, fft, freqs, M)

        return time, vaf

    def filt(self, X0, X1):
        deltax = np.abs(X1-X0)
        #print deltax
        
        xmean = deltax.mean()
        xstd = deltax.std()

        xindi = np.where(deltax<(xmean-xstd))  ##sup limit
        #print deltax[xindi] 
        deltax[xindi] = deltax[xindi]+self.lbox[0][0]
        xinds = np.where(deltax>xmean+xstd)  ## inf limit
        #print deltax[xinds]
        deltax[xinds] = -deltax[xinds]+self.lbox[0][0]  ##Assuming cubic box

        fixed = len(xindi[0])+len(xinds[0])
        #print len(deltax),fixed#, xindi
        return deltax#, fixed

    def smooth(self, vaf):
        M = len(vaf)
        sqtp = np.sqrt(2*np.pi)
        self.sigma*=M
        for i in range(M):
            vaf[i]*=np.exp(-i**2/(2*self.sigma**2))/(self.sigma*sqtp)
            #for m in range (1,M):
            #vaf/=vaf[0]
        return vaf

    def comi(self, x):
        #x  = np.arange(-lbox,lbox+1, dtype=float)
        #x  = np.arange(0,lbox+1, dtype=float)
        #x = np.random.rand(100)*100
        #x = np.array([16,17,18,19,20])
        xmax = x.max()
        map = 2.*np.pi/xmax
        teta = x*map
        s = np.sin(teta)
        c = np.cos(teta)
        ps = xmax*c
        xi = xmax*s
        psm = ps.mean()
        xim = xi.mean()

        #ps =ps-(psm) 
        #xi =xi-(xim)

        tetam =np.arctan2(xim,psm)#+np.pi  #If the interval goes from -lbox to lbox, it is not necesary to add pi.
        xcom = 1/map*tetam

        #tetaShifted = np.arctan2(xi,ps) - tetam
        #xShifted =  1/map*tetaShifted
        #print xcom
        #plt.plot(ps,xi)
        #plt.show()
        #print xShifted
        return xcom

    def test(self, folderID, first):
        folderID   = folderID #sys.argv[1]
        first = first #sys.argv[2]
        v = MSD(folderID, first)
        v.calculateZ()
