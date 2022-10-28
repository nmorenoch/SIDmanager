"""
04-13-2013
@author: Nicolas Moreno Ch.
e-mail: nicolas.morenochaparro@kaust.edu.sa



CalculateZ routine is based on the code presented in: http://www.timteatro.net/2010/09/29/velocity-autocorrelation-and-vibrational-spectrum-calculation/
It computes the velocity autocorrelation funtion of a system, defined as c(t)=<v(0) * v(t)>.  c(0)=<v^2>=d k_b T /m

input: Array of dump velocity files to compute the time correlation function



"""
import numpy as np
import sys

typeAnalysis = "traj"

class analysis():

    def __init__(self, data):
        self.name = ["time", "velocity autocorrelation function"]
        self.label = ['t', '$\left\langle V_tV_0 \right\rangle$']
        self.outType = ['x','y']
        self.a, self.m, self.t, self.vx, self.vy, self.vz, self.lbox = data
        self.vol = self.lbox[0][0]*self.lbox[0][1]*self.lbox[0][2]
        self.nofp = float(len(self.a[0]))
        self.sigma = 1./2.5



    def run(self):
        vaf = []
        time = []
        bin = 1.
        max = 1000
        M=len(self.a)
        velvel = 0
        i=0
        #points =0
        #for dumpFile in sorted(glob.glob(os.path.join(self.folderID+'/vel/', '*.dump')), key=os.path.getmtime):
        #    i+=1.
        #    r = i/bin
        #    if (int(r)-r) == 0 and points<max:
        #        points+=1
                #print dumpFile
        #         a,m, t, vx, vy, vz, lbox = rd.read(dumpFile)
        #        index = np.argsort(a)   ##sorted index according to the ID name fro each particle
        #        vx = vx[index]
        #        vy = vy[index]
        #        vz = vz[index]
        #        V.append([vx,vy,vz])
        #M=points
        for m in range(M):
            C = 0.
            for n in range(M-m-1):
                ## looping over all the time steps looking forward only for correlation.
                C += self.vx[n+m]*self.vx[n]+self.vy[n+m]*self.vy[n]+self.vz[n+m]*self.vz[n]     ## Dot product of the velocities at distance m
            C = 1./self.nofp * C.sum() ##Removed indentetation this should be checked.
            C/=(M-m)
            vaf.append(C)


        pad = np.zeros(max*5)
        padTime = (np.array(time)+max)

        
        vaf = np.append(vaf, pad)
        #time = np.append(time, padTime)

        time = np.arange(len(vaf))
        
        #vaf = self.smooth(vaf)
        
        fft =ftb.fft3DBox(vaf)
        #print (freqs)

        F = np.fft.fftn(vaf)
        F = np.fft.fftshift( F )
        F= (np.abs(F))**2
        F=np.log(F)
        
        vv = np.argsort(np.array(vaf))
        tt= np.array(time)[vv]
        
        freqs = (np.fft.fftfreq(fft.size, 1.))**2
        #plt.loglog(freqs, fft)
        #plt.plot(freqs, F)
        #plt.show()
        #self.plot(vaf, time, fft, freqs, M)

        return time, vaf
    
    def calculate(self):
        vaf = []
        time = []
        bin = 1.
        max = 500
        velvel = 0
        i=0
        points =0
        for dumpFile in sorted(glob.glob(os.path.join(self.folderID+'/vel/', '*.dump')), key=os.path.getmtime):
            i+=1.
            r = i/bin
            if (int(r)-r) == 0 and points<max:
                points+=1
                print dumpFile
                a,m, t, vx, vy, vz, lbox = rd.read(dumpFile)
                #index = np.argsort(a)   ##sorted index according to the ID name fro each particle
                #vx = vx[index]
                #vy = vy[index]
                #vz = vz[index]
                C = self.vx*vx + self.vy*vy + self.vz*vz
                C = 1./self.nofp * C.sum()
                #time.append(int(dumpFile.split(".")[1]))   # position 1 should be the time step.
                time.append(i)
                vaf.append(C)
                if velvel ==1:  ###This is not working so far, since the ordering of the dump files is affecting the sequence of computations, therefore the previous time step may not be the correct one for a given dump file.
                    self.vx = vx
                    self.vy = vy
                    self.vz = vz
        
        pad = np.zeros(max*5)
        padTime = (np.array(time)+max)

        
        vaf = np.append(vaf, pad)
        #time = np.append(time, padTime)

        time = np.arange(len(vaf))
        
        #vaf = self.smooth(vaf)
        
        fft =ftb.fft3DBox(vaf)
        #print (freqs)


        F = np.fft.fftn(vaf)
        F = np.fft.fftshift( F )
        F= (np.abs(F))**2
        F=np.log(F)
        
        vv = np.argsort(np.array(vaf))
        tt= np.array(time)[vv]


        
        freqs = (np.fft.fftfreq(fft.size, 1.))**2
        #plt.loglog(freqs, fft)
        #plt.plot(freqs, F)
        #plt.show()
        
        self.plot(vaf, time, -fft, freqs, max)
        #return vaf, time


    def smooth(self, vaf):
        M = len(vaf)
        sqtp = np.sqrt(2*np.pi)
        self.sigma*=M
        for i in range(M):
            vaf[i]*=np.exp(-i**2/(2*self.sigma**2))/(self.sigma*sqtp)
            #for m in range (1,M):
            #vaf/=vaf[0]

        return vaf
    
    def plot(self, vaf, t, fft, freqs, max):
        self.id = "vacf"+self.id
        fftID = "FFTvacf"+self.id

        pa.plotTrajectories(vaf, t, 't',"C(t)", self.id, max, True, 'o')
        pa.plotTrajectories(fft, freqs, 'freq',"fft(C)", fftID,0, True, 'o')


    def test(self, fID, fisrt):

        folderID   = fID #sys.argv[1]
        first = first #sys.argv[2]
        v = vAcF(folderID, first)
        v.calculateZ()
