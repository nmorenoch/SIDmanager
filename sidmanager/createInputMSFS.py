import os, re, sys
import numpy as np

kbt = 0.1
kbtbuf = 0.1
m=0.04
vther = (kbt/m)**0.5


scriptSource = "lmpScripts/"
varToChange = np.round(vther*np.arange(1,50,10),3)
#varToChange = np.round(vther*np.linspace(0.25,10,10),3)
varToChange = np.round(vther*np.linspace(0.25,20,20),3)


lcore = 8.0 #varToChange[0]/gamma
lbuf  = 2.2
lcon  = 0.8

gamma = np.linspace(0.005,2.0,10)
varToChange = np.round(gamma*(lcore/2.),3)



Lcore = np.arange(3,27,3)  #for constant vx = 1.7vther
lbuf = 1.2
lcon = 0.8
#gamma = 0.05
varToChange = np.round(1.7*vther*np.ones(len(Lcore)),3)



lastsid =int(sys.argv[1])
pairT = sys.argv[2]
sourcBC = sys.argv[3]
sourcFine = sys.argv[4]
cores = int(sys.argv[5])

for sid in range(1+lastsid,len(varToChange)+1+lastsid):
    ### Creating static files and folders
    
    line = "cp lmpScripts/regionsDef.lmp %s/" %sid
    os.system(line) 
    line = "cp lmpScripts/atomsData %s/" %sid
    os.system(line) 
    line = "cp lmpScripts/mdpd %s/" %sid
    os.system(line)
    line = "cp lmpScripts/meso %s/" %sid
    os.system(line)
    line = "cp lmpScripts/computeFix %s/" %sid
    os.system(line)

    if sourcBC=='bcDefinition':
        f = open("lmpScripts/%s"%sourcBC, 'r')
        output = open("%s/%s"%(sid,sourcBC), 'w')
        for line in f:
            l = line
            if re.match('variable vx1', l) != None:
                output.write("variable vx1 equal %s \n"%varToChange[sid-1-lastsid])
            elif re.match('variable vx2', l) != None:
                output.write("variable vx2 equal %s \n"%varToChange[sid-1-lastsid])
            elif re.match('variable vx3', l) != None:
                output.write("variable vx3 equal -%s \n"%varToChange[sid-1-lastsid])
            elif re.match('variable vx4', l) != None:
                output.write("variable vx4 equal -%s \n"%varToChange[sid-1-lastsid])    
            else:    
                output.write(l.rstrip() + '\n')
        f.close()
        output.close()

    elif sourcBC=='extensionBC':
        f = open("lmpScripts/%s"%sourcBC, 'r')
        output = open("%s/%s"%(sid,sourcBC), 'w')
        for line in f:
            l = line
            if re.match('variable vx1', l) != None:
                output.write("variable vx1 equal -%s \n"%varToChange[sid-1-lastsid])
            elif re.match('variable vx2', l) != None:
                output.write("variable vx2 equal %s \n"%varToChange[sid-1-lastsid])
            elif re.match('variable vx3', l) != None:
                output.write("variable vx3 equal %s \n"%varToChange[sid-1-lastsid])
            elif re.match('variable vx4', l) != None:
                output.write("variable vx4 equal -%s \n"%varToChange[sid-1-lastsid])
            elif re.match('variable vy1', l) != None:
                output.write("variable vy1 equal -%s \n"%varToChange[sid-1-lastsid])    
            elif re.match('variable vy2', l) != None:
                output.write("variable vy2 equal -%s \n"%varToChange[sid-1-lastsid])
            elif re.match('variable vy3', l) != None:
                output.write("variable vy3 equal %s \n"%varToChange[sid-1-lastsid])
            elif re.match('variable vy4', l) != None:
                output.write("variable vy4 equal %s \n"%varToChange[sid-1-lastsid])
            else:    
                output.write(l.rstrip() + '\n')
        f.close()
        output.close()

    
    f = open("lmpScripts/%s"%sourcFine, 'r')
    output = open("%s/%s"%(sid,sourcFine), 'w')
    for line in f:
                l = line
                if re.match('variable pairStyle', l) != None:
                    output.write("variable pairStyle string '%s' \n"%pairT)
                elif re.match('variable T ',l) != None:
                    output.write("variable T equal %s \n"%kbt)
                elif re.match('variable Tbuf',l) != None:
                    output.write("variable Tbuf equal %s \n"%kbtbuf)
                elif re.match('variable coreL',l) != None:
                    output.write("variable coreL equal %s \n"%Lcore[sid-1-lastsid])
                elif re.match('variable consL',l) != None:
                    output.write("variable consL equal %s \n"%lcon)
                elif  re.match('variable buffL',l) != None:
                    output.write("variable buffL equal %s \n"%lbuf)
                elif  re.match('include bcDefinition',l) != None:
                    output.write("include %s \n"%sourcBC)    
                else:    
                    output.write(l.rstrip() + '\n')
    f.close()
    output.close()
    
    print(sid)

