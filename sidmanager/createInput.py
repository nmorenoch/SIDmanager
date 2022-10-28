'''
Created on Feb 7, 2012

@author: morenon

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>
'''
import re 
import sys, os
import shutil
import numpy as np

i=0

def setIndex(first, last, source, variables, names):
    path = os.getcwd()+"/"
    var = len(names)
    counter = 0
   
    for i in range(first,last):
        ##shutil.copyfile('/Users/morenon/Documents/workspaceJava/SimulationSettings/44', '/Users/morenon/Documents/workspaceJava/SimulationSettings/'+i.__str__())     
        #        f = open('/Users/morenon/Documents/workspaceJava/SimulationSettings/'+(source).__str__()+'/micelle.in', 'r')
        f = open(path+source, 'r')
        output = open(path+i.__str__()+'/colloid.in', 'w')
        for line in f:
            l = line
            if re.match('read', l) != None:
                #print variables[var-1][counter]
                output.write("read_data    "+variables['dataName'].T[counter]+'\n')  #Here It is being assummed the las value in the variables array is the name of the data file 
            elif re.match('#VARIABLES', l) != None:
                for V in range(var):
                    try: float(variables.T[counter][V])
                    except: new = variableStringCode(variables.T.index[V], variables.T[counter][V])     
                    else: new = variableCode(variables.T.index[V], variables.T[counter][V])
                    print >> output, new
            else:    
                output.write(l.rstrip() + '\n')
        f.close()
        output.close()
        counter+=1
    print "fin"

def variableCode(varName, val):
    line = "variable %s equal (%s)" % (varName, val)
    return line

def variableStringCode(varName, val):
    line = "variable %s string '%s'" % (varName, val)
    return line

def evaporationCode(rate, initFrac, finFrac, relaxTime):
    return 0

def addWallCode():
    return 0

def addRampAijCode():
    return 0 
    

def setIndexGN(first, last, source, aii, aht, ahs, ats):
    dumpt   = '1000'
    dumpv   = '500000'
    dumpf   = '500000'
    dumpEnd = '500000'
    dumpG   = '1000'

    path = os.getcwd()+"/"
    q = 0   ##counter for variables to be change in each input
    for i in range(first,last):
        
        ##shutil.copyfile('/Users/morenon/Documents/workspaceJava/SimulationSettings/44', '/Users/morenon/Documents/workspaceJava/SimulationSettings/'+i.__str__())     
        #        f = open('/Users/morenon/Documents/workspaceJava/SimulationSettings/'+(source).__str__()+'/micelle.in', 'r')
        f = open(path+(source).__str__()+'/grafted.in', 'r')
        output = open(path+i.__str__()+'/grafted.in', 'w')
        for line in f:
            l = line
            if re.match('read', l) != None:
                output.write("read_data    data_"+i.__str__()+".gn"+'\n')     
            elif re.match('dump 1', l) != None:
                output.write("dump 1 micelle custom "+dumpt+" ./dumps/gn"+i.__str__()+"_1.*.dump id mol type x y z" + '\n')
            elif re.match('dump 2', l) != None:
                output.write("dump 2 micelle custom "+dumpv+" ./vel/mVel"+i.__str__()+"_1.*.dump id mol type vx vy vz" + '\n')
            elif re.match('dump 3', l) != None:
                output.write("dump 3 micelle custom "+dumpf+" ./force/mForce"+i.__str__()+"_1.*.dump id mol type fx fy fz" + '\n')
            elif re.match('dump 4', l) != None:
                output.write("dump 4 all custom "+dumpEnd+" gnEnd"+i.__str__()+"_1.dump id mol type x y z" + '\n')
            elif re.match('compute $i', l) != None:
                output.write('compute $i micelle gyration/molecule'+'\n')     
            elif re.match('fix rofGy1', l) != None:
                output.write("fix rofGy1 micelle ave/time 100 10 "+dumpG+" c_1 file rog1.rog mode vector"+'\n') 
            elif re.match('variable aii', l) != None:
                text = ('variable aii equal (%f)'+'\n')  % aii[q]
                output.write(text)     
            elif re.match('variable aht', l) != None:
                text = ('variable aht equal (%f)'+'\n')  % aht[q]
                output.write(text)
            elif re.match('variable ahs', l) != None:
                text = ('variable ahs equal (%f)'+'\n')  % ahs[q]
                output.write(text)
            elif re.match('variable ats', l) != None:
                text = ('variable ats equal (%f)'+'\n')  % ats[q]
                output.write(text)
            else:    
                output.write(l.rstrip() + '\n')

        q+=1
        f.close()
    print "fin"




if sys.argv[0] == "/Users/morenon/mypythoncodes/createInput.py":
    first = int(sys.argv[1])
    last = int(sys.argv[2])
    source = int(sys.argv[3])
    type = sys.argv[4]
    if type=='gn':
       aii = np.ones(last-first)*25.0
       aht = aii
       ats = aii
       ahs = [15.0, 20.0, 25.0, 26.0, 27.0, 27.5, 28.0, 29.0, 30.0, 35.0]
       setIndexGN(first, last, source, aii, aht, ahs, ats)
    if type=='micelle':
        setIndex(first, last, source)
else: print "createInput: It suppose to be running from another script?" 


        
