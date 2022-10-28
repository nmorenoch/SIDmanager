'''
Created on Jun 8, 2020

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

def setIndex(first, last, source, variables,names,folderScripts):
    path = os.getcwd()+"/"
    var = len(names)
    counter = 0
    for i in range(first,last):

        f = open(path+folderScripts+source, 'r')
        output = open(path+i.__str__()+'/%s'%source, 'w')
        for line in f:
            l = line
         #   if re.match('read', l) != None:
         #       #print variables[var-1][counter]
         #       output.write("vari    "+variables['dataName'].T[counter]+'\n') 
            if re.match('#VARIABLES', l) != None:
                for V in range(var):
                    try: float(variables.T[counter][V])
                    except: new = variableStringCode(variables.T.index[V], variables.T[counter][V])     
                    else: new = variableCode(variables.T.index[V], variables.T[counter][V])
                    output.write(new+'\n')
            else:    
                output.write(l.rstrip() + '\n')
        f.close()
        output.close()
        counter+=1
    print ("fin")

def variableCode(varName, val):
    line = "variable %s equal %s" % (varName, val)
    return line

def variableStringCode(varName, val):
    line = "variable %s string '%s'" % (varName, val)
    return line

    
if sys.argv[0] == "/Users/morenon/mypythoncodes/createBC.py":   #####This is not doing anything need to be updated to be consistent with this script
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
else: print ("createInput: It suppose to be running from another script?") 


        
