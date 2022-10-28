import sys
import os
import numpy as np

"""
This script create folder with the name as an integer from init to fin
"""


def createFolder(init, fin,variables,folderScripts):
    
    for i in np.arange(init,fin):
        try:
            #Main folder 
            os.mkdir(i.__str__())
            #Internal folders: dumps, force, vel and restarts
            os.mkdir(i.__str__()+"/dumps")
            os.mkdir(i.__str__()+"/restart")
            #os.mkdir(i.__str__()+"/boundary")
            #os.mkdir(i.__str__()+"/restart")
            #os.mkdir(i.__str__()+"/vis")

        except Exception as e:
            print (e)
            continue
    counter = 0
    for i in np.arange(init,fin):

            ####This is assuming the following definition are static. Changes on this lines can be done to include specific creations.
            sid = i
            line = "cp %s%s %s/" %(folderScripts, variables['regionsFile'].T[counter],sid)
            os.system(line) 
            line = "cp %s%s %s/" %(folderScripts, variables['atomsFile'].T[counter], sid)
            os.system(line) 
            line = "cp %s%s %s/" %(folderScripts, variables['pairStyle'].T[counter], sid)
            os.system(line)
            line = "cp %s%s %s/" %(folderScripts, variables['computesFile'].T[counter],sid)
            os.system(line)

            counter+=1

if sys.argv[0] == "/Users/morenon/mypythoncodes/createFolders.py":
    print (sys.argv[0])
    init = int(sys.argv[1])
    fin = int(sys.argv[2])
    #createFolder(init, fin)
else: print ("create folder: It suppose to be running from another script?" )
