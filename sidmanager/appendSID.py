"""
03-29-2013
@author: nmorenoc
         nicolas.morenochaparro@kasut.edu.sa

It append the last SID to a plain text.
If flag in header is setted equal to 1, then the first line with the data being stored will be created again

"""



import sys, os
import numpy as np
import datetime
from tempfile import mkstemp
from shutil import move
import pandas as pd


def copyFromSID(SID, file):
            #data = []
            #f = open(file, 'r')
            #for i in f:
            #    data.append(i.split('\t'))
            #f.close()
            #data  = np.array(data)
            
            data = pd.read_table(file,sep='\t', header=0, index_col=0)

            line = data[data[:,0]==SID]
            if len(line)!=1:   ##If there is not SID to take as reference, then take the last one, and increase the SID, I must remove the new line command by the spliting of the last parameter. Durign the storage I include an space before the new line. 
                line = data[-1]
                line[0]  = (int(line[0])+1).__str__()
                line[-1] =  line[-1].split(" ")[0] 
            if len(line)==1:
                line[0][0]  = (int(data[-1][0])+1).__str__()   ##here take the info SID of the last line stored in the details file, and increase the counter.
                #print line[0][-1]    ###Due to the way the array is stored I should specify 2 indices, however is only 1D
                line[0][-1] = line[0][-1].split(" ")[0] 
            return line 

def header(logfile, flag, head):
    if hasattr(head, "__len__"):
        h=head
    else:
        h =  ['SID', 'date', 'num_dim','cut_off', 'max_phys','num_part','coll_shape', 'coll_radius', 'coll_x', 'wall_noslip','shear_type','bcdef']


    all = np.array(h)

    if flag =='1':
        new_file = open(logfile, 'w')
        all.tofile(new_file, sep='\t', format='%s')
        new_file.write('\n')
        new_file.close()
        print ("First time created logfile")

    """
    # If flag 2, this mean an exitent file will be modified in the header, this is only if new variables to track have been added. This flag has not been externally
    associated yet to any script
    #####Code adapted from : http://stackoverflow.com/questions/39086/search-and-replace-a-line-in-a-file-in-python 05-05-2013
    #Create temp file
    """
    if flag =='2':  
        #fh, abs_path = mkstemp()
        #new_file = open(abs_path,'w')
        #old_file = open(logfile)
        #all.tofile(new_file, sep='\t', format='%s')
        #new_file.write('\n')
        #for line in old_file:
        #    new_file.write(line)
            #close temp file
        #new_file.close()
        #os.close(fh)
        #old_file.close()
        #Remove original file
        #os.remove(logfile)
        #Move new file
        #move(abs_path, logfile)

        fh, abs_path = mkstemp()
        new_file = open(abs_path,'w')
        data = pd.read_table(logfile,sep='\t', header=0) #Read the table and formatted as a DataFrame. In this case I dont indicate the first column as index so I can handle with all the columns as variables
        DF=pd.DataFrame(data)
        newData = pd.DataFrame(columns=all)#creating a new data frame with the new set of modified variables, and combine them.
        newData = newData.append(DF)

        newData.columns.tofile(new_file, sep='\t')
        print >>new_file
        for line in newData.values:
            line.tofile(new_file,sep='\t')
            print >>new_file
        new_file.close()
        os.remove(logfile)
        #Move new file
        move(abs_path, logfile)
    ####################################

def appendX(logFile, input):
        o = open(logFile,'a')
        all = np.array(input)
        all.tofile(o, sep='\t')
        o.write(' \n')
        o.close()

def append(logFile, input):
        o = open(logFile,'a')
        all = np.array(input)
        line=""
        nvars = len(all)
        for i in range(nvars):
            if isinstance(all[i], int):
               line+="%d" % all[i]
            elif isinstance(all[i],float):
               line+="%g" % all[i]
            else:
               line+="%s" % all[i]
            if i<nvars-1:
               line+="\t"
        o.write(line)
        o.write('\n')
        o.close()


def createHTML(file):
    data = pd.read_csv(file,sep='\t', header=0, index_col=0)
    DF=pd.DataFrame(data)
    text = DF.to_html(justify='left', col_space=80)
    f=open(file.replace('.sid','.html'), 'w')
    f.write(text)
    f.close()


def doSID(i, fSIDdet, flagHeader='0', copyFrom='0', head = None):

    if os.path.isfile(fSIDdet)==False:
        print ("Creating a file for SID details: "+fSIDdet)
        flagHeader='1'        
    header(fSIDdet, flagHeader, head)   ##if 1 then create a header again deleting the existent information in the file WARNING!! a safer implementation is required.

    if copyFrom =='0':
        line = i
    else:
        d = copyFromSID(copyFrom, fSIDdet)  ##extracting a list with all the info of a SID. Here SID is passed as a string, and therefore needs to be compared in this way.
        line = d
        #print line
    append(fSIDdet, line)
    createHTML(fSIDdet)



if sys.argv[0]=="appendSID.py":
    now = datetime.datetime.now()
    date = "%d-%d-%d" % (now.month, now.day, now.year)
    i = ['ALLSID.sid', 2, date, 50, 317, 16, 0.15, 1000, 0.25, 50, 0.75, 5, 10, 121331, 27, 35, 35, 25, 0, 0, 0, 0, 0, 0, 0, 3, 180, 1, 'TACC', '--', 0, 0, 0.15]
    flag = sys.argv[1]
    
    sid = SID(i)
    sid.header(sys.argv[1])   ##if 1 then include the header
    
    d = sid.copyFromSID('5', 'ALLSID.sid', sys.argv[1])  ##extracting a list with all the info of a SID. Here SID is passed as a string, and therefore needs to be compared in this way.
    
    if sys.argv[1]=='0':
        #    print d.tolist()
        line = np.append('ALLSID.sid', d)
        print (line)
        sid = SID(line)
    sid.append()
else: print ("appendSID: It suppose to be running from another script?" )
