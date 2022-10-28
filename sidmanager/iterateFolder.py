import os
import glob
 
path = '/Users/morenon/Desktop'
print path
for infile in glob.glob( os.path.join(path, '*.png') ):
    f = (infile.split('/')[4]).split('.')[0]
    print f
