import sys, os

for i in range(1,len(sys.arvg)):
    os.system('cd '+int(sys.argv[i])+'\n'+'qsub script.lone')

print 'done'
