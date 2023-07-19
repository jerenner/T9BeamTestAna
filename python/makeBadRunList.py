#!/usr/bin/python

import os, sys

# we assume files with file zie in kilo bytes (K) are too small and problematic:

badntuples = os.popen("cd output ; ls -Slh *.root | grep K | awk '{print($9);}' ")


badruns = []
for xbadntuple in badntuples.readlines():
    badntuple = xbadntuple[:-1]
    srun = badntuple.split('_')[1].replace('.root','')
    strun = ""
    run = -1
    #print(srun)
    try:
        run = int(srun)
    except:
        print('Could not convert the run number for {}'.format(badntuple))
    if run > 0:
        badruns.append(run)

print('#!/usr/bin/python')
print('badruns = ', badruns)
