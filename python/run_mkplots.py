#!/usr/bin/python
# jk

from data_runs import *
from data_badruns import *

import os, sys

os.system('mkdir -p output')

dryrun = True # default
#dryrun = False # Careful!!

if not dryrun:
    print("WARNING, this will MakeDataPlots for all runs in your output/, do you really wish to continue? Y/n")
    a = input()
    if not (a == 'Y' or a == 'y'):
        exit(1)
    
for xlistname in  os.popen('cd output/ ; ls ntuple_*.root'):
    
    print('######################################################')

    rfilename = xlistname[:-1]
    base = rfilename.replace('.root','')
    #print('base: ', base)
    srun = ''
    tokens = base.split('_')
    for token in tokens:
        if '00' in token:
            srun = token.replace('000','')

    run = -1
    try:
        run = int(srun)
    except:
        print('# ERROR getting run number from the file name!')
    if run > 0:
        if run in badruns:
            print('# SKIPPING bad run {}'.format(run))
            continue

    momentum = getMomentum(srun)

    cmd='root -l -b -q "scripts/MakeDataPlots.C(\\"output/{}\\", {})"'.format(rfilename,momentum)
    print(cmd)
    #if not dryrun:
    #    os.system(cmd)

    hfilename = rfilename.replace('.root','_plots.root')
    cmd='./python/quickPlots1d.py histos/{}'.format(hfilename)
    print(cmd)
    if not dryrun:
        cmd = cmd + ' -b'
        os.system(cmd)
        





