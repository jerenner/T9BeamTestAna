#!/usr/bin/python
# jk

from data_runs import *

import os, sys

os.system('mkdir -p output')

dryrun = True # default
#dryrun = False # Careful!!

if not dryrun:
    print("WARNING, this will run the waveform analysis and MakeDataPlots for all runs in your lists/, do you really wish to continue? Y/n")
    a = input()
    if not (a == 'Y' or a == 'y'):
        exit(1)
    

print("===> Waveforms analysis:")

for xlistname in  os.popen('cd lists/ ; ls list*.txt'):
    
    print('######################################################')
    listname = xlistname[:-1]
    base = listname.replace('.txt','')
    #print('base: ', base)
    srun = base.replace('list_root_run_000', '')
    momentum = getMomentum(srun)
    
    cmd = "./bin/waveform_analysis.app -i lists/{} -o output/output_{}.root -c config/config.json".format(listname,base)
    print(cmd)
    if not dryrun:
        os.system(cmd)

    cmd='root -l -b -q "scripts/MakeDataPlots.C(\\"output/output_{}.root\\", {})"'.format(base,momentum)
    print(cmd)
    if not dryrun:
        os.system(cmd)
        
    cmd='./python/quickPlots1d.py histos/output_{}_plots.root'.format(base)
    print(cmd)
    #if not dryrun:
    #    os.system(cmd)
        





