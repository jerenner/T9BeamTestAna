#!/usr/bin/python
# jk

from data_runs import *

from data_badruns import *

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

#for xlistname in  os.popen('cd lists/ ; ls list*.txt'):
for xlistname in  os.popen('cd data/ ; ls root_run_*.root'):
    
    print('######################################################')

    # c++ version:
    #listname = xlistname[:-1]
    #base = listname.replace('.txt','')
    ##print('base: ', base)
    #srun = base.replace('list_root_run_000', '')

    rfilename = xlistname[:-1]
    base = rfilename.replace('.root','')
    #print('base: ', base)
    srun = base.replace('root_run_000', '')

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

    # C++ verision depricated:
    #cmd = "./bin/waveform_analysis.app -i lists/{} -o output/output_{}.root -c config/config.json".format(rfilename,base)

    # python version by Nick:
    cmd = 'python ./python/new_analysis/process_waveform_analysis.py data/{} config/config.json output/output_{}.root'.format(rfilename, base)
    
    print(cmd)
    if not dryrun:
        os.system(cmd)

    cmd='root -l -b -q "macros/MakeDataPlots.C(\\"output/output_{}.root\\", {})"'.format(base,momentum)
    print(cmd)
    if not dryrun:
        os.system(cmd)
        
    cmd='./python/quickPlots1d.py histos/output_{}_plots.root'.format(base)
    print(cmd)
    if not dryrun:
        cmd = cmd + ' -b'
        os.system(cmd)
        





