#!/usr/bin/python
# jk

from data_runs import *

import os, sys

os.system('mkdir -p output')

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

    cmd='root -l -b -q "scripts/MakeDataPlots.C(\\"output/output_{}.root\\", {})"'.format(base,momentum)
    print(cmd)
    cmd='./python/quickPlots1d.py histos/output_{}_plots.root'.format(base)
    print(cmd)





