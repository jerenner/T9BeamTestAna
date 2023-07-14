#!/usr/bin/python
# jk

import os, sys

os.system('mkdir -p output')

print("===> Waveforms analysis:")
for xlistname in  os.popen('cd lists/ ; ls list*.txt'):
    listname = xlistname[:-1]
    base = listname.replace('.txt','')
    cmd = "./bin/waveform_analysis.app -i lists/{} -o output/output_{}.root -c config/config.json".format(listname,base)
    print(cmd)

# TBC
momentum="1000"

print("===> ToF analysis, histogramming:")
for xlistname  in os.popen('cd lists/ ; ls list*.txt'):
    listname = xlistname[:-1]
    base=listname.replace('.txt', '')
    cmd='root -l "scripts/MakeDataPlots.C(\\"output/output_{}.root\\", {})"'.format(base,momentum)
    print(cmd)


