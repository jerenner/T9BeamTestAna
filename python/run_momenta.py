#!/usr/bin/python

import data_runs

import os, sys

RunAll = False
#RunAll = True

argv = os.sys.argv
print(argv)
if len(argv) > 1:
         if argv[1] == '1' or argv[1] == '1' or argv[1] == 'y' or argv[1] == 'Y':
                  RunAll = True

# make a link to you data directory
datadir = 'data/'
listdir = 'lists/'

#################################################

for p in data_runs.Runs:
         print('*** {} ***'.format(p))
         runspills = data_runs.Runs[p]
         runs = runspills[0]
         ip = int(p.replace('n','').replace('p',''))
         if 'n' in p:
                  ip = -ip
         cmds = []
         
         cmds.append('ls {}*clean*.root | egrep "{}" > {}list_{}.txt'.format(datadir,runs,listdir,p))
         cmds.append('./bin/waveform_analysis.app -i {}list_{}.txt -o output_{}.root -c config/config.json'.format(listdir,p,p))
         cmds.append('root -l -q "scripts/MakeDataPlots.C(\\"output_{}.root\\", {})"'.format(p,ip))
         cmds.append('root -l -b -q "scripts/FitTOF.C(\\"output_{}_plots.root\\",{})"'.format(p, ip))
         
         for cmd in cmds:
                  print(cmd)
                  if RunAll:
                           os.system(cmd)


print('DONE!')
print('Consider:')
print('python ./python/plotFromAscii.py')
