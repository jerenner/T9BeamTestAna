#!/usr/bin/python

import os, sys

RunAll = False
#RunAll = True

# make a link to you data directory
datadir = 'data/'
listdir = 'lists/'

Runs = {

         # jk
         '320p' : '142|166|189',
	 '340p' : '168',
	 '300n' : '148|156',
	 '320n' : '149|158|188',
       }
#################################################

for p in Runs:
         print('*** {} ***'.format(p))
         runs = Runs[p]
         cmds = []
         cmds.append('ls {}*clean*.root | egrep "{}" > {}list_{}.txt'.format(datadir,runs,listdir,p))
         cmds.append('./bin/waveform_analysis.app -i {}list_{}.txt -o output_{}.root -c config/config.json'.format(listdir,p,p))
         ip = int(p.replace('n','').replace('p',''))
         if 'n' in p:
                  ip = -ip
         cmds.append('root -l -q "scripts/MakeDataPlots.C(\\"output_{}.root\\", {})"'.format(p,ip))
         cmds.append('root -l -b -q "scripts/FitTOF.C(\\"output_300n_plots.root\\",{})"'.format(ip))
         for cmd in cmds:
                  print(cmd)
                  if RunAll:
                           os.system(cmd)



