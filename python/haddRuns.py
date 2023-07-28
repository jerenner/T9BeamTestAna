#!/usr/bin/python
#jk 25.7.2023

import os, sys

from data_runs_dicts import *


dryRun = True
#dryRun = False

for p in pnrunsDict:
    for n in pnrunsDict[p]:
        sp = str(abs(p))
        if p < 0:
            sp = sp + 'n'
        else:
            sp = sp + 'p'
        tag = '{}_n{}'.format(sp,n).replace('.','p')
        cmd = 'hadd -f output/merged_ntuple_{}.root '.format(tag)
        for run in pnrunsDict[p][n]:
            cmd = cmd + f'output/ntuple_000{run}.root '
        print(cmd)
        if not dryRun:
            os.system(cmd)
