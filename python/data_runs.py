#!/usr/bin/python

# jiri kvita 2023

# from collections import OrderedDict

from data_runs_dicts import *


####################################################################

def getTarget(run):
    if run >= 537:
        return 'Target3'
    else:
        return 'Target1'

####################################################################

def getMergedMomentum(srun):
    # we assume srun is actually an expression like 240p or 1000n:
    momentum = 0
    if srun[-1:] == 'p' or srun[-1:] == 'n':
        try:
            momentum = int(srun[:-1])
            if srun[-1:] == 'n':
                momentum = momentum*(-1)
            return momentum
        except:
            print('still an issue to deduce momenyum from ', srun)
    return momentum
####################################################################

def getMomentum(srun):
    momentum = 0
    run = -999
    try:
        run = int(srun)
    except:
        print(f'ERROR converting run {srun} to int!')
    try:
        momentum = runsDict[run]
    except:
        print(f'ERROR getting momentum information for run {run} from python/data_runs.py !')
    return momentum

####################################################################
def getListOfRuns(momentum):
    runs = []
    try:
        runs = momentaDict[momentum]
    except:
        print(f'ERROR getting the list of runs for momentum {momentum}')
    return runs

####################################################################
