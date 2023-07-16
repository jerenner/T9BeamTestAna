#!/usr/bin/python

# jiri kvita 2023

# from collections import OrderedDict

from data_runs_dicts import *

####################################################################

def getMomentum(srun):
    momentum = 1000
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
