#!/usr/bin/python

import data_runs

def getRuns(tag):
    Runs = {}
    if tag == 'low':
        print('OK, will run over low momenta runs!')
        Runs = data_runs.LowRuns
    elif tag == 'high':
        print('OK, will run over higher momenta runs!')
        Runs = data_runs.HighRuns
    elif tag == 'p':
        print('OK, will run over proton runs!')
        Runs = data_runs.ProtonRuns
    else:
        print('Wrongly specified momentum range! Must be low/high/p')
    return Runs
