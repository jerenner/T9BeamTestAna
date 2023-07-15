#!/usr/bin/python

# jiri kvita 2023

# from collections import OrderedDict

####################################################################

# numbers from
# ./python/parseWcteDaqRunsHtml.py
# using share/wcte-daq.html
# which is htmp page saved from the wcte daq online page

momentaDict = {300: [283], 340: [282], 380: [281], 410: [280], 460: [279], 500: [278], 540: [277], 600: [276], -660: [274], -830: [273], -600: [272], -690: [271], -740: [270, 269], -800: [268], -860: [267], -890: [266], -940: [265], 1000: [264], -1000: [263], 660: [262], 830: [261], 940: [260], 890: [259], 690: [258], 740: [257], 860: [256], 800: [255]}
runsDict = {283: 300, 282: 340, 281: 380, 280: 410, 279: 460, 278: 500, 277: 540, 276: 600, 274: -660, 273: -830, 272: -600, 271: -690, 270: -740, 269: -740, 268: -800, 267: -860, 266: -890, 265: -940, 264: 1000, 263: -1000, 262: 660, 261: 830, 260: 940, 259: 890, 258: 690, 257: 740, 256: 860, 255: 800}


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
