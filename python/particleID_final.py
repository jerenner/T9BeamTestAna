import sys
import os

import json5

import numpy as np

import lowMomentum_analysis as lm

config_file_name = sys.argv[1]

with open(config_file_name) as file:
    config = json5.load(file)["LowMomentumAnalysis"]
print(config)

particleID_analysis = lm.LowMomentumAnalysis(config)
particleID_analysis.openDataFile()

print(particleID_analysis.getBranchList(0))
if config["runNumber"]<579:
    sumDownsteamACTs = 0
    #adding all the data from the downstream ACTs
    for i in [4, 5, 6, 7]:
        sumDownsteamACTs+=particleID_analysis.getDataFrameDetector(i)["matchedHit0_WindowIntPE"]

    particleID_analysis.addBranchToAllDetectors("ACT2L+ACT2R+ACT3L+ACT3R_WindowIntPE", sumDownsteamACTs)

if particleID_analysis.BoundsAreAvailable:
    print("yes")
    particleID_analysis.addOperationBranchToAllDetectors("PMThit0-WindowCentralTime", "peakHit0_SignalTimeCorrected", "-", "matchedHit0_WindowCentralTimeCorrected")


print(particleID_analysis.getBranchList(0))

#plotting below:
subtitle = "ACTs: -16 to 45 ns - TOFs -5 to 15 ns - LG -16 to 55 - HD -26 to 45"
extraLabel = "v5"

if config["runNumber"]!=502 and config["runNumber"]!=676:
    #dark rate runs -> not wanted
    particleID_analysis.plotSomeDetectorHist1DfromData(["matchedHit0_TOF", "matchedHit1_TOF"], [0], "TOF_%s_%s"%(config["channelNames"][0], extraLabel), False, 10, None, subtitle, "log", 4, 400)

    if config["runNumber"]<579:
        #only dpo these plots for LM runs
        particleID_analysis.plot2DHistFromBranches(2, "matchedHit0_WindowIntPE", 4, "ACT2L+ACT2R+ACT3L+ACT3R_WindowIntPE", "(PE)", "(PE)", subtitle, True, [100, 100], [[0, 40], [0,150]])
        
        #Make 2D scatter branches from the branches that we have, with different detectors The z is a third variable 
        particleID_analysis.plot2DScatterFromBranches(0, "matchedHit0_TOF", -1, "peakHit0_IntPE", 4, "ACT2L+ACT2R+ACT3L+ACT3R_WindowIntPE", "(ns)", "(PE)", "(PE)", subtitle, True)

particleID_analysis.nCoincidenceSelection()

upperBounds = [10, 10, 10, 10, 150, 150, 150, 150, 80, 80, 80, 80, 80, 80, 80, 80, 10, 10, 20]
for detectorID in np.arange(0, len(config["channelNames"]), 1):
    if config["channelNames"][detectorID] != None:
        particleID_analysis.plotSomeDetectorHist1DfromData(["peakHit0_IntPE", "matchedHit0_WindowIntPE", "WholeWaveformIntPE"], [detectorID], "%s_charge_all_%s"%(config["channelNames"][detectorID], extraLabel), False, -1, 200, subtitle, "log", 4, 200)

        particleID_analysis.plotSomeDetectorHist1DfromData(["PMThit0-WindowCentralTime"], [detectorID], "%s_time_all_%s"%(config["channelNames"][detectorID], extraLabel), False, -20, 20, subtitle, "", 4, 200, True)

#plotting below:
subtitle = "NON-PROTONS ONLY ACTs: -16 to 45 ns - TOFs -5 to 15 ns - LG -16 to 55 - HD -26 to 45"

particleID_analysis.selectBasedOnCondition(0, "matchedHit0_TOF", '<', 15)

for detectorID in np.arange(0, len(config["channelNames"]), 1):
    if config["channelNames"][detectorID] != None:
        particleID_analysis.plotSomeDetectorHist1DfromData(["matchedHit0_WindowIntPE", "WholeWaveformIntPE"], [detectorID], "%s_charge_nonProtons_%s"%(config["channelNames"][detectorID], extraLabel), False, -1, 200, subtitle, "log", 4, 200)

        particleID_analysis.plotSomeDetectorHist1DfromData(["PMThit0-WindowCentralTime"], [detectorID], "%s_time_nonProtons_%s"%(config["channelNames"][detectorID], extraLabel), False, None, None, subtitle, "", 4, 200, True)