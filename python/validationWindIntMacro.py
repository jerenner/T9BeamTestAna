import sys
import os

import json5

import numpy as np

import lowMomentum_analysis as lm

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

config_file_name = sys.argv[1]

with open(config_file_name) as file:
    config = json5.load(file)["BeamDataAnalysis"]
print(config)

leadGlassID = config["channelNames"].index("PbGlass")

#make the analysis class
BeamDataAna = lm.LowMomentumAnalysis(config)

#open the data file 
BeamDataAna.openDataFile()

#make the helper, higher level variables  
BeamDataAna.makeSumTS()
BeamDataAna.makeSumTSwindow2()
BeamDataAna.makeSumDownstreamACTs()
BeamDataAna.makeSumACT1()

#select events with exactly one coincidence (set in the config file)
BeamDataAna.nCoincidenceSelection()

#check the distribution of WholeWaveformIntPE against sumTS to check if cut is sensible
if "matchedHit0_Window2IntPE" in BeamDataAna.getBranchList(0):
    BeamDataAna.plot2DHistFromBranches((config["channelNames"].index("PbGlass")), "matchedHit0_Window2IntPE", 0, "sumTSwindow2", "(PE)", "(PE)", "sumTSwindow2_Window2IntPe_vsMatchedHit0_Window2IntPE_nCoincidence_noTSselection", True, [200, 300], [[0, 60], [0, 3000]])

#select low energy deposition in the TS but high eneough that it doesn't look like failed coincidence
BeamDataAna.TStotalChargeSelection()

#This can be useful to have a look out: timing of hits
#add a new branch with the delay between the centre of the window and the peak time. 
if BeamDataAna.WindowBoundsAreAvailable:
    BeamDataAna.addOperationBranchToAllDetectors("delayBetweenCentreWindowAndPeakTime", "peakHit0_SignalTimeCorrected", "-", "matchedHit0_WindowCentralTimeCorrected")

#based on the -16ns to 45ns window selection, ID all particles that we have in this run, values in the config file
BeamDataAna.makeAllParticleSelection()

#using all the available particles measure the beam momentum with bins of 0.1 
BeamDataAna.measureMomentumUsingTOF(0.1)

#Useful to know: config["channelNames"].index("PbGlass") gives the index of the lead glass detector.

#for all the particles, plot the 1D histogram corresponding to certain high level branches that can be useful
BeamDataAna.plotBranchHistForAllParticles(0, "sumACT1", 5, True)

if BeamDataAna.thereIsSecondWindow:
    #only if we have calculated the second integration window can look at the TS distribution 
    BeamDataAna.plotBranchHistForAllParticles(0, "sumTSwindow2", 25, True)

#Look at the windowIntPE for all the particles for some interesting detectors 
for detector in ["ACT0L", "ACT0R", "ACT1L", "ACT1R", "ACT2L", "ACT2R", "ACT3L", "ACT3R", "PbGlass"]:
    BeamDataAna.plotBranchHistForAllParticles(config["channelNames"].index(detector), "matchedHit0_WindowIntPE", 0.2, True, [0, 50])


BeamDataAna.plotBranchHistForAllParticles(0, "sumDownstreamACTs", 5, True)


#Make n= bins equally populated in terms of the Trigger scintillator 10 wholeWaveformIntPE charge 
BeamDataAna.measureElTOFresolutionFunctionOfTScharge(10)

#Output the results as a csv file
BeamDataAna.outputResults()


#If interestested, one can study the weird electrons
muonLikeWeirdElectron = BeamDataAna.muonLikeWeirdElectronArray

pionLikeWeirdElectron = BeamDataAna.pionLikeWeirdElectronArray


#raise end
fig, ax = plt.subplots(1, 1, figsize = (16, 9))
title = "TimeOfFlight"

_, muonWEmeanTOF, muonWEstdTOF = BeamDataAna.calculateTOF(muonLikeWeirdElectron[leadGlassID]["matchedHit0_TOF"])

BeamDataAna.plot1DHist(muonLikeWeirdElectron[leadGlassID]["matchedHit0_TOF"], 0.1, "TOF (ns)", "Occurences/0.1ns", "Muon-looking weird electrons\nfitted TOF: %.2f +/- %.2f"%(muonWEmeanTOF, muonWEstdTOF), "Time of flight of weird electrons", ax, fig)

_, pionWEmeanTOF, pionWEstdTOF = BeamDataAna.calculateTOF(pionLikeWeirdElectron[leadGlassID]["matchedHit0_TOF"])

BeamDataAna.plot1DHist(pionLikeWeirdElectron[leadGlassID]["matchedHit0_TOF"], 0.1, "TOF (ns)", "Occurences/0.1ns", "Pion-looking weird electrons\nfitted TOF: %.2f +/- %.2f"%(pionWEmeanTOF, pionWEstdTOF), "Time of flight of weird electrons", ax, fig)

_, muonmeanTOF, muonstdTOF = BeamDataAna.calculateTOF(BeamDataAna.muonArray[leadGlassID]["matchedHit0_TOF"])

BeamDataAna.plot1DHist(BeamDataAna.muonArray[leadGlassID]["matchedHit0_TOF"], 0.1, "TOF (ns)", "Occurences/0.1ns", "Muons\nfitted TOF: %.2f +/- %.2f"%(muonmeanTOF, muonstdTOF), "Time of flight of weird electrons", ax, fig)

_, pionmeanTOF, pionstdTOF = BeamDataAna.calculateTOF(BeamDataAna.pionArray[leadGlassID]["matchedHit0_TOF"])

BeamDataAna.plot1DHist(BeamDataAna.pionArray[leadGlassID]["matchedHit0_TOF"], 0.1, "TOF (ns)", "Occurences/0.1ns", "Pions\nfitted TOF: %.2f +/- %.2f"%(pionmeanTOF, pionstdTOF), "Time of flight of weird electrons", ax, fig)

ax.grid()

plt.savefig("../new_pdf_results/WeirdElectrons_%s_Run%i.pdf"%(title, BeamDataAna.runNumber))
plt.savefig("../new_png_results/WeirdElectrons_%s_Run%i.png"%(title, BeamDataAna.runNumber))

fig, ax = plt.subplots(1, 1, figsize = (16, 9))
title = "sumTS"

BeamDataAna.plot1DHist(muonLikeWeirdElectron[leadGlassID]["sumTS"], 10, "SumTS", "Occurences/0.1(PE)", "Muon-looking weird electrons", "Total charge deposited in TS for weird electrons", ax, fig)

BeamDataAna.plot1DHist(pionLikeWeirdElectron[leadGlassID]["sumTS"], 10, "sumTS", "Occurences/0.1(PE)", "Pion-looking weird electrons", "Total charge deposited in TS for weird electrons", ax, fig)

BeamDataAna.plot1DHist(BeamDataAna.muonArray[leadGlassID]["sumTS"], 10, "SumTS", "Occurences/10(PE)", "Muons", "Total charge deposited in TS for weird electrons", ax, fig)

BeamDataAna.plot1DHist(BeamDataAna.pionArray[leadGlassID]["sumTS"], 10, "sumTS", "Occurences/10(PE)", "Pions", "Total charge deposited in TS for weird electrons", ax, fig)

ax.grid()

plt.savefig("../new_pdf_results/WeirdElectrons_%s_Run%i.pdf"%(title, BeamDataAna.runNumber))
plt.savefig("../new_png_results/WeirdElectrons_%s_Run%i.png"%(title, BeamDataAna.runNumber))


for i in range(10):
    plt.close()

# branchList = ["matchedHit0_WindowIntPE", "WholeWaveformIntPE", "sumACT1", "sumDownstreamACTs", "delayBetweenCentreWindowAndPeakTime"]

branchList = ["WholeWaveformIntPE", "matchedHit0_WindowIntPE", "sumTS", "sumTSwindow2"]

# list_xlimits = [[0,  10], [0, 50], [0, 25], None, None]
list_xlimits = [[0,  30],[0,  30], None, None ]

# binsList = [0.5, 0.5, 0.3, 1, 2]
binsList = [0.5, 0.5, 10, 10]
 
for detectorName in ["PbGlass"]:
    for branchID in range(len(branchList)):
        fig, ax = plt.subplots(1, 1, figsize = (16, 9))
        branch = branchList[branchID]
        bin = binsList[branchID] 
        xlimits = list_xlimits[branchID]
        detector = config["channelNames"].index(detectorName)
        title = "%s_%s_MH0WIcut"%(detectorName, branch)

        # BeamDataAna.plot1DHist(muonLikeWeirdElectron[detector][branch], bin, "%s %s"%(detectorName, branch),  "Occurences/%.2f"%bin, "Muon-looking weird electrons", "%s in %s for weird electrons, muons, pions"%(branch, detectorName), ax, fig)

        # BeamDataAna.plot1DHist(pionLikeWeirdElectron[detector][branch], bin, "%s %s"%(detectorName, branch),  "Occurences/%.2f"%bin, "Pion-looking weird electrons", "%s in %s for weird electrons, muons, pions"%(branch, detectorName), ax, fig)

        BeamDataAna.plot1DHist(BeamDataAna.muonArray[detector][branch], bin, "%s %s"%(detectorName, branch),  "Occurences/%.2f"%bin, "Muons", "%s in %s for weird electrons, muons, pions"%(branch, detectorName), ax, fig)

        BeamDataAna.plot1DHist(BeamDataAna.pionArray[detector][branch], bin, "%s %s"%(detectorName, branch), "Occurences/%.2f"%bin, "Pions", "%s in %s for weird electrons, muons, pions"%(branch, detectorName), ax, fig, True)

        ax.grid()
        if xlimits != None:
            ax.set_xlim(xlimits)

        plt.savefig("../new_pdf_results/WeirdElectrons_%s_Run%i.pdf"%(title, BeamDataAna.runNumber))
        plt.savefig("../new_png_results/WeirdElectrons_%s_Run%i.png"%(title, BeamDataAna.runNumber))

        plt.close()