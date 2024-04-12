#This is the helper code doing the particleID_final analysis, class based instead of the previous messy multiple functions approach

import uproot as ur
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys


class LowMomentumAnalysis:
    "This class sets up the analysis tools for the low momentum set-up"
    def __init__(self, config):
        self.channelNames = config["channelNames"]
        self.dataFile = config["dataFile"]
        self.runNumber = config["runNumber"]
        self.runMomentum = config["runMomentum"]
        self.runRefractiveIndex = config["runRefractiveIndex"]
        self.referenceSelectionCutFile = config["referenceSelectionCutFile"]
        self.batchMode = config["batchMode"]
        self.saveRootFile = config["saveRootFile"]

        #Selection of events based on the number of coincident hits
        self.nCoincidenceSelectionBool = config["nCoincidenceSelectionBool"]
        self.nCoincidenceSelectionValue = config["nCoincidenceSelectionValue"]

        #Selection of events based on the dE/dx of particles in the trigger scintillator
        self.triggerScintillatorSelectionBool = config["triggerScintillatorSelection"]
        self.triggerScintillatorSelectionValue = config["triggerScintillatorSelectionValue"]


        self.openedFile = None
        #looking only at the window integrated charge
        #arrays holding the dataframes
        self.arrayDataFramesCharge = []
        self.arrayDataFramesTime = []

        #Good event selection, by default
        self.nCoincidenceSelectionPassed = None
        self.triggerScintillatorSelectionPassed = None


        #Saving the efficiencies of the different cuts





    def openDataFile(self):
        #only open the file once
        if self.openedFile is None:
            print("Opening the data root file: %s"%self.dataFile)
            self.openedFile = ur.open(self.dataFile)
            for channelNumber in range(len(channelNames)):
                print(f"Reading channel {self.channelNames[channelNumber]}...")
                df_temp = self.openedFile[self.channelNames[channelNumber]].arrays(library = "pd")
                charge = pd.DataFrame(df_temp['WindowIntPE'].values.tolist())
                time = pd.DataFrame(df_temp['SignalTimeCorrected'].values.tolist())
                df_charge = pd.concat([df_temp, charge], axis =1)
                df_time = pd.concat([df_temp, time], axis = 1)

                self.arrayDataFramesCharge.append(df_charge)
                self.arrayDataFramesTime.append(df_time)
                print(f" ... done", end = "", flush=True)
        else:
            raise Exception("The file already seems to be open.")

    def nCoincidenceSelection(self):
        print()
        if self.openedFile is None:
            openDataFile(self)

        #only if we want the nCoincidenceSelection
        if self.nCoincidenceSelectionBool :
            print(f"Performing nCoincidenceSelection with n = {self.nCoincidenceSelectionValue} ...")
            #the branches SignalTimeMatchedTOF1/0 have the format of the coincidence, we can just have a
            #read over them in one detector to have the general selection
            self.nCoincidenceSelectionPassed = self.arrayDataFramesCharge[0]["SignalTimeMatchedTOF1"].map(len)==self.nCoincidenceSelectionValue

        else:
             print(f"In the config file, checking the coincidence is set to {self.nCoincidenceSelectionBool}, we are therefore not processing the change", end="", flush=True)
             self.nCoincidenceSelectionPassed =




    for channelNumber in range(len(channelNames)):
                self.arrayDataFramesCharge[channelNumber] =










