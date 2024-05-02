#This is the helper code doing the particleID_final analysis, class based instead of the previous messy multiple functions approach

import uproot as ur
import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd
import sys
import awkward as ak


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
        self.triggerScintillatorSelectionBool = config["triggerScintillatorSelectionBool"]
        self.triggerScintillatorSelectionValue = config["triggerScintillatorSelectionValue"]


        self.openedFile = None
        #looking only at the window integrated charge
        #arrays holding the dataframes
        self.arrayData = []

        #Good event selection, by default
        self.nCoincidenceSelectionPassed = None
        self.triggerScintillatorSelectionPassed = None

        #define some basic numbers
        self.numberOfTOFpmt = 4
        self.flag = -9999
        self.BoundsAreAvailable = False



    def openOneBranch(self, channelNumber):
        df_temp = self.openedFile[self.channelNames[channelNumber]].arrays(library = "pd")
        charge = pd.DataFrame(df_temp['WindowIntPE'].values.tolist())
        time = pd.DataFrame(df_temp['SignalTimeCorrected'].values.tolist())
        
        #Calculate the max number of hits to form sensible column name
        maxNbMatchedEvents = np.array(charge.shape)[1]
        maxNbPeakEvents = np.array(time.shape)[1]
        #matched events informations
        names_windowIntPE=["matchedHit%i_WindowIntPE"%i for i in range(maxNbMatchedEvents)]
        names_TOF=["matchedHit%i_TOF"%i for i in range(maxNbMatchedEvents)]
        names_WindowUpperBound = ["matchedHit%i_WindowUpperTime"%i for i in range(maxNbMatchedEvents)]
        names_WindowLowerBound = ["matchedHit%i_WindowLowerTime"%i for i in range(maxNbMatchedEvents)]
        names_WindowCentralTime = ["matchedHit%i_WindowCentralTime"%i for i in range(maxNbMatchedEvents)]
        names_WindowCentralTimeCorrected = ["matchedHit%i_WindowCentralTimeCorrected"%i for i in range(maxNbMatchedEvents)]

        #peak information
        names_IntPE=["peakHit%i_IntPE"%i for i in range(maxNbPeakEvents)]
        names_IntCharge=["peakHit%i_IntCharge"%i for i in range(maxNbPeakEvents)]
        names_SignalTimeCorrected=["peakHit%i_SignalTimeCorrected"%i for i in range(maxNbPeakEvents)]
        
        

        
        #Have a new column for TOF that will be filled either with the matched hit times or a flag at -9999
        #if the matching did not happen
        if 'SignalTimeMatchedTOF1' in df_temp.columns:
            print("SignalTimeMatched is available, using it to compute TOF")
            TOF = pd.DataFrame((ak.Array(df_temp['SignalTimeMatchedTOF1'].values.tolist())-ak.Array(df_temp['SignalTimeMatchedTOF0'].values.tolist())).tolist(), columns=names_TOF)
        
        else:
            print("SignalTimeMatched is not available, TOF will be set to %i"%self.flag)
            TOF = pd.DataFrame((df_temp['WindowIntPE']).values.tolist(), columns=names_TOF)
            for col in TOF.columns:
                TOF[col].values[:] = self.flag

        WindowIntPE = pd.DataFrame(df_temp['WindowIntPE'].values.tolist(), columns = names_windowIntPE)

        #peak information
        IntPE = pd.DataFrame(df_temp['IntPE'].values.tolist(), columns = names_IntPE)
        IntCharge = pd.DataFrame(df_temp['IntCharge'].values.tolist(), columns = names_IntCharge)
        SignalTimeCorrected = pd.DataFrame(df_temp['SignalTimeCorrected'].values.tolist(), columns = names_SignalTimeCorrected)

        #if we have information about the bounds
        if "WindowUpperTime" in df_temp.columns:
            self.BoundsAreAvailable = True
            WindowUpperBound = pd.DataFrame(df_temp['WindowUpperTime'].values.tolist(), columns = names_WindowUpperBound)
            WindowLowerBound = pd.DataFrame(df_temp['WindowLowerTime'].values.tolist(), columns = names_WindowLowerBound)
            WindowCentralTime = pd.DataFrame(df_temp['WindowCentralTime'].values.tolist(), columns = names_WindowCentralTime)
            WindowCentralTimeCorrected = pd.DataFrame(df_temp['WindowCentralTimeCorrected'].values.tolist(), columns = names_WindowCentralTimeCorrected)

            df_temp = pd.concat([df_temp, WindowUpperBound], axis = 1)
            df_temp = pd.concat([df_temp, WindowLowerBound], axis = 1)
            df_temp = pd.concat([df_temp, WindowCentralTime], axis = 1)
            df_temp = pd.concat([df_temp, WindowCentralTimeCorrected], axis = 1)

    
    
        df_temp = pd.concat([df_temp, WindowIntPE], axis = 1)
        df_temp = pd.concat([df_temp, TOF], axis =1)
        df_temp = pd.concat([df_temp, IntPE], axis = 1)
        df_temp = pd.concat([df_temp, IntCharge], axis = 1)
        df_temp = pd.concat([df_temp, SignalTimeCorrected], axis =1)


        return df_temp

    def openDataFile(self):
        #only open the file once
        if self.openedFile is None:
            print("Opening the data root file: %s"%self.dataFile)
            self.openedFile = ur.open(self.dataFile)
            availableChannels = [channel[:-2] for channel in self.openedFile.keys()]
            for channelNumber in range(len(self.channelNames)):
                print(f"Reading channel {self.channelNames[channelNumber]}...")
                if self.channelNames[channelNumber] in availableChannels:
                    df_temp = self.openOneBranch(channelNumber)
                    self.arrayData.append(df_temp)

                    print(f"... done\n", end = "", flush=True)
                else:
                    #instead make a copy of the previous branch and set everything 
                    #to a flag value
                    df_temp = self.openOneBranch(channelNumber-1) * 0 + -9999
                    self.arrayData.append(df_temp)
                    print(f"... skipping channel")
        else:
            raise Exception("The file already seems to be open.")

    
    def getDataFrameAllDetectors(self):
        print("The data has ", len(self.arrayData),"entries which hold the following column in the dataframe: ", self.arrayData[0].columns)
        return self.arrayData
    
    def getDataFrameDetector(self, detectorID):
        print("Returning the data frame for detector %i: %s"%(detectorID, self.channelNames[detectorID]))
        return self.arrayData[detectorID]
    
    def getColumnDataFrameDetector(self, branch_name, detectorID):
        print("Returning the column %s in the data frame for detector %i: %s"%(branch_name, detectorID, self.channelNames[detectorID]))
        return self.arrayData[detectorID][branch_name]
    
    def addBranchToDataFrameDetector(self, new_branch_name, detectorID, value):
        print("Adding the branch %s to the data frame for detector %i: %s"%(new_branch_name,detectorID, self.channelNames[detectorID]))
        self.arrayData[detectorID][new_branch_name] = value
        return 0 
    
    def addBranchToAllDetectors(self, new_branch_name, value):
        print("Adding the branch %s to the data frame for  all detectors"%(new_branch_name))
        for i in range(len(self.arrayData)):
            self.arrayData[i][new_branch_name] = value
        return 0 
    
    def addOperationBranchToAllDetectors(self, new_branch_name,  branch1_name, operation, branch2_name):
        for detectorID in np.arange(0, len(self.arrayData), 1):
            if operation=="+":
                res = self.getColumnDataFrameDetector(branch1_name, detectorID) + self.getColumnDataFrameDetector(branch2_name, detectorID)
            elif operation == '-':
                res = self.getColumnDataFrameDetector(branch1_name, detectorID) - self.getColumnDataFrameDetector(branch2_name, detectorID)
            elif operation == '*':
                res = self.getColumnDataFrameDetector(branch1_name, detectorID) * self.getColumnDataFrameDetector(branch2_name, detectorID)
            elif operation == '/':
                res = self.getColumnDataFrameDetector(branch1_name, detectorID) / self.getColumnDataFrameDetector(branch2_name, detectorID)
            else:
                raise Exception("Wrong operator, use '+' or '-' or '*' or '/' please.")
            self.addBranchToDataFrameDetector(new_branch_name, detectorID, res)

    
    
    # def plot2dHist(self, x, y):
    
    def plotHist1DfromData(self, array_columns, targetDetectors, plot_saving, normalise = False, MinBound = None, MaxBound = None, additionalLabel = None, yscale=None, nMaxColumns = 4, nbBins = 100, statsBin= False):
        print("Plotting 1D histogram of columns: ", array_columns, " normalisation is set to: ", normalise)
        fig, ax = plt.subplots(max(1, math.ceil(len(targetDetectors)/nMaxColumns)), nMaxColumns, figsize = (16, 9))
        range_min = []
        range_max = []
        binwidth = []

        if MinBound!=None:
            #in case we gave an int
            if type(MinBound) == int or type(MinBound) == float:
                MinBound = [MinBound for i in range(len(targetDetectors))]
            #in case we gave an array of the wrong size
            elif len(MinBound) != len(targetDetectors):
                print("\nWarning: you have not given the correct number of min Bounds, asusming some of them")
                MinBound = [MinBound[0] for i in range(len(targetDetectors))]

        if MaxBound!=None:
            #in case we gave an int
            if type(MaxBound) == int or type(MaxBound) == float:
                MaxBound = [MaxBound for i in range(len(targetDetectors))]
            #in case we gave an array of the wrong size
            elif len(MaxBound) != len(targetDetectors):
                print("\nWarning: you have not given the correct number of Max Bounds, asusming some of them")
                MaxBound = [MaxBound[0] for i in range(len(targetDetectors))]

        for variableID, column in enumerate(array_columns):
            #label the plot depending on the column we are looking at
            if (column.find('Time') != -1):
                xaxis_label = 'Time (ns)'
            elif (column.find('TOF') != -1):
                xaxis_label = 'Time of Flight (ns)'
            elif (column.find('PE') != -1):
                xaxis_label = 'Charge (PE)'
            else:
                xaxis_label = column
            for indexInTargetDetector, channelID in enumerate(targetDetectors):
                #channelID is the detector of interest, indexInTargetDetector is only for detectors we
                #want to plot
                if (variableID == 0):
                    #so the same variables will have the same bins 
                    if MinBound == None:
                        possible_min = min(self.arrayData[channelID][column]) * 0.8
                        if np.isnan(possible_min) or min(self.arrayData[channelID][column]) == -9999:
                            possible_min = -1

                        range_min.append(possible_min)
                    else:
                        possible_min = max(MinBound[indexInTargetDetector],min(self.arrayData[channelID][column]) * 0.8)
                        range_min.append(possible_min)
                    

                    if MaxBound == None:
                        possible_max = max(self.arrayData[channelID][column]) * 1.2
                        if np.isnan(possible_max) or max(self.arrayData[channelID][column]) == -9999:
                            possible_max = possible_min + 2
                        range_max.append(possible_max)
                    else:
                        range_max.append(min(MaxBound[indexInTargetDetector], max(self.arrayData[channelID][column]) * 1.2))

                    

                    binwidth.append((range_max[-1]-range_min[-1])/nbBins)

                if len(targetDetectors) == 1:
                    axis = ax
                elif len(targetDetectors)>nMaxColumns:
                    #if we have multiple subplots: need to put into correct one
                    plot_row = math.floor(indexInTargetDetector/nMaxColumns)
                    plot_column = indexInTargetDetector - (plot_row * nMaxColumns)

                    axis = ax[plot_row, plot_column]
                else:
                    axis = ax[indexInTargetDetector]

                if statsBin:
                    label = "%s: %i entries\n Mean: %.2f Std: %.2f"%(column, len(self.arrayData[channelID][column]), self.arrayData[channelID][column].mean(), self.arrayData[channelID][column].std())
                else:
                    label = "%s: %i entries"%(column, len(self.arrayData[channelID][column]))
                
                axis.hist(self.arrayData[channelID][column], bins = nbBins, range = (range_min[indexInTargetDetector], range_max[indexInTargetDetector]), label = label, density=normalise, histtype="step")
        #cosmetics titles etcs
        for indexInTargetDetector, channelID in enumerate(targetDetectors):
            if len(targetDetectors) == 1:
                axis = ax
            elif len(targetDetectors)>nMaxColumns:
                    plot_row = math.floor(indexInTargetDetector/nMaxColumns)
                    plot_column = indexInTargetDetector - (plot_row * nMaxColumns)
                    axis = ax[plot_row, plot_column]
            else:
                axis = ax[indexInTargetDetector]
            
            if yscale == 'log':
                axis.set_yscale('log')
                axis.set_ylim([0.5, None])
            
            if nMaxColumns >= len(targetDetectors):
                axis.legend(fontsize = 15)
            else:
                axis.legend(fontsize = 10)
            axis.set_title("%s"%(self.channelNames[channelID]), fontsize = 18, weight = "bold")
            axis.set_xlabel("%s"%xaxis_label, fontsize = 18)
            if binwidth[indexInTargetDetector] > 10e3:
                axis.set_ylabel("Occurences/%.3e"%(binwidth[indexInTargetDetector]), fontsize = 18)
            if binwidth[indexInTargetDetector] > 1:
                axis.set_ylabel("Occurences/%.1f"%(binwidth[indexInTargetDetector]), fontsize = 18)
            elif binwidth[indexInTargetDetector] > 0.1:
                axis.set_ylabel("Occurences/%.2f"%(binwidth[indexInTargetDetector]), fontsize = 18)
            elif binwidth[indexInTargetDetector] > 0.01:
                axis.set_ylabel("Occurences/%.3f"%(binwidth[indexInTargetDetector]), fontsize = 18)
            elif binwidth[indexInTargetDetector] > 0.001:
                axis.set_ylabel("Occurences/%.3e"%(binwidth[indexInTargetDetector]), fontsize = 18)
            axis.grid()

        if normalise:
            fig.suptitle('WCTE Beamtest - Run %i, p = %i MeV/c n = %s - density \n%s'%(self.runNumber, self.runMomentum, self.runRefractiveIndex, additionalLabel), fontsize=18, weight ='bold')
        else:
            fig.suptitle('WCTE Beamtest - Run %i, p = %i MeV/c \n%s'%(self.runNumber, self.runMomentum, additionalLabel), fontsize=18, weight ='bold')
        plt.tight_layout()
        plt.savefig("../new_png_results/%s_Run%i.png"%(plot_saving, self.runNumber))
        plt.savefig("../new_pdf_results/%s_Run%i.pdf"%(plot_saving, self.runNumber))
        

    def plotAllHist1DfromData(self, array_columns, plot_saving, normalise = False, MinBound = None, MaxBound = None, additionalLabel = None, yscale = None, nMaxColumns = 4, nbBins = 100):
        targetDetectors = range(len(self.arrayData))
        self.plotHist1DfromData(array_columns, targetDetectors, plot_saving, normalise, MinBound, MaxBound, additionalLabel, yscale, nMaxColumns, nbBins)
    
    def plotSomeDetectorHist1DfromData(self, array_columns, targetDetectors, plot_saving,  normalise = False, MinBound = None, MaxBound = None, additionalLabel = None, yscale = None, nMaxColumns = 4, nbBins = 100, statsBin = False):
        #if we only want to plot a few detectors
        if len(targetDetectors)<=4:
            nMaxColumns = len(targetDetectors)
        print(nMaxColumns)
        self.plotHist1DfromData(array_columns, targetDetectors, plot_saving, normalise, MinBound, MaxBound, additionalLabel, yscale, nMaxColumns, nbBins, statsBin)

    def plot2DScatterFromBranches(self, detector_x,  branch_name_x, detector_y, branch_name_y, detector_z = None, branch_name_z = None, units_x = "", units_y = "", units_z = "", additionalLabel = None, logz = False):
        fig, ax = plt.subplots(1, 1, figsize = (16, 9))
        x = self.arrayData[detector_x][branch_name_x]
        y = self.arrayData[detector_y][branch_name_y]

        range = [[max(x.mean()*0.6, min(x)*0.8), max(x)*1.04], [min(y)*0.8, max(y)*1.02]]
        if logz:
            hist = ax.scatter(x, y, c= self.arrayData[detector_z][branch_name_z], s = 3, norm = 'log', cmap = "viridis")
        else:
            hist = ax.scatter(x, y, c= self.arrayData[detector_z][branch_name_z],s = 3,  cmap = "viridis")

        fig.colorbar(hist, ax=ax, label = '%s %s %s'%(self.channelNames[detector_z],branch_name_z, units_z))
        ax.legend()
        ax.grid()
        # ax.set_xlim(range[0])
        # ax.set_ylim(range[1])
        ax.set_xlabel("%s %s %s"%(self.channelNames[detector_x], branch_name_x, units_x), fontsize=18)
        ax.set_ylabel("%s %s %s"%(self.channelNames[detector_y], branch_name_y, units_y), fontsize=18)
        fig.suptitle('WCTE Beamtest - Run %i, p = %i MeV/c n = %s\n%s'%(self.runNumber, self.runMomentum, self.runRefractiveIndex, additionalLabel), fontsize=18, weight ='bold')
        plt.savefig("../new_png_results/Scatter_%s_%s_Run%i_viridis.png"%(self.channelNames[detector_z],branch_name_z, self.runNumber))
        plt.savefig("../new_pdf_results/Scatter_%s_%s_Run%i_viridis.pdf"%(self.channelNames[detector_z], branch_name_z, self.runNumber))

    def plot2DHistFromBranches(self, detector_x,  branch_name_x, detector_y, branch_name_y, units_x = "", units_y = "", additionalLabel = None, logz = False, bins = [100, 100], range = None):
        fig, ax = plt.subplots(1, 1, figsize = (16, 9))
        x = self.arrayData[detector_x][branch_name_x]
        y = self.arrayData[detector_y][branch_name_y]

        if range == None:
            range = [[max(x.mean()*0.6, min(x)*0.8), max(x)*1.04], [min(y)*0.8, max(y)*1.02]]
        
        if logz:
            hist = ax.hist2d(x, y, bins = bins, norm = 'log', range = range, cmap = "viridis")
        else:
            hist = ax.hist2d(x, y, bins = bins, range = range, cmap = "viridis")

        clb = fig.colorbar(hist[3], ax=ax)
        clb.ax.set_title(label = 'Occurences/%.2f %s %.2f %s'%((range[0][1]-range[0][0])/bins[0], units_x, (range[1][1]-range[1][0])/bins[1], units_y), fontsize=18)
        # ax.legend()
        ax.grid()
        ax.set_xlim(range[0])
        ax.set_ylim(range[1])
        ax.set_xlabel("%s %s %s"%(self.channelNames[detector_x], branch_name_x, units_x), fontsize=18)
        ax.set_ylabel("%s %s %s"%(self.channelNames[detector_y], branch_name_y, units_y), fontsize=18)
        fig.suptitle('WCTE Beamtest - Run %i, p = %i MeV/c n = %s\n%s'%(self.runNumber, self.runMomentum, self.runRefractiveIndex, additionalLabel), fontsize=18, weight ='bold')
        plt.savefig("../new_png_results/Hist2d_%s%s_%s%s_Run%i_viridis.png"%(self.channelNames[detector_x],branch_name_x, self.channelNames[detector_y],branch_name_y, self.runNumber))
        plt.savefig("../new_png_results/Hist2d_%s%s_%s%s_Run%i_viridis.pdf"%(self.channelNames[detector_x],branch_name_x, self.channelNames[detector_y],branch_name_y, self.runNumber))


    def getAcceptancenCoincidenceCut(self):
        return sum(self.nCoincidenceSelectionPassed)/len(self.nCoincidenceSelectionPassed)

    def nCoincidenceSelection(self):
        if self.openedFile is None:
            self.openDataFile(self)

        #only if we want the nCoincidenceSelection
        if self.nCoincidenceSelectionBool and "SignalTimeMatchedTOF1" in self.arrayData[0].columns:
            print(f"Performing nCoincidenceSelection with n = {self.nCoincidenceSelectionValue} for run {self.runNumber}...")
            #the branches SignalTimeMatchedTOF1/0 have the format of the coincidence, we can just have a
            #read over them in one detector to have the general selection
            self.nCoincidenceSelectionPassed = self.arrayData[0]["SignalTimeMatchedTOF1"].map(len)==self.nCoincidenceSelectionValue
            self.applyCut(self.nCoincidenceSelectionPassed)
        
        elif self.nCoincidenceSelectionBool:
            print("\nCannot perform a selection on the number of coincidence found: no coincidence was processed for this file, please make sure that the branch SignalTimeMatchedTOF1 (and0) exists in the dataframe, might need to change the config_(no)Hodoscope.json...\n")
            #this will always be true, but it is usefulme sure keep the same output
            # self.nCoincidenceSelectionPassed = self.arrayData[0]["DigitimingOffset"]!= -9998
            # self.applyCut(self.nCoincidenceSelectionPassed)

        else:
             print(f"In the config file, checking the coincidence is set to {self.nCoincidenceSelectionBool}, we are therefore not processing the change", end="", flush=True)
            #  self.nCoincidenceSelectionPassed = self.arrayData[0]["DigitimingOffset"]!= -9998

    def applyCut(self, cut_boolean):
        print("... this cut keeps %i events (%.2f percent), relative to df before this cut\n"%(sum(cut_boolean), sum(cut_boolean)/len(cut_boolean) * 100))
        for i in range(len(self.arrayData)):
            self.arrayData[i] = self.arrayData[i][cut_boolean]
        return self.arrayData

    def getBranchList(self, detectorID):
        return self.arrayData[detectorID].columns
    
    def selectBasedOnCondition(self, detectorID, branch_name, operation, value):
        print(f"Performing condition based selection, keeping {self.channelNames[detectorID]} {branch_name} {operation} {value}...")
        if operation == ">=":
            passedSelection = self.arrayData[detectorID][branch_name] >= value
        if operation == ">":
            passedSelection = self.arrayData[detectorID][branch_name] > value
        if operation == "<":
            passedSelection = self.arrayData[detectorID][branch_name] < value
        if operation == "==":
            passedSelection = self.arrayData[detectorID][branch_name] == value
        self.applyCut(passedSelection)




            

    #
    # for channelNumber in range(len(channelNames)):
    #             self.arrayDataCharge[channelNumber] =










