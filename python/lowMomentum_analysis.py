#This is the helper code doing the particleID_final analysis, class based instead of the previous messy multiple functions approach

import uproot as ur
import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd
import sys
import awkward as ak
from scipy.optimize import curve_fit
from collections import OrderedDict
import scipy.integrate as integrate
import os.path

import csv



#the mass of the particles in MeV/c
ms = OrderedDict()
ms['electron'] = 0.511
ms['muon'] = 105.658
ms['pion'] = 139.6
ms['proton'] = 938.3
ms['deuterium'] = 1876.

conv = 1.e9 #convert back to ns
c = 299792458 # m/s


def gaussian(x, amplitude, mean, std_dev):
    return amplitude * np.exp(-0.5 * ((x - mean) / std_dev) ** 2)

def fitGaussian(counts, bins):
    y = np.array(counts)
    x = (bins[1]-bins[0])/2 + np.array(bins)[:-1]
    #amplitude, mean, std
    initial_guess = [max(y),  x[y.argmax(axis = 0)], 1]
    # Fit the Gaussian function to the data points
    params, covariance = curve_fit(gaussian, x, y, p0=initial_guess)

    return params, covariance


def oneOverSqrtN(x, A, B):
        return A/np.sqrt(x) + B

class LowMomentumAnalysis:
    "This class sets up the analysis tools for either the low momentum or tagged gamma setups"
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
        self.TStotalChargeSelectionBool = config["TStotalChargeSelectionBool"]
        
        if self.TStotalChargeSelectionBool:
            self.TStotalChargeSelectionValue = config["TStotalChargeSelectionValue"]
            self.TSwindow2totalChargeSelectionValue = config["TSwindow2totalChargeSelectionValue"]
        
        else:
            self.TStotalChargeSelectionBool = None
            self.TSwindow2totalChargeSelectionValue = None

        self.isLowMomentum = config["isLowMomentum"]

        self.downstreamACTs = ["ACT2L", "ACT2R", "ACT3L", "ACT3R"]


        self.openedFile = None
        #looking only at the window integrated charge
        #arrays holding the dataframes
        self.arrayData = []
        self.totalNumberOfEvents = None

        #Good event selection, by default
        #Has the number of coincidence that we want
        self.nCoincidenceSelectionPassed = None
        #doesn't deposit a lot of energy in the trigger scintillator
        #except if if it is a fast particle 
        self.TStotalChargeSelectionPassed = None
        self.TSwindow2totalChargeSelectionPassed = None


        #define some basic numbers
        self.numberOfTOFpmt = 4
        self.flag = -9999
        self.WindowBoundsAreAvailable = False

        #time of flight related items
        self.distanceTOF1toTOF0 = config["distanceTOF1toTOF0"]
        self.electronMeanTOF = None

        #make the dataframes for the particles
        self.protonArray = None
        self.electronArray = None
        self.muonArray = None
        self.pionArray = None
        self.deuteriumArray = None

        #make the arry of bools, with the same shape as the original dataset with 
        #each particle type, try to apply the best possible cuts
        self.isProton = None
        self.isDeuterium = None
        self.isElectron = None
        self.isMuon = None
        self.isPion = None


        self.thereIsSecondWindow = False


        self.nProtons, self.nPions, self.nMuons, self.nDeuterium, self.nElectrons = None, None, None, None, None

        #get the selection items
        self.protonsTOFCut, self.protonsTOFMax = self.getProtonTOFSelectionBounds()
        

        if self.isLowMomentum:
            self.deuteriumTOFcut, self.deuteriumTOFmax = self.getDeuteriumTOFSelectionBounds()
    
        if self.isLowMomentum:
            self.ACTlinearA = config["ACTlinearA"]
            self.ACTlinearB = config["ACTlinearB"]

            self.horizontal_el = config["horizontal_el"]
            self.weirdElectronLGcut = config["weirdElectronLGcut"]

        #for deuterium/proton selection we only want to try to fit if we have a large number of events
        self.minNbEventsToFitTOF = 150

        #make dictionaries for the outputs
        #Tof
        self.dictTOFMean = {"electron": None, 
                            "muon": None,
                            "pion": None,
                            "proton": None,
                            "deuterium": None}

        self.dictTOFfitErrOnMean = {"electron": None, 
                            "muon": None,
                            "pion": None,
                            "proton": None,
                            "deuterium": None}
        
        self.dictTOFStd = {"electron": None, 
                            "muon": None,
                            "pion": None,
                            "proton": None,
                            "deuterium": None}
        
        self.dictTOFfitErrOnStd = {"electron": None, 
                            "muon": None,
                            "pion": None,
                            "proton": None,
                            "deuterium": None}
        

        #momentum
        self.dictMomentumMean = {"muon": None,
                            "pion": None,
                            "proton": None,
                            "deuterium": None}
        
        self.dictMomentumStatError = {"muon": None,
                            "pion": None,
                            "proton": None,
                            "deuterium": None}
        
        #Number of particle fitted in the TOF dimension (useful for protons)
        self.dictTOFfittedNparticles = {"electron": None, 
                            "muon": None,
                            "pion": None,
                            "proton": None,
                            "deuterium": None}


        if self.isLowMomentum:
            #for the selection, check which particles we are looking for
            self.ACTLowerCut = config["ACTLowerCut"]
            self.piMuBorderACT = config["piMuBorderACT"]
            
            if self.runMomentum > 300:
                self.particleNamesList = ["electron", "muon", "pion", "proton",  "deuterium"]
            else:
                self.particleNamesList = ["electron", "muon", "pion"]
                
        else:
            self.particleNamesList = ["electron", "proton"]


        self.outputFileName = config["outputFileName"]

        #Dean's number of particle in each bunch study to remove deadtime, for now empty 
        self.probaBunch = 1
            
    def setMinNbEventsToFitTOF(self,value):
        #set the minimum number of events that need to have passed any selection for the TOF-estimated momentum to be calculated, keeping a failry high number ensures that the statistical error is small
        self.minNbEventsToFitTOF = value

    def getProtonTOFSelectionBounds(self, p=None):
        #for plotting, it is useful to be able to get these bounds for any momentum. 
        if p == None:
            p = self.runMomentum 
        TOFresolution = 0.35 #ns
        fiveSigmaOfPionTOF = self.momentumToTOF(p, 'pion') + 5 * TOFresolution
        protonTOFminus3ns = self.momentumToTOF(p, 'proton') - 3 
        protonTOFplus3ns = self.momentumToTOF(p, 'proton') + 3 
        deuteriumTOFminus5ns = self.momentumToTOF(p, 'deuterium') - 5 
        if self.isLowMomentum:
            #worry about deuterium
            return max(fiveSigmaOfPionTOF, protonTOFminus3ns), min(protonTOFplus3ns, (deuteriumTOFminus5ns+protonTOFplus3ns)/2)

        else:
            #no need to worry about deuterium
            return max(fiveSigmaOfPionTOF, protonTOFminus3ns), protonTOFplus3ns

        
    
    def getDeuteriumTOFSelectionBounds(self, p=None):
        #for plotting, it is useful to be able to get these bounds for any momentum. 
        if p == None:
            p = self.runMomentum 
        TOFresolution = 0.35 #ns
        tenSigmaOfProtonTOF = self.momentumToTOF(p, 'proton') + 10 * TOFresolution
        deuteriumTOFminus5ns = self.momentumToTOF(p, 'deuterium') - 5 
        deutriumTOFplus5ns = self.momentumToTOF(p, 'deuterium') + 5 
        protonTOFplus5ns = self.momentumToTOF(p, 'proton') + 5 
        return max((deuteriumTOFminus5ns+protonTOFplus5ns)/2, deuteriumTOFminus5ns), deutriumTOFplus5ns
    
    
    def correctForTOFwrtElectron(self, tof):
        """"here we need to offset the TOF collected for each particle by the difference between the measured electron TOF and the physical elctron TOF, to account for any leftover cable length offset or calibration issues"""
        if self.electronMeanTOF == None:
            print("Electron TOF not measured, will have to measure later")
            TOF_array = self.getColumnDataFrameDetector("matchedHit0_TOF", 0, "electron") 
            _, self.electronMeanTOF, _ = self.calculateTOF(TOF_array)

        #if there is an offset in the electron flight time compared to the expected distance then we need to correct for that (cable lengths)
        electronTimeOffset = self.electronMeanTOF - self.distanceTOF1toTOF0 * conv / c

        print("The time offset between the electron mean TOF and the expected photon TOF is %.3fns"%electronTimeOffset)
        
        tof = tof - electronTimeOffset * self.distanceTOF1toTOF0

    
        return tof


    def TofToMomentum(self, tof, particle_label):
        #the tof needs to be the absolute flying time
        m = ms[particle_label]
        #account for the non-zero travel time of electrons, assume they travel at c
        # tof = tof - (self.distanceTOF1toTOF0) * conv /c
        
        p = m/np.sqrt(pow((tof) * c / (conv * self.distanceTOF1toTOF0), 2) - 1)
        return p
    
    def momentumToTOF(self, p, particle_label):
        m = ms[particle_label]
        TOF = (conv * self.distanceTOF1toTOF0) / c * np.sqrt((m*m)/(p*p) + 1)
        return TOF

    def openOneBranch(self, channelNumber):
        df_temp = self.openedFile[self.channelNames[channelNumber]].arrays(library = "pd")
        charge = pd.DataFrame(df_temp['WindowIntPE'].values.tolist())
        time = pd.DataFrame(df_temp['SignalTimeCorrected'].values.tolist())
        if 'Window2IntPE' in df_temp.columns:
            self.thereIsSecondWindow = True
            chargeWindow2 = pd.DataFrame(df_temp['Window2IntPE'].values.tolist())
            maxNbWindow2Events = np.array(chargeWindow2.shape)[1]


        
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

        if self.thereIsSecondWindow:
            #second windowInformation, might not have the same format as the first window due to different size
            names_window2IntPE=["matchedHit%i_Window2IntPE"%i for i in range(maxNbWindow2Events)]
            names_Window2UpperBound = ["matchedHit%i_Window2UpperTime"%i for i in range(maxNbWindow2Events)]
            names_Window2LowerBound = ["matchedHit%i_Window2LowerTime"%i for i in range(maxNbWindow2Events)]
            

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
            self.WindowBoundsAreAvailable = True
            WindowUpperBound = pd.DataFrame(df_temp['WindowUpperTime'].values.tolist(), columns = names_WindowUpperBound)
            WindowLowerBound = pd.DataFrame(df_temp['WindowLowerTime'].values.tolist(), columns = names_WindowLowerBound)
            WindowCentralTime = pd.DataFrame(df_temp['WindowCentralTime'].values.tolist(), columns = names_WindowCentralTime)
            WindowCentralTimeCorrected = pd.DataFrame(df_temp['WindowCentralTimeCorrected'].values.tolist(), columns = names_WindowCentralTimeCorrected)

            df_temp = pd.concat([df_temp, WindowUpperBound], axis = 1)
            df_temp = pd.concat([df_temp, WindowLowerBound], axis = 1)
            df_temp = pd.concat([df_temp, WindowCentralTime], axis = 1)
            df_temp = pd.concat([df_temp, WindowCentralTimeCorrected], axis = 1)

        if self.thereIsSecondWindow:
            Window2UpperBound = pd.DataFrame(df_temp['Window2UpperTime'].values.tolist(), columns = names_Window2UpperBound)
            Window2LowerBound = pd.DataFrame(df_temp['Window2LowerTime'].values.tolist(), columns = names_Window2LowerBound)
            Window2IntPE = pd.DataFrame(df_temp['Window2IntPE'].values.tolist(), columns = names_window2IntPE)

            df_temp = pd.concat([df_temp, Window2UpperBound], axis = 1)
            df_temp = pd.concat([df_temp, Window2LowerBound], axis = 1)
            df_temp = pd.concat([df_temp, Window2IntPE], axis = 1)


    
    
        df_temp = pd.concat([df_temp, WindowIntPE], axis = 1)
        df_temp = pd.concat([df_temp, TOF], axis =1)
        df_temp = pd.concat([df_temp, IntPE], axis = 1)
        df_temp = pd.concat([df_temp, IntCharge], axis = 1)
        df_temp = pd.concat([df_temp, SignalTimeCorrected], axis =1)


        return df_temp

    def openDataFile(self):
        #only open the file once
        #if the plot saving folders do not already exist, create them
        os.makedirs("../new_png_results", exist_ok=True)
        os.makedirs("../new_pdf_results", exist_ok=True)

        if self.openedFile is None:
            print("Opening the data root file: %s"%self.dataFile)
            self.openedFile = ur.open(self.dataFile)
            availableChannels = [channel[:-2] for channel in self.openedFile.keys()]
            for channelNumber in range(len(self.channelNames)):
                print(f"Reading channel {self.channelNames[channelNumber]}...")
                if self.channelNames[channelNumber] in availableChannels:
                    df_temp = self.openOneBranch(channelNumber)
                    self.arrayData.append(df_temp)
                    self.totalNumberOfEvents = len(df_temp)
                    print(f"... done\n", end = "", flush=True)
                else:
                    #instead make a copy of the previous branch and set everything 
                    #to a flag value
                    df_temp = self.openOneBranch(channelNumber-1) * 0 + -9999
                    self.arrayData.append(df_temp)
                    print(f"... skipping channel")
            
        else:
            raise Exception("The file already seems to be open.")
        
    def makeSumDownstreamACTs(self):
        "Sum of the waveform integrated charge (PE) for all the ACT2, 3 PMTs"
        if self.isLowMomentum:
            sumDownsteamACTs = 0
            for downstreamACTpmt in self.downstreamACTs:
                if downstreamACTpmt in self.channelNames:
                    detectorID = self.channelNames.index(downstreamACTpmt)
                    sumDownsteamACTs += self.getDataFrameDetector(detectorID)["matchedHit0_WindowIntPE"]
                else:
                    print(f"Careful, PMT {downstreamACTpmt} is not in your dataset")
                    return 0
            self.addBranchToAllDetectors("sumDownstreamACTs", sumDownsteamACTs)
            
        else:
            print("This is the tagged photon set-up, not low momentum, not calculating the sum of downstream ACT light")

    def getTotalNumberEvents(self):
        return self.totalNumberOfEvents

    def makeSumACT1(self):
        "Sum of the waveform integrated charge (PE) for all the ACT1 PMTs"
        sumACT1 = 0
        for ACTpmt in ["ACT1L", "ACT1R"]:
            if ACTpmt in self.channelNames:
                detectorID = self.channelNames.index(ACTpmt)
                sumACT1 += self.getDataFrameDetector(detectorID)["matchedHit0_WindowIntPE"]
            else:
                print(f"Careful, PMT {ACTpmt} is not in your dataset")
                return 0
        self.addBranchToAllDetectors("sumACT1", sumACT1)

    def makeSumTS(self):
        "Sum of the window integrated charge (PE) for all the Trigger Scintillators PMTs"
        sumTS = 0
        for TSpmt in ["TOF00", "TOF01", "TOF02", "TOF03", "TOF10", "TOF11", "TOF12", "TOF13"]:
            if TSpmt in self.channelNames:
                detectorID = self.channelNames.index(TSpmt)
                sumTS += self.getDataFrameDetector(detectorID)["matchedHit0_WindowIntPE"]
            else:
                print(f"Careful, PMT {TSpmt} is not in your dataset")
                return 0
        
        self.addBranchToAllDetectors("sumTS", sumTS)
    
    def makeSumTSwindow2(self):
        "Sum of the second window integrated charge (PE) for all the Trigger Scintillators PMTs"
        if self.thereIsSecondWindow:
            sumTSwindow2 = 0
            for TSpmt in ["TOF00", "TOF01", "TOF02", "TOF03", "TOF10", "TOF11", "TOF12", "TOF13"]:
                if TSpmt in self.channelNames:
                    detectorID = self.channelNames.index(TSpmt)
                    sumTSwindow2 += self.getDataFrameDetector(detectorID)["matchedHit0_Window2IntPE"]
                else:
                    print(f"Careful, PMT {TSpmt} is not in your dataset")
                    return 0
            
            self.addBranchToAllDetectors("sumTSwindow2", sumTSwindow2)

    
    def getDataFrameAllDetectors(self, particle = None):
        # print("The data has ", len(self.arrayData),"entries which hold the following column in the dataframe: ", self.arrayData[0].columns)
        try:
            if particle == None:
                return self.arrayData
            if particle == "proton":
                return self.protonArray
            if particle == "electron":
                return self.electronArray
            if particle == "muon":
                return self.muonArray
            if particle == "pion":
                return self.pionArray
            if particle == "deuterium":
                return self.deuteriumArray
        except:
            raise Exception("Wrong particle name, leave empty for all particles, or electron, muon, pion, proton, deuterium")
    
    def getDataFrameDetector(self, detectorID, particle = None):
        print("Returning the data frame for detector %i: %s"%(detectorID, self.channelNames[detectorID]))
        try:
            if particle == None:
                return self.arrayData[detectorID]
            if particle == "proton":
                return self.protonArray[detectorID]
            if particle == "electron":
                return self.electronArray[detectorID]
            if particle == "muon":
                return self.muonArray[detectorID]
            if particle == "pion":
                return self.pionArray[detectorID]
            if particle == "deuterium":
                return self.deuteriumArray[detectorID]
        except:
            raise Exception("Wrong particle name, leave empty for all particles, or electron, muon, pion, proton, deuterium")
    
    
    def getColumnDataFrameDetector(self, branch_name, detectorID, particle = None):
        print("Returning the column %s in the data frame for detector %i: %s"%(branch_name, detectorID, self.channelNames[detectorID]))
        try:
            if particle == None:
                return self.arrayData[detectorID][branch_name]
            if particle == "proton":
                return self.protonArray[detectorID][branch_name]
            if particle == "electron":
                return self.electronArray[detectorID][branch_name]
            if particle == "muon":
                return self.muonArray[detectorID][branch_name]
            if particle == "pion":
                return self.pionArray[detectorID][branch_name]
            if particle == "deuterium":
                return self.deuteriumArray[detectorID][branch_name]
        except:
            raise Exception("Wrong particle name, leave empty for all particles, or electron, muon, pion, proton, deuterium")
    
    def addBranchToDataFrameDetector(self, new_branch_name, detectorID, value, particle = None):
        print("Adding the branch %s to the data frame for detector %i: %s"%(new_branch_name,detectorID, self.channelNames[detectorID]))
        try:
            if particle == None:
                self.arrayData[detectorID][new_branch_name] = value
            if particle == "proton":
                self.protonArray[detectorID][new_branch_name] = value
            if particle == "electron":
                self.electronArray[detectorID][new_branch_name] = value
            if particle == "muon":
                self.muonArray[detectorID][new_branch_name] = value
            if particle == "pion":
                self.pionArray[detectorID][new_branch_name] = value
            if particle == "deuterium":
                self.deuteriumArray[detectorID][new_branch_name] = value
        except:
            raise Exception(f"Branch adding has failed, possibly because the particle name: {particle} was not correct")
        
    def getArraysParticleData(self, particleIndex):
        return self.getDataFrameAllDetectors(self.particleNamesList[particleIndex])
    
    def makeAllParticleSelection(self):
        
        for particleName in self.particleNamesList:
            print("Making selection for %s"%particleName)
            
            if particleName == 'electron':
                self.makeElectronSelection()

            elif particleName == 'proton':  
                self.makeProtonSelection()

            elif particleName == 'muon':
                self.makeMuonSelection()

            elif particleName == 'pion':
                self.makePionSelection()

            elif particleName == 'deuterium':
                self.makeDeuteriumSelection()
        

    def measureMomentumUsingTOF(self, binWidth = 0.5):
        """For each particle type, if we have identified it, fit the TOF peak to extract the momentum, only do it with protons if the other ones haven't been identified yet"""
        fig, ax = plt.subplots(1, 1, figsize = (16, 9))
        for particle in self.particleNamesList:
            if self.getDataFrameAllDetectors(particle) == None:
                if particle == "proton":
                    self.makeProtonSelection()
                elif particle == "electron":
                    self.makeElectronSelection()
            if self.getDataFrameAllDetectors(particle) != None:
                #for low momentum set-up sometimes we do not have any sensible number of particle, when there is not enough e.g. deuterium or protons at low momentum 
                if len(self.getDataFrameDetector(0, particle)) > self.minNbEventsToFitTOF or not(self.isLowMomentum):
                    calc_mom, stat_err_mom = self.calculateMomentumUsingTOF(particle, fig, ax, binWidth)
                    self.dictMomentumMean[particle] = calc_mom
                    self.dictMomentumStatError[particle] = stat_err_mom
                    
        ax.legend(fontsize = 18)
        ax.grid()
        particleName = "AllParticles"
        ax.set_yscale('log')
        ax.set_ylim([0.5, None])
        ax.set_xlabel("Time of Flight (ns)", fontsize=15)
        ax.set_ylabel("Occurences/%.2f ns"%binWidth, fontsize=15)
        ax.set_title("Particle: %s"%particleName, fontsize=25)
        ax.tick_params(axis='both', which='major', labelsize=15)
        fig.suptitle('WCTE Beamtest - Run %i, p = %i MeV/c'%(self.runNumber, self.runMomentum), fontsize=22, weight ='bold')
        fig.savefig("../new_png_results/TOFHisto_%s_Run%i.png"%(particleName, self.runNumber))
        fig.savefig("../new_pdf_results/TOFHisto_%s_Run%i.pdf"%(particleName, self.runNumber))
        
                
    def calculateTOF(self, TOF_array, binWidth = 0.5, output_error_std = False):
        
        nBins = int((np.ceil(TOF_array.max() - np.floor(TOF_array.min())))/binWidth)
        
        counts, bins = np.histogram(TOF_array, nBins)

        params, covariance = fitGaussian(counts, bins)

        amplitude, mean, std = params[0], params[1], abs(params[2])

        if not(output_error_std):
            return amplitude, mean, std
        else:
            return amplitude, mean, std, covariance[2][2]
    
    def calculateMomentumUsingTOF(self, particleName, fig = None, ax = None, binWidth = 0.5):
        """here calculate the momentum by fitting the histogram of the particle TOF momentum, note we can use detector 0 because all have it now and it is a straight copy"""

        if "matchedHit0_TOF" not in (self.getBranchList(0, particleName)): 
            raise Exception("The matchedHit0_TOF branch is not available, make sure you are using the coincidence in your pre-processing of the data so that the TOF measurement is accurate ")
            
        TOF_array = self.getColumnDataFrameDetector("matchedHit0_TOF", 0, particleName) 
        nBins = int((np.ceil(TOF_array.max() - np.floor(TOF_array.min())))/binWidth)

        fig_particle, ax_particle = plt.subplots(1, 1, figsize = (16, 9))
        
        counts, bins, _ = ax_particle.hist(TOF_array, nBins, histtype = 'step', label = '# events passing the selection: %i\n(%.1f percent of all events)'%(len(TOF_array), (len(TOF_array)/self.totalNumberOfEvents) * 100))

        if ax != None:
            ax.hist(TOF_array, nBins, histtype = 'step', label = '# events passing the %s selection: %i\n(%.1f percent of all events)'%(particleName, len(TOF_array), (len(TOF_array)/self.totalNumberOfEvents) * 100))


        params, covariance = fitGaussian(counts, bins)

        amplitude, mean, std = params[0], params[1], abs(params[2])

        #need to correct for cable length offsets, especially important for TG set-up
        if particleName != "electron":
            mean = self.correctForTOFwrtElectron(mean)

        self.dictTOFMean[particleName] = mean
        self.dictTOFfitErrOnMean[particleName] = covariance[1][1]
        self.dictTOFStd[particleName] = std
        self.dictTOFfitErrOnStd[particleName] = covariance[2][2]

        p_pred = self.TofToMomentum(mean, particleName)

        n_events_inPeak = integrate.quad(gaussian, 0, 50, args = (amplitude, mean, std))

        n_events_inPeak = n_events_inPeak[0]/binWidth

        self.dictTOFfittedNparticles[particleName] = n_events_inPeak
        #error on the p_error is (proportional to) sigma_TOF/sqrt(n) where n is the number of events in the peak 
        p_error = (p_pred * std/mean)  * (p_pred /(ms[particleName]) * (mean/((self.distanceTOF1toTOF0 * conv)/c))) ** 2 / np.sqrt(n_events_inPeak) 
        #/ np.sqrt(n_events_inPeak) 

        print(f"Predicted momentum for {particleName} is {p_pred} MeV/c, true is {self.runMomentum} MeV/c.")

        x = np.linspace(bins[0], bins[-1], 1000)
        ax_particle.plot(x, gaussian(x, *params), '--', label = '# %s in fitted peak %.1f \n (%.1f percent of all events)\nTOF = %.2f +/- %.2f ns\nMomentum: %.1f +/- %.1f(stat) MeV/c'%(particleName, n_events_inPeak, (n_events_inPeak/self.totalNumberOfEvents) * 100, mean, std, p_pred, p_error))

        if ax != None:
            if particleName != 'electron':
                ax.plot(x, gaussian(x, *params), '--', label = 'TOF = %.2f +/- %.2f ns\nMomentum: %.1f +/- %.1f(stat) MeV/c'%(mean, std, p_pred, p_error))
            else:
               ax.plot(x, gaussian(x, *params), '--', label = 'TOF = %.2f +/- %.2f ns'%(mean, std)) 

        ax_particle.legend(fontsize = 18)
        ax_particle.grid()
        ax_particle.set_xlabel("Time of Flight (ns)", fontsize=15)
        ax_particle.set_ylabel("Occurences/%.2f ns"%binWidth, fontsize=15)
        ax_particle.set_title("Particle: %s"%particleName, fontsize=25)
        ax_particle.tick_params(axis='both', which='major', labelsize=15)
        fig_particle.suptitle('WCTE Beamtest - Run %i, p = %i MeV/c'%(self.runNumber, self.runMomentum), fontsize=22, weight ='bold')
        fig_particle.savefig("../new_png_results/TOFHisto_%s_Run%i.png"%(particleName, self.runNumber))
        fig_particle.savefig("../new_pdf_results/TOFHisto_%s_Run%i.pdf"%(particleName, self.runNumber))

        return p_pred, p_error
    

    def plotBranchHistForAllParticles(self, detectorID, branch_name, binwidth = 1, logScale = False, lim = None):
        
        """Plot for all the particles the branch in the given detector """
        fig, ax = plt.subplots(1, 1, figsize = (16, 9)) 
        ax.grid()

        if lim !=None:
            ax.set_xlim(lim)
        
        for particle in self.particleNamesList:
            if branch_name in self.getBranchList(detectorID, particle):
                self.plot1DHist(self.getColumnDataFrameDetector(branch_name, detectorID, particle), binwidth, "%s %s"%(self.channelNames[detectorID], branch_name), "Number of occurences", particle, "%s_%s_allParticles"%(self.channelNames[detectorID], branch_name), ax, fig)
            else:
                break
        if logScale:
            if branch_name in self.getBranchList(detectorID):
                self.plot1DHist(self.getColumnDataFrameDetector(branch_name, detectorID), binwidth,"%s %s"%(self.channelNames[detectorID], branch_name), "Number of occurences", "All particles", "%s_%s_allParticles"%(self.channelNames[detectorID], branch_name), ax, fig, True)

        else:
            if branch_name in self.getBranchList(detectorID):
                self.plot1DHist(self.getColumnDataFrameDetector(branch_name, detectorID), binwidth, "%s %s"%(self.channelNames[detectorID], branch_name), "Number of occurence", "All particles", "%s_%s_allParticles"%(self.channelNames[detectorID], branch_name), ax, fig)
        
        
        

    def measureElTOFresolutionFunctionOfTScharge(self, nBins = 10):
        """Plot the standard deviation of the fitted TOF distribution for the electron sample after having split it accounding to bins (of equal number of events) in the Trigger Scintillator summed charge (given as WholeWaveformIntPE)"""
        electronSample = self.getColumnDataFrameDetector("sumTS", 0, "electron")
        ESsorted = sorted(electronSample)
        nbEventsPerBin = int(len(ESsorted)/nBins)
        array_mean, array_std, array_meanX, array_stdX = [], [], [], []

        fig, ax1 = plt.subplots(1, 1, figsize = (16, 9))

        for bin in range(nBins):
            binMinIndex = bin * nbEventsPerBin
            binMaxIndex = min((bin+1) * nbEventsPerBin, len(ESsorted))
            valuesInBin = np.array(ESsorted[binMinIndex:binMaxIndex])
            binEdgeMin = min(valuesInBin)
            binEdgeMax = max(valuesInBin)
            
            meanBinValue = valuesInBin.mean()
            stdBin = valuesInBin.std()
            medianBinValue = valuesInBin[int(len(valuesInBin)/2)]

            print("Bin %i ranges between %.2f (ID: %i) and %.2f(ID: %i), holds %i values with a mean of %.2f and an std of %.2f with a median of %.2f"%(bin, binEdgeMin, binMinIndex, binEdgeMax,binMaxIndex, len(valuesInBin), meanBinValue, stdBin, medianBinValue))

            isAboveLowerBinBound = (electronSample>=binEdgeMin)
            isBelowUpperBinBound = (electronSample<=binEdgeMax)
            isInBin = isAboveLowerBinBound & isBelowUpperBinBound
            e_sampleInBin = self.makeNewDataFrameFromSelection(self.getDataFrameAllDetectors("electron"), isInBin)
            _, TOF_mean, TOF_std, TOF_std_err = self.calculateTOF(e_sampleInBin[0]["matchedHit0_TOF"], 0.5, True)

            print("Mean time of flight in bin: %.2f std: %.3f"%(TOF_mean, TOF_std))

            array_mean.append(TOF_mean)
            array_std.append(TOF_std)
            array_meanX.append(meanBinValue)
            array_stdX.append(stdBin)

        plot_saving = 'stdElectronsPerTSbin_WholeWaveformIntPE'
        ax1.errorbar(array_meanX, array_std, xerr = array_stdX, yerr = TOF_std_err, fmt = 'x')
        
        params, covariance = curve_fit(oneOverSqrtN, array_meanX[:-1], array_std[:-1], p0 = [array_mean[0], 0.3])

        print("Run %i covariance matrix for a fit A/sqrt(n) + b:\n"%self.runNumber, covariance)

        x = np.linspace(array_meanX[0], array_meanX[-2], 100)
        ax1.plot(np.array(x), np.array(oneOverSqrtN(x, params[0], params[1])), 'r--', label = 'Fit: %.3f/sqrt(x) + %.3f'%(params[0],params[1]) )
        ax1.legend(fontsize = 18)

        ax1.set_xlabel("Mean TS whole waveform int charge in bin +/- std (PE)", fontsize = 18)
        ax1.set_ylabel("Std of electron time of flight (ns)", fontsize = 18)
        ax1.grid()
        fig.suptitle('WCTE Beamtest - Run %i, p = %i MeV/c n = %s'%(self.runNumber, self.runMomentum, self.runRefractiveIndex), fontsize=18, weight ='bold')
        ax1.tick_params(axis='both', which='major', labelsize=15)
        fig.savefig("../new_pdf_results/%s_%iTSchargeBins_Run%i.pdf"%(plot_saving, nBins, self.runNumber))
        fig.savefig("../new_png_results/%s_%iTSchargeBins_Run%i.png"%(plot_saving, nBins, self.runNumber))


    
    def addBranchToAllDetectors(self, new_branch_name, value, particle = None):
        print("Adding the branch %s to the data frame for  all detectors"%(new_branch_name))
        for dectector in range(len(self.arrayData)):
            self.addBranchToDataFrameDetector(new_branch_name, dectector, value, particle)
        return 0 
    
    def addOperationBranchToAllDetectors(self, new_branch_name,  branch1_name, operation, branch2_name, particle = None):
        dataFrameOfInterest = self.getDataFrameAllDetectors(particle)
        for detectorID in np.arange(0, len(dataFrameOfInterest), 1):
            if operation=="+":
                res = self.getColumnDataFrameDetector(branch1_name, detectorID, particle) + self.getColumnDataFrameDetector(branch2_name, detectorID, particle)
            elif operation == '-':
                res = self.getColumnDataFrameDetector(branch1_name, detectorID, particle) - self.getColumnDataFrameDetector(branch2_name, detectorID, particle)
            elif operation == '*':
                res = self.getColumnDataFrameDetector(branch1_name, detectorID, particle) * self.getColumnDataFrameDetector(branch2_name, detectorID, particle)
            elif operation == '/':
                res = self.getColumnDataFrameDetector(branch1_name, detectorID, particle) / self.getColumnDataFrameDetector(branch2_name, detectorID, particle)
            else:
                raise Exception("Wrong operator, use '+' or '-' or '*' or '/' please.")
            self.addBranchToDataFrameDetector(new_branch_name, detectorID, res)


    
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
        self.plotHist1DfromData(array_columns, targetDetectors, plot_saving, normalise, MinBound, MaxBound, additionalLabel, yscale, nMaxColumns, nbBins, statsBin)

    def plot1DHist(self, data1D, binWidth = 1, xlabel = 'x', ylabel = 'y', legend = None, title = None, ax = None, fig = None, ylogscale = False):

        bins = int((max(data1D)-min(data1D))/binWidth)
        if ax == None or fig == None:
            fig, ax = plt.subplots(1, 1, figsize = (16, 9))
        if legend != None:
            ax.hist(data1D, bins = bins, label = '%s\n%i events'%(legend, len(data1D)), histtype= 'step')
            ax.legend(fontsize = 18)
        else:
            ax.hist(data1D, bins = bins, histtype= 'step', label = '%i events'%(len(data1D)))
            ax.legend(fontsize = 18)

        # ax.grid()
        ax.set_xlabel("%s"%xlabel, fontsize = 18)
        ax.set_ylabel("%s/%.2f"%(ylabel, binWidth), fontsize = 18)
        ax.tick_params(axis='both', which='major', labelsize=15)
        if title != None:
            fig.suptitle("%s Run %i - %.0f MeV/c"%(title, self.runNumber, self.runMomentum), fontsize = 20, weight = 'bold')
        else:
            fig.suptitle("Run %i - %.0f MeV/c"%(self.runNumber, self.runMomentum), fontsize = 20, weight = 'bold')

        if ylogscale:
            ax.set_yscale('log')

        plt.savefig("../new_pdf_results/Hist1d_%s_Run%i.pdf"%(title, self.runNumber))
        plt.savefig("../new_png_results/Hist1d_%s_Run%i.png"%(title, self.runNumber))



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
        plt.savefig("../new_png_results/Scatter_%s_%s_Run%i.png"%(self.channelNames[detector_z],branch_name_z, self.runNumber))
        plt.savefig("../new_pdf_results/Scatter_%s_%s_Run%i.pdf"%(self.channelNames[detector_z], branch_name_z, self.runNumber))

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
        plt.savefig("../new_png_results/Hist2d_%s%s_%s%s_Run%i.png"%(self.channelNames[detector_x],branch_name_x, self.channelNames[detector_y],branch_name_y, self.runNumber))
        plt.savefig("../new_pdf_results/Hist2d_%s%s_%s%s_Run%i.pdf"%(self.channelNames[detector_x],branch_name_x, self.channelNames[detector_y],branch_name_y, self.runNumber))


    def getAcceptancenCoincidenceCut(self):
        return sum(self.nCoincidenceSelectionPassed)/len(self.nCoincidenceSelectionPassed)
    
    def getArrayData(self):
        return self.arrayData

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

    def TStotalChargeSelection(self):
        """"Apply a cut on the minimum energy deposited in the Trigger scintillator (sum of all of the windowIntPE in TOFxy PMTs) to get rid of the poor coincidence matches to improve the TOF resolution"""
        if self.openedFile is None:
            self.openDataFile(self)
        
        if self.TStotalChargeSelectionBool:
            if "sumTS" not in self.getBranchList(0):
                self.makeSumTS()
            depositsEnoughEnergy = self.getSelectionBasedOnCondition(0, "sumTS", ">=", self.TStotalChargeSelectionValue)
            isFast = self.getSelectionBasedOnCondition(0, "matchedHit0_TOF", ">=", self.protonsTOFCut)

            self.TStotalChargeSelectionPassed =  depositsEnoughEnergy
            #apply that cut to the main dataframe with all of the data
            self.applyCut(self.TStotalChargeSelectionPassed)

        if self.thereIsSecondWindow:
            if "sumTSwindow2" not in self.getBranchList(0):
                self.makeSumTSwindow2()
            depositsNotTooMuchEnergy = self.getSelectionBasedOnCondition(0, "sumTSwindow2", "<=", self.TSwindow2totalChargeSelectionValue)
            isFast = self.getSelectionBasedOnCondition(0, "matchedHit0_TOF", ">=", self.protonsTOFCut)
            self.TSwindow2totalChargeSelectionPassed = depositsNotTooMuchEnergy | isFast

            self.applyCut(self.TSwindow2totalChargeSelectionPassed)


    def applyCut(self, cut_boolean):
        print("... this cut keeps %i events (%.2f percent), relative to df before this cut\n"%(sum(cut_boolean), sum(cut_boolean)/len(cut_boolean) * 100))
        for i in range(len(self.arrayData)):
            self.arrayData[i] = self.arrayData[i][cut_boolean]
        return self.arrayData
    
    def makeNewDataFrameFromSelection(self, initial_df_array, cut_boolean):
        print("... this cut keeps %i events (%.2f percent), relative to df before this cut\n"%(sum(cut_boolean), sum(cut_boolean)/len(cut_boolean) * 100))
        final_df_array = initial_df_array.copy()
        for detector in range(len(initial_df_array)):
            final_df_array[detector] = initial_df_array[detector][cut_boolean]
        return final_df_array


    def getBranchList(self, detectorID, particle = None):
        "Default we get the main datafram but then we can have others"
        if particle == None:
            return self.arrayData[detectorID].columns
        if particle == "proton":
            if self.protonArray == None:
                raise Exception("proton branch not yet defined, run makeAllParticleSelection or makeProtonSelection first")
            return self.protonArray[detectorID].columns
        if particle == "electron":
            if self.electronArray == None:
                raise Exception("electron branch not yet defined, run makeAllParticleSelection or makeElectronSelection first")
            return self.electronArray[detectorID].columns
        if particle == "muon":
            if self.muonArray == None:
                raise Exception("muon branch not yet defined, run makeAllParticleSelection or makeMuonSelection first")
            return self.muonArray[detectorID].columns
        if particle == "pion":
            if self.pionArray == None:
                raise Exception("pion branch not yet defined, run makeAllParticleSelection or makePionSelection first")
            return self.pionArray[detectorID].columns
        if particle == "deuterium":
            if self.deuteriumArray == None:
                raise Exception("deuterium branch not yet defined, run makeAllParticleSelection or makeDeuteriumSelection first")
            return self.deuteriumArray[detectorID].columns
        
    
    def cutBasedOnCondition(self, detectorID, branch_name, operation, value):
        """directly applies the cut on the main datafram, a bit dangerous, not recommended to use"""
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

    def getSelectionBasedOnCondition(self, detectorID, branch_name, operation, value, initial_df=None):
        "Get for a given array the array of bools being true if the event passes selection, false otherwise"
        if initial_df == None:
            initial_df = self.arrayData

        print(f"Performing condition based selection, keeping {self.channelNames[detectorID]} {branch_name} {operation} {value}...")
        if operation == ">=":
            passedSelection = initial_df[detectorID][branch_name] >= value
        if operation == "<=":
            passedSelection = initial_df[detectorID][branch_name] <= value
        if operation == ">":
            passedSelection = initial_df[detectorID][branch_name] > value
        if operation == "<":
            passedSelection = initial_df[detectorID][branch_name] < value
        if operation == "==":
            passedSelection = initial_df[detectorID][branch_name] == value
        return passedSelection

    def makeProtonSelection(self, timingCut = None, maxTimingCut = None):
        if timingCut == None:
            #use the default timing cut from the config file
            timingCut = self.protonsTOFCut
            
        if  maxTimingCut == None:
            maxTimingCut = self.protonsTOFMax 
        print(f"Selecting as protons all the events that have TOF larger than {timingCut} and smaller than {maxTimingCut}")
        if self.protonArray != None:
            print("This is overwriting the existing proton selection")
        isSlow = self.getSelectionBasedOnCondition(0, "matchedHit0_TOF", '>=', timingCut)
        isNotDeuterium = self.getSelectionBasedOnCondition(0, "matchedHit0_TOF", '<=', maxTimingCut)
        self.isProton = isSlow & isNotDeuterium

        if self.isLowMomentum:
            isNotElectron = self.getSelectionBasedOnCondition(0, "sumDownstreamACTs", "<=", self.horizontal_el)
            self.isProton = self.isProton & isNotElectron
        self.protonArray = self.makeNewDataFrameFromSelection(self.arrayData, self.isProton)

    def makeDeuteriumSelection(self,  timingCut = None, maxTimingCut = None):
        if  maxTimingCut == None:
            maxTimingCut = self.deuteriumTOFmax
        if  timingCut == None:
            timingCut = self.deuteriumTOFcut
        print(f"Selecting as deuterium all the events that have TOF larger than {maxTimingCut}")
        isFasterThanProtons = self.getSelectionBasedOnCondition(0, "matchedHit0_TOF", '>', timingCut)
        isNotTooFast = self.getSelectionBasedOnCondition(0, "matchedHit0_TOF", '<', maxTimingCut)

        self.isDeuterium = isFasterThanProtons & isNotTooFast

        
        if self.isLowMomentum:
            #deuterium events should not deposit a lot of light in 
            #the downstream ACT, these could be electrons with coincidence issues. 
            isNotElectron = self.getSelectionBasedOnCondition(0, "sumDownstreamACTs", "<=", self.horizontal_el)
            self.isDeuterium = self.isDeuterium & isNotElectron

        self.deuteriumArray = self.makeNewDataFrameFromSelection(self.arrayData, self.isDeuterium)

    def makePionSelection(self):
        if "sumDownstreamACTs" not in self.getBranchList(0):
                self.makeSumDownstreamACTs()

        print(f"Selecting as pions all the events that have sumDownstreamACTs larger than {self.ACTLowerCut} and smaller than {self.piMuBorderACT}")
        #I am not sure whether we should keep that, we are actually expecting to have 0 charge there, maybe we are loosing quite a bit of pions, but actually it was less than percentage
        isAboveACTLower = self.getSelectionBasedOnCondition(0, "sumDownstreamACTs", ">=", self.ACTLowerCut)

        isBelowPiMuCutLine = self.getSelectionBasedOnCondition(0, "sumDownstreamACTs", "<", self.piMuBorderACT)

        isNotWeirdElectron = self.getSelectionBasedOnCondition(self.channelNames.index("PbGlass"), "matchedHit0_WindowIntPE", "<", self.weirdElectronLGcut)
    
        isSlow = self.getSelectionBasedOnCondition(0, "matchedHit0_TOF", "<", self.protonsTOFCut)

        ACTSelection = isAboveACTLower & isBelowPiMuCutLine

        isPotentialPion = ACTSelection & isSlow

        self.isPion  = isPotentialPion & isNotWeirdElectron
        # self.nPions = sum(self.isPion)
        self.pionArray = self.makeNewDataFrameFromSelection(self.arrayData, self.isPion)

        self.isPionLikeWeirdElectron = isPotentialPion & np.logical_not(isNotWeirdElectron)

        self.pionLikeWeirdElectronArray = self.makeNewDataFrameFromSelection(self.arrayData, self.isPionLikeWeirdElectron)


    def makeMuonSelection(self):
        if "sumACT1" not in self.getBranchList(0):
                    self.makeSumACT1()
        if "sumDownstreamACTs" not in self.getBranchList(0):
                self.makeSumDownstreamACTs()

        isNotWeirdElectron = self.getSelectionBasedOnCondition(self.channelNames.index("PbGlass"), "matchedHit0_WindowIntPE", "<", self.weirdElectronLGcut)

        isBelowMuElCut = self.getSelectionBasedOnCondition(0, "sumDownstreamACTs", "<=", self.getDataFrameDetector(0)["sumACT1"] * self.ACTlinearA + self.ACTlinearB)

        isBelowHorizontalElCut = self.getSelectionBasedOnCondition(0, "sumDownstreamACTs", "<=", self.horizontal_el)

        isAbovePiMuCutLine = self.getSelectionBasedOnCondition(0, "sumDownstreamACTs", ">=", self.piMuBorderACT)

        isSlow = self.getSelectionBasedOnCondition(0, "matchedHit0_TOF", "<", self.protonsTOFCut)

        isNotElectron = isBelowMuElCut | isBelowHorizontalElCut

        isMuonOrElectron = isSlow & isAbovePiMuCutLine

        isPotentialMuon = isMuonOrElectron & isNotElectron

        self.isMuon = isPotentialMuon & isNotWeirdElectron

        self.muonArray = self.makeNewDataFrameFromSelection(self.arrayData, self.isMuon)

        self.isMuonLikeWeirdElectron = isPotentialMuon & np.logical_not(isNotWeirdElectron)

        self.muonLikeWeirdElectronArray = self.makeNewDataFrameFromSelection(self.arrayData, self.isMuonLikeWeirdElectron)


    def makeElectronSelection(self):
        if self.isLowMomentum:
            if "sumACT1" not in self.getBranchList(0):
                    self.makeSumACT1()
            if "sumDownstreamACTs" not in self.getBranchList(0):
                    self.makeSumDownstreamACTs()

            isAboveMuElCut = self.getSelectionBasedOnCondition(0, "sumDownstreamACTs", ">", self.getDataFrameDetector(0)["sumACT1"] * self.ACTlinearA + self.ACTlinearB)

            isAboveHorizontalElCut = self.getSelectionBasedOnCondition(0, "sumDownstreamACTs", ">", self.horizontal_el)

            isSlow = self.getSelectionBasedOnCondition(0, "matchedHit0_TOF", "<", self.protonsTOFCut)

            self.isElectron = isAboveMuElCut | isAboveHorizontalElCut

            self.isElectron = self.isElectron & isSlow

            self.electronArray = self.makeNewDataFrameFromSelection(self.arrayData, self.isElectron)

        else:
            self.isElectron = self.getSelectionBasedOnCondition(0, "matchedHit0_TOF", "<", self.protonsTOFCut)
            self.electronArray = self.makeNewDataFrameFromSelection(self.arrayData, self.isElectron)

    def getNumberEventsPassingSelection(self, particle = None):
        if particle == None:
            return len(self.getDataFrameDetector(0)["spillNumber"])
        else:
            if self.getDataFrameAllDetectors(particle)!=None:
                return len(self.getDataFrameAllDetectors(particle)[0]["spillNumber"])
            else:
                return None
            
    def plotTOFbounds(self):
        '''Simple function to plot the TOF bounds used to select the proton and deuterium like events'''
        p_range = np.arange(200, 1200, 10)
        upper = []
        lower = []
        expected_proton = []
        upperD = []
        lowerD = []
        expected_deuterium = []
        for p in p_range:
            TOFprotonLower,TOFprotonUpper =  self.getProtonTOFSelectionBounds(p)
            upper.append(TOFprotonUpper)
            lower.append(TOFprotonLower)
            expected_proton.append(self.momentumToTOF(p, "proton"))
            
            TOFdeuteriumLower,TOFdeuteriumUpper = self.getDeuteriumTOFSelectionBounds(p)
            upperD.append(TOFdeuteriumUpper)
            lowerD.append(TOFdeuteriumLower)
            expected_deuterium.append(self.momentumToTOF(p, "deuterium"))
        fig, ax = plt.subplots(1, 1, figsize = (16, 9))
        ax.fill_between(p_range, lower, upper, color = 'blue', alpha = 0.3, label = 'Proton TOF selection Bounds')
        ax.plot(p_range, expected_proton, 'b--', label = 'Nominal proton TOF')
        # ax.fill_between(p_range, lowerD, upperD, color = 'red', alpha = 0.3, label = 'Deuterium TOF selection Bounds')
        # ax.plot(p_range, expected_deuterium, 'r--', label = 'Nominal deuterium TOF')
        ax.set_xlabel("Nominal beam momentum (MeV/c)",fontsize = 18)
        ax.set_ylabel("Time of flight (ns)",fontsize = 18)
        ax.tick_params(axis='both', which='major', labelsize=15)
        
        ax.grid()
        ax.legend(fontsize = 18)
        if self.isLowMomentum:
            ax.fill_between(p_range, lowerD, upperD, color = 'red', alpha = 0.3, label = 'Deuterium TOF selection Bounds')
            ax.plot(p_range, expected_deuterium, 'r--', label = 'Nominal deuterium TOF')
            ax.set_title("Time of flight selection for\n protons and deuterium - Low momentum configuration",fontsize = 20, weight = 'bold')
            fig.savefig("../new_png_results/SelectionTOFBoundsLM.png")
            fig.savefig("../new_pdf_results/SelectionTOFBoundsLM.pdf")
        else:
            ax.set_title("Time of flight selection for\n protons and deuterium - Tagged photon configuration",fontsize = 20, weight = 'bold')
            fig.savefig("../new_png_results/SelectionTOFBoundsTG.png")
            fig.savefig("../new_pdf_results/SelectionTOFBoundsTG.pdf")




    def getProbaParticleInBunch(self):
        """TODO: complete with the calculation, as described in the beam flux paper by Dean of the calculation of the probability of having a particle in a given bunch, based on the delay between particles arrival times"""
        return self.probaBunch
    
    def plotBeamSpot(self, particle=None):
        """TODO: write a function that plots for each particle type (or all particles when None) where those particles hit the TS based on the weighted sum of their windowIntPE signal in the TS, see Jiri's study"""


            

    def outputResults(self):
        """Output all of the relevant information into a csv file"""

        entriesNames = ["runNumber",
                        "runMomentum", 
                        "runRefractiveIndex",
                        "totalNumberSpills",  
                        "probabilityToHaveParticleInBunch", "totalNumberOfEvents", 
                        "numberOfEventsPassingBaseSelections", 
                        "numberOfEventsPassingElectronSelection", 
                        "numberOfEventsPassingMuonSelection", 
                        "numberOfEventsPassingPionSelection", 
                        "numberOfEventsPassingProtonSelection",
                        "numberOfEventsPassingDeuteriumSelection",
                        "meanTOFelectron",
                        "meanTOFmuon",
                        "meanTOFpion",
                        "meanTOFproton",
                        "meanTOFdeuterium",
                        "stdTOFelectron",
                        "stdTOFmuon",
                        "stdTOFpion",
                        "stdTOFproton",
                        "stdTOFdeuterium",
                        "meanMomentumMuon",
                        "meanMomentumPion",
                        "meanMomentumProton",
                        "meanMomentumDeuterium",
                        "meanMomentumStatErrMuon",
                        "meanMomentumStatErrPion",
                        "meanMomentumStatErrProton",
                        "meanMomentumStatErrDeuterium",
                        "numberOfTOFfittedElectron",
                        "numberOfTOFfittedMuon",
                        "numberOfTOFfittedPion",
                        "numberOfTOFfittedProton",
                        "numberOfTOFfittedDeuterium",
                        ]

        entries = {"runNumber": self.runNumber,
                   "runMomentum": self.runMomentum,
                   "runRefractiveIndex": self.runRefractiveIndex,
                   "totalNumberSpills": max(self.getDataFrameDetector(0)["spillNumber"]),
                   "probabilityToHaveParticleInBunch": self.getProbaParticleInBunch(),
                   "totalNumberOfEvents": self.totalNumberOfEvents,
                   #nCoincidence, TS charge cut, etc... for each perticle
                   "numberOfEventsPassingBaseSelections": self.getNumberEventsPassingSelection(),
                   "numberOfEventsPassingElectronSelection": self.getNumberEventsPassingSelection("electron"),
                   "numberOfEventsPassingMuonSelection": self.getNumberEventsPassingSelection("muon"),
                   "numberOfEventsPassingPionSelection": self.getNumberEventsPassingSelection("pion"),
                   "numberOfEventsPassingProtonSelection": self.getNumberEventsPassingSelection("proton"),
                   "numberOfEventsPassingDeuteriumSelection":self.getNumberEventsPassingSelection("deuterium"),
                   #TOF estimates
                   "meanTOFelectron": self.dictTOFMean["electron"],"meanTOFmuon": self.dictTOFMean["muon"],
                   "meanTOFpion": self.dictTOFMean["pion"],
                   "meanTOFproton": self.dictTOFMean["proton"],
                   "meanTOFdeuterium": self.dictTOFMean["deuterium"],
                   "stdTOFelectron": self.dictTOFStd["electron"],
                   "stdTOFmuon": self.dictTOFStd["muon"],
                   "stdTOFpion": self.dictTOFStd["pion"],
                   "stdTOFproton": self.dictTOFStd["proton"],
                   "stdTOFdeuterium": self.dictTOFStd["deuterium"],
                   #momentum estimates
                   "meanMomentumMuon": self.dictMomentumMean["muon"],
                   "meanMomentumPion": self.dictMomentumMean["pion"],
                   "meanMomentumProton": self.dictMomentumMean["proton"],
                   "meanMomentumDeuterium": self.dictMomentumMean["deuterium"],
                   "meanMomentumStatErrMuon": self.dictMomentumStatError["muon"],
                   "meanMomentumStatErrPion": self.dictMomentumStatError["pion"],
                   "meanMomentumStatErrProton": self.dictMomentumStatError["proton"],
                   "meanMomentumStatErrDeuterium": self.dictMomentumStatError["deuterium"],
                   "numberOfTOFfittedElectron": self.dictTOFfittedNparticles["electron"],
                   "numberOfTOFfittedMuon": self.dictTOFfittedNparticles["muon"],
                   "numberOfTOFfittedPion": self.dictTOFfittedNparticles["pion"],
                   "numberOfTOFfittedProton": self.dictTOFfittedNparticles["proton"],
                   "numberOfTOFfittedDeuterium": self.dictTOFfittedNparticles["deuterium"],

                   
        }


        
        if not(os.path.exists(self.outputFileName)):
            with open(self.outputFileName, "w") as f:
                writer = csv.DictWriter(f, fieldnames = entriesNames)
                writer.writeheader()
                writer.writerow(entries)
        else:
            with open(self.outputFileName, "a") as f:
                writer = csv.DictWriter(f, fieldnames = entriesNames)
                writer.writerow(entries)
            
            
        
        







        



