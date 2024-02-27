import uproot as ur

import matplotlib.pyplot as plt

import numpy as np

import pandas as pd
import sys

import matplotlib.colors as colors

from scipy.optimize import curve_fit
from collections import OrderedDict

from scipy.integrate import quad

from scipy.optimize import curve_fit

############################### READ ME ###############################################
#started by acraplet
#This is the selection of particles and ID
#Later on can add the momentum measurement using the TOF inoformation for muons and pions - DONE
#######################################################################################

#the selection variables are global
global ACTlinearA, ACTlinearB, piMuBorderACT, ACTlower, thereIsProtons, protonsTOFCut, horizontal_el, LGupper, protonsTOFMax

global nProtons, nPions, nMuons, nElectrons, nParticles, nDeuterium



#the mass of the particles in MeV/c
ms = OrderedDict()
ms['Electrons'] = 0.511
ms['Muons'] = 105.658
ms['Pions'] = 139.6
ms['Protons'] = 938.3
ms['Deuterium'] = 1876.

conv = 1.e9 #convert back to ns
#TB2022 l = 2.90 # m
#TB2023
l = 3.49 # m
c = 299792458 # m/s

def gaussian(x, amplitude, mean, std_dev):
    return amplitude * np.exp(-0.5 * ((x - mean) / std_dev) ** 2)

def noHoleHit(file, conditionInitial):
    #TODO: this has to be modified to add Deesha's cut on hole counter amplitude (30mV)
    df = file['Hole0'].arrays(library="pd")
    condition = (df['nPeaks']==0) & conditionInitial
    df = file['Hole1'].arrays(library="pd")
    condition = (df['nPeaks']==0) & condition
    return condition

def nPeakInToF(file, n):
    #TODO This will have to be changed for Dean's 2 particle veto selection based on integration of the first 200ns of the waveform.
    if n == False:
        df = file['TOF00'].arrays(library="pd")
        #in the case we want to see everything
        return (df['nPeaks']!= -1)

    else:
        df = file['TOF00'].arrays(library="pd")
        nPeaksInToF = (df['nPeaks']== n)
        df = file['TOF01'].arrays(library="pd")
        nPeaksInToF = (df['nPeaks']==n) & nPeaksInToF
        df = file['TOF02'].arrays(library="pd")
        nPeaksInToF = (df['nPeaks']==n) & nPeaksInToF
        df = file['TOF03'].arrays(library="pd")
        nPeaksInToF = (df['nPeaks']==n) & nPeaksInToF
        df = file['TOF10'].arrays(library="pd")
        nPeaksInToF = (df['nPeaks']==n) & nPeaksInToF
        df = file['TOF11'].arrays(library="pd")
        nPeaksInToF = (df['nPeaks']==n) & nPeaksInToF
        df = file['TOF12'].arrays(library="pd")
        nPeaksInToF = (df['nPeaks']==n) & nPeaksInToF
        df = file['TOF13'].arrays(library="pd")
        nPeaksInToF = (df['nPeaks']==n) & nPeaksInToF

        return nPeaksInToF

def printSelection(ACTlinearA, ACTlinearB, piMuBorderACT, ACTlower, thereIsProtons, protonsTOFCut, horizontal_el, LGupper):
     print("\n#############  Summary of the selection rule ############# \n",
          "thereIsProtons = ", thereIsProtons, "\n ",
          "protonsTOFCut = ", protonsTOFCut, "ns (if above is True) \n ",
          "protonsTOFMax = ", protonsTOFMax, " ns\n",
          "ACTlinearA = ", ACTlinearA, "\n ",
          "ACTlinearB = ", ACTlinearB, "\n ",
          "piMuBorderACT = ", piMuBorderACT, "\n ",
          "ACTlower = ", ACTlower, "\n ",
          "horizontal_el = ", horizontal_el, "\n",
          "LGupper = ", LGupper, "\n",
          "############# End of the selection summary #############\n")
     return 0

def plotSelection(df_all, ACTlinearA, ACTlinearB, piMuBorderACT, ACTlower, thereIsProtons, protonsTOFCut, horizontal_el, max_ACT23 = 200, max_ACT1 = 100, plot_events = True, ax = None):

    #make the avearge of the PMTs
    ACT23 = (df_all[4][0]  + df_all[5][0] + df_all[6][0]  + df_all[7][0])
    ACT1 = (df_all[2][0]  + df_all[3][0])

    xmin, xmax = min(ACT1) * 0.8,  min(max_ACT1, max(ACT1)* 0.8)
    ymin, ymax = min(ACT23) * 0.8, min(max_ACT23, max(ACT23)* 0.8)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    if plot_events:
        plt.plot(np.array(ACT1), np.array(ACT23), "x", color = "lightgrey")
        plt.plot([xmin, xmax], [piMuBorderACT, piMuBorderACT], "g--", label ="Pi-Mu separation %.3f PE" %(piMuBorderACT))
        plt.plot([xmin, xmax], [ACTlower, ACTlower], "r-", label ="Limit noise %.3f PE"%(ACTlower))
        plt.plot([xmin, (horizontal_el - ACTlinearB)/ACTlinearA], [xmin * ACTlinearA + ACTlinearB, ((horizontal_el - ACTlinearB)/ACTlinearA) * ACTlinearA + ACTlinearB], "k--", label ="Diagonal electron veto A = %.3f, B = %.3f" %(ACTlinearA, ACTlinearB))

        plt.plot([(horizontal_el - ACTlinearB)/ACTlinearA, xmax], [horizontal_el, horizontal_el], "k--", label ="Horizontal electron veto %.3f PE"%(horizontal_el))
        plt.legend()

        plt.xlabel("ACT1 (PE)")
        plt.ylabel("ACT23 (PE)")
        # plt.close()
        plt.show()

    else: # just plot the selection lines
        ax.plot([xmin, xmax], [piMuBorderACT, piMuBorderACT], "g--")
        ax.plot([xmin, xmax], [ACTlower, ACTlower], "r-")
        ax.plot([xmin, (horizontal_el - ACTlinearB)/ACTlinearA], [xmin * ACTlinearA + ACTlinearB, ((horizontal_el - ACTlinearB)/ACTlinearA) * ACTlinearA + ACTlinearB], "k--")
        ax.plot([(horizontal_el - ACTlinearB)/ACTlinearA, xmax], [horizontal_el, horizontal_el], "k--")

def plotTOFSelection(tof, ACT23, thereIsProtons, protonsTOFCut, max_ACT23=200):
    xmin, xmax = min(tof) * 0.8, max(tof)* 0.8
    ymin, ymax = min(ACT23) * 0.8, min(max_ACT23, max(ACT23)* 0.8)
    plt.plot(np.array(tof), np.array(ACT23), "x", color = "lightgrey")
    if thereIsProtons:
        plt.plot([protonsTOFCut, protonsTOFCut], [ymin, ymax], "r--", label ="Proton TOF cut %.3f ns" %(protonsTOFCut))
        plt.ylabel("ACT23 (PE)")
        plt.xlabel("TOF (ns)")
        plt.legend()
        plt.grid()
        # plt.close()
        plt.show()

def distToLinearCut(x1, y1, a, b):
    #simple trigonometry
    yA = a*x1 + b
    xA = (y1-b)/a
    theta = np.arctan((yA-y1)/(x1-xA))
    return np.sin(theta) * (x1-xA)


def makeSelection(df_all, df_all_intPE, df_all_times, tof):
    '''The selection of the particles is as follows:
    electrons: evereything that is above a diagonal line or the horizontal_el line
    this later condition is helping with getting pions and muons that might have scintillated in ACT1
    muons: what is below the diagonal line and above the horizontal piMuBorderACT line
    pions: what is below both the diagonal line and the piMuBorderACT line & above ACTlower (noise)
    protons: if there is protons then anything that is arriving later than the TOF time (ns) given
    The position of the horizontal piMuBorderACT line depends on how much scintillation there is
    This function outputs 4 arrays of bools telling which particle is which, need to add a safety
    check that no particle is simultaneously two different particles
    '''

    #Initialise the number of particles
    global ACTlinearA, ACTlinearB, piMuBorderACT, ACTlower, thereIsProtons, protonsTOFCut, horizontal_el, LGupper, protonsTOFMax
    global nProtons, nPions, nMuons, nElectrons, nParticles, nDeuterium

    nProtons, nPions, nMuons, nElectrons, nParticles = 0, 0, 0, 0, 0


    #some sanity checks:
    printSelection(ACTlinearA, ACTlinearB, piMuBorderACT, ACTlower, thereIsProtons, protonsTOFCut, horizontal_el, LGupper)


    #for now, not including the scintillation cuts, but later on can do that
    #that is, if any signal is arriving too late, we do not use its charge in the selection
    multiplicator_right = 1.
    multiplicator_left = 1.

    #make the avearge of the PMTs
    ACT23 = (df_all[4][0] * multiplicator_left  + df_all[5][0] * multiplicator_right + df_all[6][0] * multiplicator_left + df_all[7][0] * multiplicator_right)

    ACT1 = (df_all[2][0] * multiplicator_left + df_all[3][0] * multiplicator_right)

    ACT0 = (df_all[1][0] + df_all[0][0])

    #for calculating the number of particles, use the LG
    LG = df_all_intPE[18][0]

    #also calculate the TOF value, for proton selection
    TOF_all = pd.Series(data = tof, name = 'TOF')


    ############## check selection rules visually ################
    #maybe here we should plot the selection, to see if it looks sensible
    #and then have some user input to change the values and test again! that would be great
    plotSelection(df_all, ACTlinearA, ACTlinearB, piMuBorderACT, ACTlower, thereIsProtons, protonsTOFCut, horizontal_el, max_ACT23 = 300, max_ACT1 = 60)

    plotTOFSelection(tof, ACT23, thereIsProtons, protonsTOFCut)

    keep_selection = input("Are you happy with this selection? y/N ")
    while keep_selection != "y":
        print("Changing the selection, please input the desired value for\n")
        thereIsProtons = bool(int(input("\nthereIsProtons: True/False, previous value: %i "% int(thereIsProtons))))
        protonsTOFCut = float(input("\nprotonsTOFCut (ns), previous value: %.3f "% float(protonsTOFCut)))
        protonsTOFMax = float(input("\nupper protonsTOFCut (ns), previous value: %.3f "% float(protonsTOFMax)))
        ACTlinearA = float(input("\nACTlinearA: previous value:%.3f "% float(ACTlinearA)))
        ACTlinearB = float(input("\nACTlinearB: previous value:%.3f "% float(ACTlinearB)))
        piMuBorderACT = float(input("\npiMuBorderACT: previous value:%.3f "% float(piMuBorderACT)))
        ACTlower = float(input("\nACTlower: previous value:%.3f "% float(ACTlower)))
        horizontal_el = float(input("\nhorizontal_el: previous value:%.3f "% float(horizontal_el)))
        LGupper =  float(input("\nLGupper: previous value:%.3f "% float(LGupper)))

        print("\nThank you for updating your selection, this is your new selection\n")
        printSelection(ACTlinearA, ACTlinearB, piMuBorderACT, ACTlower, thereIsProtons, protonsTOFCut, horizontal_el, LGupper)
        plotSelection(df_all, ACTlinearA, ACTlinearB, piMuBorderACT, ACTlower, thereIsProtons, protonsTOFCut, horizontal_el)

        keep_selection = input("Are you happy with this selection? y/N ")

    nParticles = len(df_all[0][0])

    ##First identify the protons if there are any
    if thereIsProtons:
        protonOrD_selection = np.where(TOF_all>protonsTOFCut, True, False)

        print("\n\nWe have a total of %i protons for %i total number of events, i.e. %.4f percent"%(nProtons, nParticles, nProtons/nParticles * 100))
        proton_selection = np.where(TOF_all<protonsTOFMax, protonOrD_selection, False)
        deuterium_selection = np.where(TOF_all>protonsTOFMax, protonOrD_selection, False)
        nProtons = np.sum(proton_selection)
        nDeuterium = np.sum(deuterium_selection)

    else:
        proton_selection = np.where(df_all[0][0] != -9999, False, False)
        print("\n\nWe have a total of %i protons for %i total number of events, i.e. %.4f percent"%(nProtons, nParticles, nProtons/nParticles * 100))#
        protonOrD_selection = proton_selection
        nDeuterium = 0
        nProtons = 0

    #deuterium has already been identified even if it is not a proton
    hasBeenIDed = np.where(protonOrD_selection == True, True, False)




    ###################################################################
    ##Second, flag the electrons, bitwise operation....
    ###################################################################
    electron_selection = np.where((ACT23 > ACT1 * ACTlinearA + ACTlinearB), True, False)
    #do not put the electrons that are inbetween the diagonal and horozontal
    #if you are not above you might not be an electron

    horizontalCondition = np.where(ACT23>horizontal_el, True, False)
    electron_selection = np.where(electron_selection == True, horizontalCondition, False)

    #TODO:this has to be changed for the actual WCTE run!!
    electron_selection = np.where(LG>=LGupper, True, electron_selection)

    #Ignore the ones that are below the ACT23
    electron_selection = np.where(ACT23<=ACTlower, False, electron_selection)

    #cannot be another particle aleady
    electron_selection = np.where(hasBeenIDed == False, electron_selection, False)

    nElectrons = np.sum(electron_selection)
    print("We have a total of %i electrons for %i total number of events, i.e. %.4f percent"%(nElectrons, nParticles, nElectrons/nParticles * 100))
    #update the counter of particles identified
    hasBeenIDed = np.where(electron_selection == True, True, hasBeenIDed)

    ############################################################################
    ##Third, ID the muons, technically, only need to not have been identified yet
    #and need to be above the pion line
    ############################################################################
    muon_selection = np.where((ACT23 > piMuBorderACT), True, False)
    muon_selection = np.where((hasBeenIDed == False), muon_selection, False)
    nMuons = np.sum(muon_selection)
    print("We have a total of %i muons for %i total number of events, i.e. %.4f percent"%(nMuons, nParticles, nMuons/nParticles * 100))
    hasBeenIDed = np.where(muon_selection == True, True, hasBeenIDed)

    ################################################################################
    #Fourth, ID the pions, technically, should been evereything that hasn't been yet
    #IDed and is above the limit in ACT23 (noise)
    ################################################################################
    #need to go in steps otherwise overlaps
    pion_selection = np.where((ACT23 >= ACTlower), True, False)
    pion_selection = np.where(ACT23 < piMuBorderACT, pion_selection, False)

    pion_selection = np.where(hasBeenIDed == False, pion_selection, False)
    nPions = np.sum(pion_selection)
    print("We have a total of %i pions for %i total number of events, i.e. %.4f percent"%(nPions, nParticles, nPions/nParticles * 100))
    hasBeenIDed = np.where(pion_selection == True, True, hasBeenIDed)

    #################################################################################
    #Fifth, sanity check: did we miss anything?
    #################################################################################
    nParticleIDed = nProtons + nElectrons + nPions + nMuons #np.sum(hasBeenIDed)
    print("We have a total of %i IDed particles for %i total number of events, i.e. %.4f percent"%(nParticleIDed, nParticles, nParticleIDed/nParticles * 100))
    #Is what we missed indeed what we wanted to reject?
    noise_selection = np.where((ACT23 < ACTlower), True, False)
    noise_selection =  np.where((hasBeenIDed == False), noise_selection, False)
    nNoiseParticle = np.sum(noise_selection)
    print("We have rejected a total of %i particles for %i total number of events, i.e. %.4f percent"%(nNoiseParticle, nParticles, nNoiseParticle/nParticles * 100))

    print("In total, the selection has looked at %i particles out of %i events, i.e. %.4f percent \n\n" % (nParticleIDed + nNoiseParticle, nParticles, (nParticleIDed + nNoiseParticle)/nParticles * 100))

    return proton_selection, electron_selection, muon_selection, pion_selection

def getHist(array_df):
    Mean_ACT23 = (array_df[4][0] + array_df[5][0] +array_df[6][0] + array_df[7][0])
    Mean_ACT1 = (array_df[3][0] + array_df[2][0] )
    tof = array_df[0]['TOF']
    leadGlass = array_df[18][0]
    return Mean_ACT23, Mean_ACT1, tof, leadGlass

def getHist_ACT0(array_df):
    Mean_ACT0 = (array_df[0][0] + array_df[1][0])
    return Mean_ACT0


def getTimingHist(array_df):
    return array_df[2][0], array_df[3][0], array_df[4][0], array_df[5][0], array_df[6][0], array_df[7][0]

def getRefTiming(array_df):
    return (array_df[12][0] + array_df[13][0] + array_df[14][0] + array_df[15][0])/4


def TofToMomentum(tof, particle_label):
    #the tof needs to be the absolute flying time
    m = ms[particle_label]
    p = m/np.sqrt(pow((tof) * c / (conv * l), 2) - 1)
    return p

def getTof(m, momentum):
    return l/c * np.sqrt(1. + m**2/momentum**2) * conv



def singlePE(argv):
    #the selection variables are global
    global ACTlinearA, ACTlinearB, piMuBorderACT, ACTlower, thereIsProtons, protonsTOFCut, horizontal_el, LGupper, protonsTOFMax

    global nProtons, nPions, nMuons, nElectrons, nParticles, nDeuterium, signalTimeBranch

    signalTimeBranch = "SignalTimeCorrected" #can change for the old timings (SignalTime)

    #get some preliminary guesses for the parameters of the selection
    headers = ["Run", "momentum", "refIndex", "nSpill",  "probaBunch", "nParticles", "nElectrons", "nMuons", "nPions", "nProtons", "nDeuterium", "fractionPass1ParticleVeto", "fractionPassNanVeto", "ACTlinearA", "ACTlinearB", "piMuBorderACT", "ACTlower", "thereIsProtons", "protonsTOFCut", "horizontal_el", "LGupper", "bery"]

    df_ref = pd.DataFrame()

    #the root file
    root_filenames = argv[1]

    #This is for labelling the runs correctly
    run = 0
    momentum = 0
    refIndex = 0
    probaBunch = 1

    #in case we want to create an identical root file with the particle type included
    #is very slow but can be useful for comparing selections for example.
    saveRootFile = False #True

    #Read in the user input run characterisation
    if len(argv) >= 5:
        run = int(argv[2])
        #this is the reference with the cuts that I have derived -> good to use as starting point
        df_ref = pd.read_csv('referenceNumberParticles.txt', names = headers, sep = " ", index_col=None)
        momentum = float(argv[3])
        #same momentum should be behaving the same
        df_ref = df_ref[df_ref["momentum"] == momentum]
        refIndex = float(argv[4])

    isBeryllium = 1 #one is berylium, default
    targetType = ["Aluminium", "Beryllium"]

    if len(argv) >= 6:
        isBeryllium = int(argv[5])
        print("The target chosen was %s"%targetType[isBeryllium])

    if len(argv) == 7:
        probaBunch = float(argv[6])

    else:
        print("\n\n\n##########################################\nUsage: python python/particleID.py rootFilePath run momentum n isBeryllium probaParticleInBunch\n##########################################\n\n\n")



    #read in the file
    file = ur.open(root_filenames)

    #Require only one particle in all of the runs
    nPeaksInToF = nPeakInToF(file, 1)

    #Might be useful to calculate the fractions of events that pass this 1 particle cut
    fractionPass1ParticleVeto = np.sum(nPeaksInToF)/len(nPeaksInToF)
    print("\n(nParticles) There is excatly one particle in all of the TOF detectors for %i events out of %i that is %.3f percent \n" % (np.sum(nPeaksInToF), len(nPeaksInToF), fractionPass1ParticleVeto * 100))


    #for now we are only working with a cut on the ACT23 vs ACT1 plane
    df_all, df_e, df_p, df_mu, df_pi, df_total = [], [], [], [], [], []
    df_T_e, df_T_p, df_T_mu, df_T_pi, df_T_total = [], [], [], [], []
    df_complete_all = []

    #If we want that there is no hits in the hole detectors
    # nPeaksInToF = noHoleHit(file, nPeaksInToF)

    #append the data to the initial dataframe
    for detector in range(len(file.keys()[:-1])):
        key = file.keys()[detector]
        df_complete = file[key].arrays(library="pd")
        df = df_complete[nPeaksInToF]
        list_index = df.index
        all_index = df_complete.index
        df = df.reset_index()
        df_all.append(df)
        df_complete["onlyOneParticle"] = np.array(nPeaksInToF).astype(int)
        df_complete_all.append(df_complete)

    key = file.keys()[-1]
    df_enventInfo =  file[key].arrays(library="pd")


    print((nPeaksInToF == 1), df_complete["onlyOneParticle"], list_index)

    #need to store the peak integrated charge for some detectors (LG mainly)
    df_all_intPE = df_all.copy()
    df_all_times = df_all.copy()


    windPEisNan = False #remove NaN issues
    #need to split the vector of collected charge so we access each element as a column
    for pmt in range(len(file.keys()[:-1])):
        intCharge = pd.DataFrame(df_all_intPE[pmt]['IntPE'].values.tolist())
        windowCharge = pd.DataFrame(df_all[pmt]['WindowIntPE'].values.tolist())
        hitTimes = pd.DataFrame(df_all_times[pmt]['%s'%signalTimeBranch].values.tolist())

        df_all[pmt] = pd.concat([df_all[pmt],windowCharge], axis=1)
        df_all_intPE[pmt] = pd.concat([df_all_intPE[pmt],intCharge], axis=1)
        df_all_times[pmt] = pd.concat([df_all_times[pmt],hitTimes], axis=1)

        #neeed at least one entry, if we have any PMTs in which we did not collect
        #any data in ACT23, technically, shouldn't need it
        if pmt <= 7:
            windPEisNan = windPEisNan | df_all[pmt].isna()[0]

    #check how much we loose after removing the nans
    print("(Nans) At least one of the ACTs has NaN in windPE for %i events out of %i that is we keep %.3f percent of events \n" % (np.sum(windPEisNan), len(windPEisNan), (len(windPEisNan)-np.sum(windPEisNan))/len(windPEisNan) * 100))

    fractionPassNanVeto = (len(windPEisNan)-np.sum(windPEisNan))/len(windPEisNan)

    if np.sum(windPEisNan) !=0:
        #Remove the NaNs if necessary
        windPEisNotNan = np.where(windPEisNan == True, False, True)


        for pmt in range(len(file.keys()[:-1])):
            df_all[pmt] = df_all[pmt][windPEisNotNan]
            df_all_intPE[pmt] = df_all_intPE[pmt][windPEisNotNan]
            df_all_times[pmt] = df_all_times[pmt][windPEisNotNan]
            #need to reset index, useful for later
            df_all[pmt] = df_all[pmt].reset_index()
            df_all_intPE[pmt] = df_all_intPE[pmt].reset_index()
            df_all_times[pmt] = df_all_times[pmt].reset_index()


    # #Check the total number of spills in this run
    TOTAL_spill = max(df_all[0]['spillNumber'])
    print("Number of spills:", TOTAL_spill)

    #calculate the particle travel time from the mean flare time of the TOF detectors
    #later on we will append that to the dataframes for each detector, sueful to have
    mean_first_hit_TOF1 = (df_all_times[12][0] + df_all_times[13][0] + df_all_times[14][0] + df_all_times[15][0])/4
    mean_first_hit_TOF0 = (df_all_times[8][0] + df_all_times[9][0] + df_all_times[10][0] + df_all_times[11][0])/4

    tof = np.array(np.array(mean_first_hit_TOF1)-np.array(mean_first_hit_TOF0))
    TOF = pd.Series(data = tof, name = 'TOF')


    #Now extract the particle ID
    #2d cut in ACT1 ACT23 and TOF for the protons
    #THESE ARE GLOBAL VARIABLES, NEED TO BE UPDATED AS WE GO
    ACTlinearA = -6#-3 #the diagonal parameters
    ACTlinearB = 73#73#73#90
    piMuBorderACT = 30#17 #30 #25 #38 #pe
    ACTlower = 1
    protonsTOFCut = 15 # 30 #20 #ns
    protonsTOFMax = 35#32 #20 #ns
    thereIsProtons = (momentum>= 0) # False #True
    horizontal_el = 35#35#40 #70
    LGupper = 3#2.5#2.5 #4

    if len(df_ref != 0):
        #we have reference values, caluculated previously, use those instead
        ACTlinearA, ACTlinearB, piMuBorderACT = float(df_ref["ACTlinearA"]), float(df_ref["ACTlinearB"]), float(df_ref["piMuBorderACT"])
        ACTlower, protonsTOFCut,  LGupper = float(df_ref["ACTlower"]), float(df_ref["protonsTOFCut"]), float(df_ref["LGupper"])
        horizontal_el = float(df_ref["horizontal_el"])
        protonsTOFMax = protonsTOFCut + 8

    #Identify the particles, the limits are global variables, careful
    proton_selection, electron_selection, muon_selection, pion_selection = makeSelection(df_all, df_all_intPE, df_all_times, tof)

    #Create a new dataframe with the particle type, so we can save it with the dataframe later
    df_particleType = pd.DataFrame()

    df_particleType["isElectron"] = np.array(electron_selection).astype(int)
    df_particleType["isMuon"] = np.array(muon_selection).astype(int)
    df_particleType["isPion"] = np.array(pion_selection).astype(int)
    df_particleType["isProton"] = np.array(proton_selection).astype(int)

    #here, set the index of the hits where we had more than 1 particle in the waveform
    df_particleType = df_particleType.set_index(list_index)
    #and here, extend the dataframe, filling zeros, everywhere that doesn't have only one particle in the waveform
    df_particleType = df_particleType.reindex(index = all_index, fill_value = 0)

    #for Dean, export the particle type
    df_particleType.to_csv("Run%i_eventID_v2.csv"%run, sep = ',', index = True)



    ############################ save a copy of the root file with the particle identification
    if saveRootFile:
        # #here append the particle type to the dataframe
        for pmt in range(len(file.keys()[:-1])):
            df_complete_all[pmt] = pd.concat([df_complete_all[pmt], df_particleType], axis=1)

        #Now, export the root file with the particle identification, which can be useful
        df_enventInfo['probaBunch'] = np.array([probaBunch] * len(all_index))
        df_enventInfo['berylliumTarget'] =  np.array([isBeryllium] * len(all_index))
        df_enventInfo['beamMomentum'] =  np.array([momentum] * len(all_index))

        print("\nSaving a root file with the particle ID, please be patient, this is a slow process")
        with ur.recreate('particleIDed_%s'%(argv[1])) as file1:
            for pmt in range(len(file.keys()[:-1])):
                # tree = ur.newtree'%s'%(file.keys()[pmt]), {col: df_complete_all[pmt][col].dtype for col in df_complete_all[pmt].columns})
                file1["%s"%(file.keys()[pmt][:-2])] = df_complete_all[pmt]
                file1["%s"%(file.keys()[pmt][:-2])].show()
            file1["%s"%(file.keys()[-1][:-2])] = df_enventInfo  #need to add this at the end, otherwise it messes it up





    #for bookkepping, keep all of the particles
    total_selection = np.where(df_all[0][0] != -9999, True, True)

    #all the dataframes that we will append to

    list_data_frames = [df_total, df_e, df_mu, df_pi, df_p]
    list_selection_bools = [total_selection, electron_selection, muon_selection, pion_selection, proton_selection]
    list_particle_type = ["All events", "Electrons", "Muons", "Pions", "Protons"]
    list_particle_numbers = [nParticles, nElectrons, nMuons, nPions, nProtons]

    #we might want to check the timing of the hits -> is it scintillation?
    list_time_frames = [df_T_total, df_T_e, df_T_mu, df_T_pi, df_T_p]

    #Now fill the dataframe correctly with the time of flight
    for detector in range(len(file.keys()[:-1])):

        #add the TOF
        df_all_intPE[detector] = pd.concat([df_all_intPE[detector],TOF], axis=1)
        df_all[detector] = pd.concat([df_all[detector],TOF], axis=1)

        #add the hit timings
        for particleType in range(len(list_data_frames)):
            if (list_particle_type[particleType] == "Protons" and (not thereIsProtons)):
                    break
            list_time_frames[particleType].append(df_all_times[detector][list_selection_bools[particleType]])

        #add the charge
        if detector <= 7:
            for particleType in range(len(list_data_frames)):
                if (list_particle_type[particleType] == "Protons" and (not thereIsProtons)):
                    break
                list_data_frames[particleType].append(df_all[detector][list_selection_bools[particleType]])
        else:
            for particleType in range(len(list_data_frames)):
                if (list_particle_type[particleType] == "Protons" and (not thereIsProtons)):
                    break
                list_data_frames[particleType].append(df_all_intPE[detector][list_selection_bools[particleType]])

    #Have now a clean way to plot things? technically we already have the numbers
    #of particles, they are clean, get the histograms to check
    list_ACT23, list_ACT1, list_TOF, list_LG = [], [], [], []
    timing_ACT1L, timing_ACT1R, timing_ACT2L, timing_ACT2R, timing_ACT3L, timing_ACT3R, timing_reference = [], [], [], [], [], [], []
    list_ACT0 = []

    Mean_ACT23_mu, Std_ACT23_mu, Mean_ACT23_pi, Std_ACT23_pi, Mean_ACT23_p, Std_ACT23_p = 0,0,0,0,0,0



    #read in the data
    for particleType in range(len(list_data_frames)):

        if (list_particle_type[particleType] == "Protons" and (not thereIsProtons)):
                    break
        temp_ACT23, temp_ACT1, temp_TOF, temp_LG = getHist(list_data_frames[particleType])
        temp_ACT0 = getHist_ACT0(list_data_frames[particleType])

        list_ACT23.append(temp_ACT23)
        list_ACT1.append(temp_ACT1)
        list_ACT0.append(temp_ACT0)
        list_TOF.append(temp_TOF)
        list_LG.append(temp_LG)

        if list_particle_type[particleType] == "Pions":
                    Mean_ACT23_pi = temp_ACT23.mean()
                    Std_ACT23_pi = temp_ACT23.std()

        if list_particle_type[particleType] == "Muons":
            Mean_ACT23_mu = temp_ACT23.mean()
            Std_ACT23_mu = temp_ACT23.std()

        if list_particle_type[particleType] == "Protons":
            Mean_ACT23_p = temp_ACT23.mean()
            Std_ACT23_p = temp_ACT23.std()

    #read in the hit times
    for particleType in range(len(list_time_frames)):
        if (list_particle_type[particleType] == "Protons" and (not thereIsProtons)):
                break
        temp_ACT1L, temp_ACT1R, temp_ACT2L, temp_ACT2R, temp_ACT3L, temp_ACT3R = getTimingHist(list_time_frames[particleType])
        #the reference timing is the hit time iof the TOF (mean)
        temp_ref = getRefTiming(list_time_frames[particleType])
        timing_ACT1L.append(temp_ACT1L)
        timing_ACT1R.append(temp_ACT1R)
        timing_ACT2L.append(temp_ACT2L)
        timing_ACT2R.append(temp_ACT2R)
        timing_ACT3L.append(temp_ACT3L)
        timing_ACT3R.append(temp_ACT3R)
        timing_reference.append(temp_ref)

    ################################################################################
    ####################### Fit the TOF to get the momentum ! ######################
    ################################################################################
    fig, ax = plt.subplots(figsize = (16,9))

    selectionText = 'Selection used p: thereIsProton = %i, TOF > %.2fns\ne: ACT23 > %.2f ACT1 + %.2f or ACT23 > %.2fpe\nmu: (not e or p) & ACT23 > %.2fpe & LG < %.2fpe \npi: (not e, p or mu) & ACT23 > %.2f pe & LG < %.2fpe' % (thereIsProtons, protonsTOFCut, ACTlinearA, ACTlinearB, horizontal_el, piMuBorderACT, LGupper, ACTlower, LGupper)

    ax.annotate(selectionText, xy=(0.75, 0.7), xytext=(-20, -15), fontsize=14,
    xycoords='axes fraction', textcoords='offset points',
    bbox=dict(facecolor='white', alpha=0.8),
    horizontalalignment='center', verticalalignment='center')


    TOF_min = 10.5
    TOF_max = max(protonsTOFCut, 15.5)
    if thereIsProtons:
        TOF_max = protonsTOFMax + 5
    widthBin = 0.2#ns how finely sample we are
    nBins = int((TOF_max-TOF_min)/widthBin)

    # stacked = 0
    plt.ylim(0.5,nElectrons)

    p_pred_mu, p_pred_pi, p_pred_p, p_error_mu, p_error_pi, p_error_p = 0,0,0,0,0,0



    for particleType in range(len(list_data_frames)):
        #careful we overlay the selection on the other particles
        if particleType != 0:
            if (list_particle_type[particleType] == "Protons" and (not thereIsProtons)):
                    break

            counts, edges = np.histogram(np.array(list_TOF[particleType]), bins=nBins, range=(TOF_min, TOF_max), density=False)

            bin_centers = (edges[:-1] + edges[1:]) / 2
            plotting_x = np.linspace(min(bin_centers), max(bin_centers), 10*nBins)
            errors = np.sqrt(counts / len(np.array(list_TOF[particleType])))

            meanTOF = np.array(list_TOF[particleType]).mean()
            stdTOF = np.array(list_TOF[particleType]).std()


            if list_particle_type[particleType] == "Electrons":
                ax.errorbar(bin_centers, counts, yerr=errors, fmt='o', color = 'darkgray', label = '%s: TOF = %.2f +/- %.2f (%.3f percent of events w/ 1 particle)'%(list_particle_type[particleType], meanTOF, stdTOF,list_particle_numbers[particleType]/nParticles * 100))

                meanTOF_e = meanTOF
                # print(meanTOF_e, getTof(ms['Electrons'], momentum), getTof(ms['Muons'], momentum) -  meanTOF_e + getTof(ms['Electrons'], momentum), getTof(ms['Pions'], momentum) -  meanTOF_e + getTof(ms['Electrons'], momentum))

                ####### need to fit with a gaussian
                # Initial guess for parameters (amplitude, mean, std_dev)
                initial_guess = [len(list_TOF[particleType]), meanTOF, stdTOF]
                # Fit the Gaussian function to the data points
                params, covariance = curve_fit(gaussian, bin_centers, counts, p0=initial_guess)
                # Get the fitted parameters
                amplitude, mean, std_dev = params
                # Calculate the Gaussian values using the fitted parameters
                y_curve = gaussian(plotting_x, amplitude, mean, std_dev)
                #should plot the cuve, I probably was to put a limit, to ignore the other ones?
                ax.plot(plotting_x, y_curve, '--', color = 'lightgray')





            else:
                print(list_particle_type[particleType], ": predicted momentum is ", TofToMomentum(meanTOF - meanTOF_e + getTof(ms['Electrons'], momentum), list_particle_type[particleType]))

                ####### need to fit with a gaussian
                # Initial guess for parameters (amplitude, mean, std_dev)
                initial_guess = [len(list_TOF[particleType]), meanTOF, stdTOF]
                params, covariance = curve_fit(gaussian, bin_centers, counts, p0=initial_guess)
                amplitude, mean, std_dev = params
                # Extract parameter errors
                param_errors = np.sqrt(np.diag(covariance))

                y_curve = gaussian(plotting_x, amplitude, mean, std_dev)
                ax.plot(plotting_x, y_curve, '--', color = 'lightgray')

                print(list_particle_type[particleType], ": predicted with fit momentum is ", TofToMomentum(mean - meanTOF_e + getTof(ms['Electrons'], momentum), list_particle_type[particleType]))

                print(list_particle_type[particleType], ": predicted with fit momentum +1sigma  ", TofToMomentum(mean - meanTOF_e + param_errors[1] + getTof(ms['Electrons'], momentum), list_particle_type[particleType]))

                print(list_particle_type[particleType], ": predicted with fit momentum -1sigma  ", TofToMomentum(mean - meanTOF_e - param_errors[1] + getTof(ms['Electrons'], momentum), list_particle_type[particleType]))

                p_pred = TofToMomentum(mean - meanTOF_e + getTof(ms['Electrons'], momentum), list_particle_type[particleType])

                p_error = p_pred - TofToMomentum(mean - meanTOF_e + param_errors[1] + getTof(ms['Electrons'], momentum), list_particle_type[particleType])

                p_error = (abs(p_error) + abs(p_pred - TofToMomentum(mean - meanTOF_e + param_errors[1] + getTof(ms['Electrons'], momentum), list_particle_type[particleType])))/2

                ax.errorbar(bin_centers, counts, yerr=errors, fmt='o', label = '%s: TOF = %.2f +/- %.2f (%.3f percent of events w/ 1 particle) p_pred = %.2f +/- %.2f MeV/c'%(list_particle_type[particleType], meanTOF, stdTOF,list_particle_numbers[particleType]/nParticles * 100, p_pred, p_error))

                if list_particle_type[particleType] == "Pions":
                    p_pred_pi = p_pred
                    p_error_pi = p_error

                if list_particle_type[particleType] == "Muons":
                    p_pred_mu = p_pred
                    p_error_mu = p_error

                if list_particle_type[particleType] == "Protons":
                    p_pred_p = p_pred
                    p_error_p = p_error

    legend = ax.legend(fontsize = 14)
    ax.set_title("WCTE BeamTest23 Run %i - momentum %s MeV/c,  n = %.3f, %s, pBunch = %.3e"%(run, momentum, refIndex, targetType[isBeryllium], probaBunch), fontsize = 16, weight = 'bold')
    ax.grid()
    ax.set_xlabel("TOF (ns)", fontsize=14)
    ax.set_ylabel("Number of particles/%.2fns"%widthBin, fontsize=14)
    ax.set_yscale('log')
    plt.legend()
    plt.savefig("final_results/WCTEBeamTest23_Run%i_p%i_TOF.pdf"%(run, momentum))
    plt.savefig("final_results/WCTEBeamTest23_Run%i_p%i_TOF.png"%(run, momentum))
    plt.show()
    # plt.close()


    ################################################################################
    ####################### Check the selection in LG & TOF  ######################
    ################################################################################
    fig, ax = plt.subplots(figsize = (16,9))

    ax.set_xlim(TOF_min-7, TOF_max+20)

    ax.annotate(selectionText, xy=(0.75, 0.7), xytext=(-20, -15), fontsize=14,
    xycoords='axes fraction', textcoords='offset points',
    bbox=dict(facecolor='white', alpha=0.8),
    horizontalalignment='center', verticalalignment='center')

    for particleType in range(len(list_data_frames)):
        #careful we overlay the selection on the other particles
        if particleType == 0:
            if thereIsProtons:
                #to also plot the deuterium
                ax.scatter(np.array(list_TOF[particleType]), np.array(list_LG[particleType]), marker = "x", color = 'lightgray', label = '%s w/ 1 particle: %i (%.3f percent of all events)'%(list_particle_type[particleType],list_particle_numbers[particleType],fractionPass1ParticleVeto * fractionPassNanVeto * 100), alpha = 0.3)
            else:
                ax.scatter(np.array(list_TOF[particleType][0]), np.array(list_LG[particleType][0]), marker = "x", color = 'lightgray', label = '%s w/ 1 particle: %i (%.3f percent of all events)'%(list_particle_type[particleType],list_particle_numbers[particleType],fractionPass1ParticleVeto * fractionPassNanVeto * 100), alpha = 0.3)

        elif list_particle_type[particleType] == "Electrons":
            ax.scatter(np.array(list_TOF[particleType]), np.array(list_LG[particleType]), marker = "x", color = 'darkgray', label = '%s: %i (%.3f percent of events w/ 1 particle)'%(list_particle_type[particleType],list_particle_numbers[particleType],list_particle_numbers[particleType]/nParticles * 100), alpha = 0.3)

        else:
            if (list_particle_type[particleType] == "Protons" and (not thereIsProtons)):
                    break
            ax.scatter(np.array(list_TOF[particleType]), np.array(list_LG[particleType]), marker = "x", label = '%s: %i (%.3f percent of events w/ 1 particle)'%(list_particle_type[particleType],list_particle_numbers[particleType],list_particle_numbers[particleType]/nParticles * 100), alpha = 0.3)

    legend = ax.legend(fontsize = 14)
    ax.set_title("WCTE BeamTest23 Run %i - momentum %s MeV/c,  n = %.3f, %s, pBunch = %.3e"%(run, momentum, refIndex, targetType[isBeryllium], probaBunch), fontsize = 16, weight = 'bold')
    ax.grid()
    ax.set_xlabel("TOF (ns)", fontsize=14)
    ax.set_ylabel("Lead glass peak integrated charge (PE)", fontsize=14)
    plt.savefig("final_results/WCTEBeamTest23_Run%i_p%i_selectionTOFlg.pdf"%(run, momentum))
    plt.savefig("final_results/WCTEBeamTest23_Run%i_p%i_selectionTOFlg.png"%(run, momentum))
    # plt.close()
    plt.show()

    ################################################################################
    ##check the selection in ACTs, second figure
    ################################################################################
    fig, ax = plt.subplots(figsize = (16,9))
    ax.annotate(selectionText, xy=(0.75, 0.7), xytext=(-20, -15), fontsize=14,
    xycoords='axes fraction', textcoords='offset points',
    bbox=dict(facecolor='white', alpha=0.8),
    horizontalalignment='center', verticalalignment='center')


    plotSelection(df_all, ACTlinearA, ACTlinearB, piMuBorderACT, ACTlower, thereIsProtons, protonsTOFCut, horizontal_el, max_ACT23 = 300, max_ACT1 = 60, plot_events = False, ax = ax)
    for particleType in range(len(list_data_frames)):
        #careful we overlay the selection on the other particles
        if particleType == 0:
            ax.scatter(np.array(list_ACT1[particleType][0]), np.array(list_ACT23[particleType][0]), marker = "x", color = 'lightgray', label = '%s w/ 1 particle: %i (%.3f percent of all events)'%(list_particle_type[particleType],list_particle_numbers[particleType],fractionPass1ParticleVeto * fractionPassNanVeto * 100))

        elif list_particle_type[particleType] == "Electrons":
            ax.scatter(np.array(list_ACT1[particleType]), np.array(list_ACT23[particleType]), marker = "x", color = 'darkgray', label = '%s: %i (%.3f percent of events w/ 1 particle)'%(list_particle_type[particleType],list_particle_numbers[particleType],list_particle_numbers[particleType]/nParticles * 100))

        else:
            if (list_particle_type[particleType] == "Protons" and (not thereIsProtons)):
                    break
            ax.scatter(np.array(list_ACT1[particleType]), np.array(list_ACT23[particleType]), marker = "x", label = '%s: %i (%.3f percent of events w/ 1 particle)'%(list_particle_type[particleType],list_particle_numbers[particleType],list_particle_numbers[particleType]/nParticles * 100))

    legend = ax.legend(fontsize = 14)
    ax.set_title("WCTE BeamTest23 Run %i - momentum %sMeV/c,  n = %.3f, %s, pBunch = %.3e"%(run, momentum, refIndex, targetType[isBeryllium], probaBunch), fontsize = 16, weight = 'bold')
    ax.grid()
    ax.set_xlabel("ACT1 (PE)", fontsize=14)
    ax.set_ylabel("ACT23 (PE)", fontsize=14)
    plt.savefig("final_results/WCTEBeamTest23_Run%i_p%i_selectionACTs.pdf"%(run, momentum))
    plt.savefig("final_results/WCTEBeamTest23_Run%i_p%i_selectionACTs.png"%(run, momentum))
    plt.show()
    # plt.close()


    ################################################################################
    ##check the selection in ACTs, second figure
    ################################################################################
    fig, ax = plt.subplots(figsize = (16,9))
    ax.annotate(selectionText, xy=(0.75, 0.7), xytext=(-20, -15), fontsize=14,
    xycoords='axes fraction', textcoords='offset points',
    bbox=dict(facecolor='white', alpha=0.8),
    horizontalalignment='center', verticalalignment='center')


    plotSelection(df_all, ACTlinearA, ACTlinearB, piMuBorderACT, ACTlower, thereIsProtons, protonsTOFCut, horizontal_el, max_ACT23 = 300, max_ACT1 = 500, plot_events = False, ax = ax)
    for particleType in range(len(list_data_frames)):
        #careful we overlay the selection on the other particles
        if particleType == 0:
            ax.scatter(np.array(list_ACT0[particleType][0]), np.array(list_ACT23[particleType][0]), marker = "x", color = 'lightgray', label = '%s w/ 1 particle: %i (%.3f percent of all events)'%(list_particle_type[particleType],list_particle_numbers[particleType],fractionPass1ParticleVeto * fractionPassNanVeto * 100))

        elif list_particle_type[particleType] == "Electrons":
            ax.scatter(np.array(list_ACT0[particleType]), np.array(list_ACT23[particleType]), marker = "x", color = 'darkgray', label = '%s: %i (%.3f percent of events w/ 1 particle)'%(list_particle_type[particleType],list_particle_numbers[particleType],list_particle_numbers[particleType]/nParticles * 100))

        else:
            if (list_particle_type[particleType] == "Protons" and (not thereIsProtons)):
                    break
            ax.scatter(np.array(list_ACT0[particleType]), np.array(list_ACT23[particleType]), marker = "x", label = '%s: %i (%.3f percent of events w/ 1 particle)'%(list_particle_type[particleType],list_particle_numbers[particleType],list_particle_numbers[particleType]/nParticles * 100))

    legend = ax.legend(fontsize = 14)
    ax.set_title("WCTE BeamTest23 Run %i - momentum %sMeV/c,  n = %.3f, %s, pBunch = %.3e"%(run, momentum, refIndex, targetType[isBeryllium], probaBunch), fontsize = 16, weight = 'bold')
    ax.grid()
    ax.set_xlabel("ACT0 (PE)", fontsize=14)
    ax.set_ylabel("ACT23 (PE)", fontsize=14)
    # plt.savefig("final_results/WCTEBeamTest23_Run%i_p%i_selectionACTs.pdf"%(run, momentum))
    # plt.savefig("final_results/WCTEBeamTest23_Run%i_p%i_selectionACTs.png"%(run, momentum))
    plt.show()
    plt.close()


    ######################### ACT0 vs ACT1 ##################
    # plotSelection(df_all, ACTlinearA, ACTlinearB, piMuBorderACT, ACTlower, thereIsProtons, protonsTOFCut, horizontal_el, max_ACT23 = 300, max_ACT1 = 500, plot_events = False, ax = ax)
    ig, ax = plt.subplots(figsize = (16,9))
    ax.annotate(selectionText, xy=(0.75, 0.7), xytext=(-20, -15), fontsize=14,
    xycoords='axes fraction', textcoords='offset points',
    bbox=dict(facecolor='white', alpha=0.8),
    horizontalalignment='center', verticalalignment='center')
    for particleType in range(len(list_data_frames)):
        #careful we overlay the selection on the other particles
        if particleType == 0:
            ax.scatter(np.array(list_ACT0[particleType][0]), np.array(list_ACT1[particleType][0]), marker = "x", color = 'lightgray', label = '%s w/ 1 particle: %i (%.3f percent of all events)'%(list_particle_type[particleType],list_particle_numbers[particleType],fractionPass1ParticleVeto * fractionPassNanVeto * 100))

        elif list_particle_type[particleType] == "Electrons":
            ax.scatter(np.array(list_ACT0[particleType]), np.array(list_ACT1[particleType]), marker = "x", color = 'darkgray', label = '%s: %i (%.3f percent of events w/ 1 particle)'%(list_particle_type[particleType],list_particle_numbers[particleType],list_particle_numbers[particleType]/nParticles * 100))

        else:
            if (list_particle_type[particleType] == "Protons" and (not thereIsProtons)):
                    break
            ax.scatter(np.array(list_ACT0[particleType]), np.array(list_ACT1[particleType]), marker = "x", label = '%s: %i (%.3f percent of events w/ 1 particle)'%(list_particle_type[particleType],list_particle_numbers[particleType],list_particle_numbers[particleType]/nParticles * 100))

    legend = ax.legend(fontsize = 14)
    ax.set_title("WCTE BeamTest23 Run %i - momentum %sMeV/c,  n = %.3f, %s, pBunch = %.3e"%(run, momentum, refIndex, targetType[isBeryllium], probaBunch), fontsize = 16, weight = 'bold')
    ax.grid()
    ax.set_xlabel("ACT0 (PE)", fontsize=14)
    ax.set_ylabel("ACT1 (PE)", fontsize=14)
    # plt.savefig("final_results/WCTEBeamTest23_Run%i_p%i_selectionACTs.pdf"%(run, momentum))
    # plt.savefig("final_results/WCTEBeamTest23_Run%i_p%i_selectionACTs.png"%(run, momentum))
    plt.show()
    # plt.close()



    #######################################################################################
    ################## check the purity of the different selections  ######################
    #######################################################################################

    #Start with pion/muon separation
    fig, ax = plt.subplots(figsize = (16,9))
    ax.annotate(selectionText, xy=(0.75, 0.7), xytext=(-20, -15), fontsize=14,
    xycoords='axes fraction', textcoords='offset points',
    bbox=dict(facecolor='white', alpha=0.8),
    horizontalalignment='center', verticalalignment='center')

    ACT23_min = 0
    ACT23_max = horizontal_el + 30
    bin_width = 2 #nb of PE in ACT23
    nBins = int((ACT23_max-ACT23_min)/bin_width)

    for particleType in range(len(list_data_frames)):
        #careful we overlay the selection on the other particles
        if list_particle_type[particleType] != "All events":
            #and list_particle_type[particleType] != "Electrons" :

            counts, edges = np.histogram(np.array(list_ACT23[particleType]), bins=nBins, range=(ACT23_min, ACT23_max), density=False)

            bin_centers = (edges[:-1] + edges[1:]) / 2

            plotting_x = np.linspace(min(bin_centers), max(bin_centers), 10*nBins)

            errors = np.sqrt(counts)# / len(np.array(list_ACT23[particleType])))

            meanACT23 = np.array(list_ACT23[particleType]).mean()
            stdACT23 = np.array(list_ACT23[particleType]).std()

            ax.errorbar(bin_centers, counts, yerr=errors, fmt='o', label = '%s: ACT23 = %.2f +/- %.2f (%.3f percent of events w/ 1 particle)'%(list_particle_type[particleType], meanACT23, stdACT23,list_particle_numbers[particleType]/nParticles * 100))

            initial_guess = [len(list_ACT23[particleType]), meanACT23, stdACT23]
            params, covariance = curve_fit(gaussian, bin_centers, counts, p0=initial_guess)
            amplitude, mean, std_dev = params

            y_curve = gaussian(plotting_x, amplitude, mean, std_dev)
            ax.plot(plotting_x, y_curve, '--', color = 'lightgray')

            if list_particle_type[particleType] == "Muons":
               amplitude_mu, mean_mu, std_dev_mu =  amplitude, mean, std_dev
            if list_particle_type[particleType] == "Pions":
               amplitude_pi, mean_pi, std_dev_pi =  amplitude, mean, std_dev


    ax.set_xlabel("ACT23 (PE)", fontsize=14)
    ax.set_ylabel("Number of particles collected/%.2f PE"%bin_width, fontsize=14)
    plt.grid()
    plt.legend(fontsize=14)
    ax.set_title("WCTE BeamTest23 Run %i - momentum %s MeV/c,  n = %.3f, %s, pBunch = %.3e"%(run, momentum, refIndex, targetType[isBeryllium], probaBunch), fontsize = 16, weight = 'bold')
    plt.savefig("final_results/WCTEBeamTest23_Run%i_p%i_purityChecks.pdf"%(run, momentum))
    plt.savefig("final_results/WCTEBeamTest23_Run%i_p%i_purityChecks.png"%(run, momentum))
    plt.show()

    ###################################################################
    #Here calculate purity for pions
    ##################################################################


    ACT23_min = -10
    window_low_tot = ACT23_min
    window_high_tot = ACT23_max
    muon_efficiency = []
    muon_purity = []
    window_low = []


    target = 0.99
    target2 = 0.9999
    target3 = 0.99999
    target4 = 0.999999
    target5 = 0.9999999
    a, b, c, d, e = 0, 0, 0, 0, 0

    n_muon_tot = quad(gaussian, window_low_tot, window_high_tot, args = (amplitude_mu, mean_mu, std_dev_mu))[0]  *(nBins)/ (ACT23_max-ACT23_min)
    n_pion_tot = quad(gaussian, window_low_tot, window_high_tot, args = (amplitude_pi, mean_pi, std_dev_pi))[0]  *(nBins)/ (ACT23_max-ACT23_min)

    for window_high_pi in np.arange(window_high_tot-0.25, window_low_tot+0.25, -0.05):
        n_muon_windowM = quad(gaussian, window_low_tot, window_high_pi, args = (amplitude_mu, mean_mu, std_dev_mu))[0]  *(nBins)/ (ACT23_max-ACT23_min)

        n_pion_windowM = quad(gaussian, window_low_tot, window_high_pi, args = (amplitude_pi, mean_pi, std_dev_pi))[0]  *(nBins)/ (ACT23_max-ACT23_min)

        purity = (n_pion_windowM)/(n_pion_windowM+n_muon_windowM)
        efficiency = n_pion_windowM/n_pion_tot

        muon_efficiency.append(efficiency)
        muon_purity.append(purity)

        if purity>target and a == 0:
            print("Pion contamination from muons of %.2e corresponds to an efficiency of %.3f percent, cut at %.2f PE"%((1-target), efficiency*100, window_high_pi))
            cut_T1 = window_high_pi
            eff_T1 = efficiency
            nMuonsInSample_T1 = n_muon_windowM
            nPionsInSample_T1 = n_pion_windowM
            a = 1

        if purity>target2 and b == 0:
            print("Pion contamination from muons of %.2e corresponds to an efficiency of %.3f percent, cut at %.2f PE"%((1-target2), efficiency*100, window_high_pi))
            cut_T2 = window_high_pi
            eff_T2 = efficiency
            nMuonsInSample_T2 = n_muon_windowM
            nPionsInSample_T2 = n_pion_windowM
            b = 1

        if purity>target3 and c == 0:
            print("Pion contamination from muons of %.2e corresponds to an efficiency of %.3f percent, cut at %.2f PE"%((1-target3), efficiency*100, window_high_pi))
            cut_T3 = window_high_pi
            eff_T3 = efficiency
            nMuonsInSample_T3 = n_muon_windowM
            nPionsInSample_T3 = n_pion_windowM
            c= 1

        if purity>target4 and d == 0:
            print("Pion contamination from muons of %.2e corresponds to an efficiency of %.3f percent, cut at %.2f PE"%((1-target4), efficiency*100, window_high_pi))
            cut_T4 = window_high_pi
            eff_T4 = efficiency
            nMuonsInSample_T4 = n_muon_windowM
            nPionsInSample_T4 = n_pion_windowM
            d = 1

        if purity>target5 and e == 0:
            print("Pion contamination from muons of %.2e corresponds to an efficiency of %.3f percent, cut at %.2f PE"%((1-target5), efficiency*100, window_high_pi))
            cut_T5 = window_high_pi
            eff_T5 = efficiency
            nMuonsInSample_T5 = n_muon_windowM
            nPionsInSample_T5 = n_pion_windowM
            e = 1

        window_low.append(window_high_pi)

    fig, ax1 = plt.subplots(figsize = (16,9))
    # ax1.set_title("Muon Purity and Efficiency - Run %s"%run)
    ax1.set_xlabel('ACT23 cut (PE)', fontsize = 14)
    ax1.set_ylabel('Purity', color = 'green', fontsize = 14)
    ax1.plot(window_low, muon_purity, color = 'green')
    # ax1.plot([best_cut_muons, best_cut_muons], [0,1], 'k--', label = '99.99percent muon purity \n corresponds to %.4f muon efficiency \n ACT23 cut %.2fpe'%(best_efficiency_muons, best_cut_muons))
    ax1.grid()
    ax1.tick_params(axis ='y', labelcolor = 'green')
    ax1.set_title("WCTE BeamTest23 Run %i - momentum %s MeV/c \nn = %.3f, %s, pBunch = %.3e"%(run, momentum, refIndex, targetType[isBeryllium], probaBunch), fontsize = 16, weight = 'bold')

    ax2 = ax1.twinx()
    ax2.set_ylabel('Efficiency', color = 'red')
    if a == 1:
        ax2.plot([cut_T1, cut_T1], [0, 1], '--', color = 'lightgray', label = '%.1e contamination in pion sample by muons: efficiency %.3f percent\n nMuonsInSample: %.3e nPionsInSample: %.3e'%(1-target, eff_T1 * 100, nMuonsInSample_T1, nPionsInSample_T1))

    if b == 1:
        ax2.plot([cut_T2, cut_T2], [0, 1], '-', color = 'lightgray', label = '%.1e contamination in pion sample by muons: efficiency %.3f percent\n nMuonsInSample: %.3e nPionsInSample: %.3e'%(1-target2, eff_T2 * 100, nMuonsInSample_T2, nPionsInSample_T2))

    if c == 1:
        ax2.plot([cut_T3, cut_T3], [0, 1], '--', color = 'darkgray', label = '%.1e contamination in pion sample by muons: efficiency %.3f percent\n nMuonsInSample: %.3e nPionsInSample: %.3e'%(1-target3, eff_T3 * 100, nMuonsInSample_T3, nPionsInSample_T3))

    if d == 1:
        ax2.plot([cut_T4, cut_T4], [0, 1], '--', color = 'black', label = '%.1e contamination in pion sample by muons: efficiency %.3f percent\n nMuonsInSample: %.3e nPionsInSample: %.3e'%(1-target4, eff_T4 * 100, nMuonsInSample_T4, nPionsInSample_T4))

    if e == 1:
        ax2.plot([cut_T5, cut_T5], [0, 1], '-', color = 'black', label = '%.1e contamination in pion sample by muons: efficiency %.3f percent\n nMuonsInSample: %.3e nPionsInSample: %.3e'%(1-target5, eff_T5 * 100, nMuonsInSample_T5, nPionsInSample_T5))

    ax2.plot(window_low, muon_efficiency, color = 'red')
    ax2.tick_params(axis ='y', labelcolor = 'red')
    ax2.legend(fontsize = 14)
    plt.savefig("final_results/WCTEBeamTest23_Run%i_p%i_PionPurityCurves.pdf"%(run, momentum))
    plt.savefig("final_results/WCTEBeamTest23_Run%i_p%i_PionPurityCurves.png"%(run, momentum))
    plt.show()



    ####### now contamination between electrons and muons


    list_elH = [] #and we'll just add them

    elH_max, elH_min = 0, 0


    for particleType in range(len(list_data_frames)):
        h_ACTs_diag = distToLinearCut(np.array(list_ACT1[particleType]), np.array(list_ACT23[particleType]), ACTlinearA, ACTlinearB)
        h_ACTs_hori = distToLinearCut(np.array(list_ACT1[particleType]), np.array(list_ACT23[particleType]), 0, horizontal_el)
        #TODO: improve this !!!!
        minimal_distance = np.where(h_ACTs_diag >= h_ACTs_hori, h_ACTs_hori, h_ACTs_diag)
        list_elH.append(minimal_distance)

        if minimal_distance.min() < elH_min:
            elH_min = -10#minimal_distance.min()

        if minimal_distance.max() > elH_max:
            elH_max = minimal_distance.max()


    bin_width = 0.5#nb of PE in ACT23
    nBins = int((elH_max-elH_min)/bin_width)
    fig, ax = plt.subplots(figsize = (16,9))
    ax.annotate(selectionText, xy=(0.75, 0.7), xytext=(-20, -15), fontsize=14,
    xycoords='axes fraction', textcoords='offset points',
    bbox=dict(facecolor='white', alpha=0.8),
    horizontalalignment='center', verticalalignment='center')



    for particleType in range(len(list_data_frames)):
        #careful we overlay the selection on the other particles
        if list_particle_type[particleType] != "All events" and list_particle_type[particleType] != "Protons":

            counts, edges = np.histogram(np.array(list_elH[particleType]), bins=nBins, range=(elH_min, elH_max), density=False)

            print("Counts ", counts)
            bin_centers = (edges[:-1] + edges[1:]) / 2
            plotting_x = np.linspace(min(bin_centers), max(bin_centers), 10*nBins)
            errors = np.sqrt(counts)# / len(np.array(list_ACT23[particleType])))
            mean_elH = np.array(list_elH[particleType]).mean()
            std_elH = np.array(list_elH[particleType]).std()
            ax.errorbar(bin_centers, counts, yerr=errors, fmt='o', label = '%s: Distance to e veto line = %.2f +/- %.2f (%.3f percent of events w/ 1 particle)'%(list_particle_type[particleType], mean_elH, std_elH,list_particle_numbers[particleType]/nParticles * 100))
            initial_guess = [len(list_elH[particleType]), mean_elH, std_elH]
            params, covariance = curve_fit(gaussian, bin_centers, counts, p0=initial_guess)
            amplitude, mean, std_dev = params
            y_curve = gaussian(plotting_x, amplitude, mean, std_dev)
            ax.plot(plotting_x, y_curve, '--', color = 'darkgray')
            if list_particle_type[particleType] == "Muons":
               amplitude_mu, mean_mu, std_dev_mu =  amplitude, mean, std_dev
            if list_particle_type[particleType] == "Pions":
               amplitude_pi, mean_pi, std_dev_pi =  amplitude, mean, std_dev
            if list_particle_type[particleType] == "Electrons":
               amplitude_e, mean_e, std_dev_e =  amplitude, mean, std_dev


            # ax.hist(np.array(list_ACT1[particleType]), np.array(list_ACT23[particleType]), marker = "x", color = 'darkgray', label = '%s: %i (%.3f percent of events w/ 1 particle)'%(list_particle_type[particleType],list_particle_numbers[particleType],list_particle_numbers[particleType]/nParticles * 100))

    ax.set_xlabel("Distance to the electron cut line (PE)", fontsize=14)
    ax.set_ylabel("Number of particles collected/%.2f PE"%bin_width, fontsize=14)
    plt.grid()
    plt.legend(fontsize=14)
    ax.set_title("WCTE BeamTest23 Run %i - momentum %s MeV/c,  n = %.3f, %s, pBunch = %.3e"%(run, momentum, refIndex, targetType[isBeryllium], probaBunch), fontsize = 16, weight = 'bold')
    plt.savefig("final_results/WCTEBeamTest23_Run%i_p%i_electronCheck.pdf"%(run, momentum))
    plt.savefig("final_results/WCTEBeamTest23_Run%i_p%i_electronCheck.png"%(run, momentum))
    plt.show()


    # Here calculate purity of the pion sample

    window_low_tot = elH_min
    window_high_tot = elH_max
    muon_efficiency = []
    muon_purity = []
    window_low = []

    target = 0.99
    target2 = 0.9999
    target3 = 0.99999
    target4 = 0.999999
    target5 = 0.9999999
    a, b, c, d, e = 0, 0, 0, 0, 0

    n_muon_tot = quad(gaussian, window_low_tot, window_high_tot, args = (amplitude_mu, mean_mu, std_dev_mu))[0]  *(nBins)/ (elH_max-elH_min)
    for window_low_mu in np.arange(window_low_tot+0.25, window_high_tot-0.25, 0.25):
        n_muon_windowM = quad(gaussian, window_low_mu, window_high_tot, args = (amplitude_mu, mean_mu, std_dev_mu))[0] *(nBins)/ (elH_max-elH_min)

        n_electron_windowM = quad(gaussian, window_low_mu, window_high_tot, args = (amplitude_e, mean_e, std_dev_e))[0]  * (nBins)/ (elH_max-elH_min)
        purity = (n_muon_windowM)/(n_electron_windowM+n_muon_windowM)
        efficiency = n_muon_windowM/n_muon_tot
        muon_efficiency.append(efficiency)
        muon_purity.append(purity)

        if purity>target and a == 0:
            print("Muon contamination from pions of %.2e corresponds to an efficiency of %.3f percent, cut at %.2f PE"%((1-target), efficiency*100, window_low_mu))
            cut_T1 = window_low_mu
            eff_T1 = efficiency
            nMuonsInSample_T1 = n_muon_windowM
            nElectronsInSample_T1 = n_electron_windowM
            a = 1

        if purity>target2 and b == 0:
            print("Muon contamination from pions of %.2e corresponds to an efficiency of %.3f percent, cut at %.2f PE"%((1-target2), efficiency*100, window_low_mu))
            cut_T2 = window_low_mu
            eff_T2 = efficiency
            nMuonsInSample_T2 = n_muon_windowM
            nElectronsInSample_T2 = n_electron_windowM
            b = 1

        if purity>target3 and c == 0:
            print("Muon contamination from pions of %.2e corresponds to an efficiency of %.3f percent, cut at %.2f PE"%((1-target3), efficiency*100, window_low_mu))
            cut_T3 = window_low_mu
            eff_T3 = efficiency
            nMuonsInSample_T3 = n_muon_windowM
            nElectronsInSample_T3 = n_electron_windowM
            c= 1

        if purity>target4 and d == 0:
            print("Muon contamination from pions of %.2e corresponds to an efficiency of %.3f percent, cut at %.2f PE"%((1-target4), efficiency*100, window_low_mu))
            cut_T4 = window_low_mu
            eff_T4 = efficiency
            nMuonsInSample_T4 = n_muon_windowM
            nElectronsInSample_T4 = n_electron_windowM
            d = 1

        if purity>target5 and e == 0:
            print("Muon contamination from pions of %.2e corresponds to an efficiency of %.3f percent, cut at %.2f PE"%((1-target5), efficiency*100, window_low_mu))
            cut_T5 = window_low_mu
            eff_T5 = efficiency
            nMuonsInSample_T5 = n_muon_windowM
            nElectronsInSample_T5 = n_electron_windowM
            e = 1

        window_low.append(window_low_mu)


    fig, ax1 = plt.subplots(figsize = (16,9))
    # ax1.set_title("Muon Purity and Efficiency - Run %s"%run)
    ax1.set_xlabel('ACT23 cut (PE)', fontsize = 14)
    ax1.set_ylabel('Purity', color = 'green', fontsize = 14)
    ax1.plot(window_low, muon_purity, color = 'green')
    # ax1.plot([best_cut_muons, best_cut_muons], [0,1], 'k--', label = '99.99percent muon purity \n corresponds to %.4f muon efficiency \n ACT23 cut %.2fpe'%(best_efficiency_muons, best_cut_muons))
    ax1.grid()
    ax1.tick_params(axis ='y', labelcolor = 'green')
    ax1.set_title("WCTE BeamTest23 Run %i - momentum %s MeV/c \nn = %.3f, %s, pBunch = %.3e"%(run, momentum, refIndex, targetType[isBeryllium], probaBunch), fontsize = 16, weight = 'bold')

    if a == 1:
        ax2 = ax1.twinx()
        ax2.set_ylabel('Efficiency', color = 'red')
        ax2.plot([cut_T1, cut_T1], [0, 1], '--', color = 'lightgray', label = '%.1e contamination in muon sample by electrons: efficiency %.3f percent\n nMuonsInSample: %.3e nElectronsInSample: %.3e'%(1-target, eff_T1 * 100, nMuonsInSample_T1, nElectronsInSample_T1))

    if b == 1:
        ax2.plot([cut_T2, cut_T2], [0, 1], '-', color = 'lightgray', label = '%.1e contamination in muon sample by electron: efficiency %.3f percent\n nMuonsInSample: %.3e nElectronsInSample: %.3e'%(1-target2, eff_T2 * 100, nMuonsInSample_T2, nElectronsInSample_T2))

    if c == 1:
        ax2.plot([cut_T3, cut_T3], [0, 1], '--', color = 'darkgray', label = '%.1e contamination in muon sample by electron: efficiency %.3f percent\n nMuonsInSample: %.3e nElectronsInSample: %.3e'%(1-target3, eff_T3 * 100, nMuonsInSample_T3, nElectronsInSample_T3))

    if d == 1 :
        ax2.plot([cut_T4, cut_T4], [0, 1], '--', color = 'black', label = '%.1e contamination in muon sample by electron: efficiency %.3f percent\n nMuonsInSample: %.3e nElectronsInSample: %.3e'%(1-target4, eff_T4 * 100, nMuonsInSample_T4, nElectronsInSample_T4))

    if e == 1:
        ax2.plot([cut_T5, cut_T5], [0, 1], '-', color = 'black', label = '%.1e contamination in muon sample by electron: efficiency %.3f percent\n nMuonsInSample: %.3e nElectronsInSample: %.3e'%(1-target5, eff_T5 * 100, nMuonsInSample_T5, nElectronsInSample_T5))

    ax2.plot(window_low, muon_efficiency, color = 'red')
    ax2.tick_params(axis ='y', labelcolor = 'red')
    ax2.legend(fontsize = 14)
    plt.savefig("final_results/WCTEBeamTest23_Run%i_p%i_purityECurves.pdf"%(run, momentum))
    plt.savefig("final_results/WCTEBeamTest23_Run%i_p%i_purityECurves.png"%(run, momentum))
    plt.show()


    ###############################################################################
    #and finally save the number of particles and also the fitting paraemters
    # so we do not have to do it again
    #################################################################################
    table = np.array([run, momentum, refIndex, TOTAL_spill, probaBunch, nParticles, nElectrons, nMuons, nPions, nProtons, nDeuterium, fractionPass1ParticleVeto, fractionPassNanVeto, ACTlinearA, ACTlinearB, piMuBorderACT, ACTlower, thereIsProtons, protonsTOFCut, horizontal_el, LGupper, isBeryllium, p_pred_mu, p_pred_pi, p_pred_p, p_error_mu, p_error_pi, p_error_p, Mean_ACT23_mu, Std_ACT23_mu, Mean_ACT23_pi, Std_ACT23_pi, Mean_ACT23_p, Std_ACT23_p])


    with open("numberParticles.txt", "a") as file:
        #first the header, only once
        file.write("run momentum refIndex nSpills probaBunch nParticles nElectrons nMuons nPions nProtons nDeuterium fractionPass1ParticleVeto fractionPassNanVeto ACTlinearA ACTlinearB piMuBorderACT ACTlower thereIsProtons protonsTOFCut horizontal_el LGupper bery\n")
        file.write("%i %i %.3f %i %.3e %i %i %i %i %i %i %.5f %.5f %.2f %.2f %.2f %2f %i %.2f %.2f %.2f %i %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n" % (table[0], table[1], table[2], table[3], table[4], table[5], table[6], table[7], table[8],table[9], table[10], table[11], table[12], table[13],table[14], table[15], table[16], table[17],table[18], table[19], table[20], table[21], table[22], table[23],table[24], table[25], table[26], table[27],  table[28], table[29],table[30], table[31], table[32], table[33]) )

    #exit the program for now
    raise end




if __name__ == "__main__":
    singlePE(sys.argv)




