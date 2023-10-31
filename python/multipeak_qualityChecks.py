import uproot as ur

import matplotlib.pyplot as plt

import numpy as np

import pandas as pd
import sys

import matplotlib.colors as colors

from scipy.optimize import curve_fit

from scipy.integrate import quad

############################### READ ME #################################################
#This is a very crude code to look at multi peaks events - tailored for run 373 for now
#strated by acraplet
#reads the new ntuple format and looks at the hit time, trying to make useful plots to
#understand how best to match the peaks and calibrate for the ACT diffrent cable length
#########################################################################################

#will look at the number of hits against T0F for a given event

# plt.style(['notebook', 'grid'])

def peakMatchInTOf(detectors_pd, all_hit_times, pmt, pmt1, pmt2, pmt3, FirstEvent = 0, LastEvent = 200, matchWindow = 1.5):
    #match window is the delay in ns around the first PMT hit time that we accept other PMT's hits in
    matchedCounter = 0
    notMatchedCounter = 0
    number_of_hits_per_event = []
    certified_hit_times = [] #this is where we save the mean hit times of interest

    for event in range(FirstEvent, LastEvent):
        hit = 0
        matchedCounterPerEvent = 0
        # print(event)
        while (not np.isnan(all_hit_times[pmt][hit][event]) and hit<max(detectors_pd[pmt]['nPeaks'])-1):
            # print("Event ", event, ' pmt: ', pmt, ' hit :', hit, ' time: ', all_hit_times[pmt][hit][event], 'ns')
            MatchedInPMT1 = False
            MatchedInPMT2 = False
            MatchedInPMT3 = False
            tPMT0 = all_hit_times[pmt][hit][event]
            tPMT1 = 0
            tPMT2 = 0
            tPMT3 = 0
            hit1=0
            while (not np.isnan(all_hit_times[pmt1][hit1][event])and hit1<max(detectors_pd[pmt1]['nPeaks'])-1):
                if abs(all_hit_times[pmt1][hit1][event]-all_hit_times[pmt][hit][event])<matchWindow:
                    MatchedInPMT1 = True
                    tPMT1 = all_hit_times[pmt1][hit1][event]
                    break
                hit1 += 1

            hit2=0
            while (not np.isnan(all_hit_times[pmt2][hit2][event])and hit2<max(detectors_pd[pmt2]['nPeaks'])-1):
                if abs(all_hit_times[pmt2][hit2][event]-all_hit_times[pmt][hit][event])<matchWindow:
                    MatchedInPMT2 = True
                    tPMT2 = all_hit_times[pmt2][hit2][event]
                    break
                hit2 += 1

            hit3=0
            while (not np.isnan(all_hit_times[pmt3][hit3][event]) and hit3<max(detectors_pd[pmt3]['nPeaks'])-1):
                if abs(all_hit_times[pmt3][hit3][event]-all_hit_times[pmt][hit][event])<matchWindow:
                    MatchedInPMT3 = True
                    tPMT3 = all_hit_times[pmt3][hit3][event]
                    break
                hit3 += 1

            if (MatchedInPMT1 and MatchedInPMT2 and MatchedInPMT3 and matchedCounterPerEvent == 0):
                # print("Macthed in all detectors")

                matchedCounter += 1
                matchedCounterPerEvent += 1
            else:
                # print('Not matched in all detectors')
                # certified_hit_times.append(0)
                notMatchedCounter += 1
            hit += 1
        number_of_hits_per_event.append(matchedCounterPerEvent)
        certified_hit_times.append((tPMT0+tPMT1+tPMT2+tPMT3)/4)
    return matchedCounter, notMatchedCounter, np.array(certified_hit_times), np.array(number_of_hits_per_event)

def peakMatching(argv):
    #actually peak matching
    root_filenames = argv[1]
    file = ur.open(root_filenames)
    detectors_pd = []

    print(file.keys())
    print(file['TOF03'].keys())

    all_hit_times = []
    # entries_split_accros_branches = False
#
# #for the new format when we have the entries split up accross branches
#     if entries_split_accros_branches:
#         list_detector = ["ACT0L"]# 'ACT0R', 'ACT1L', 'ACT1R', 'ACT2L', 'ACT2R', 'ACT3L', 'ACT3R', 'TOF00', 'TOF01', 'TOF02', 'TOF03', 'TOF10', 'TOF11', 'TOF12', 'TOF13', 'Hole0', 'Hole1', 'PbGlass']
#
#         for key in list_detector:
#             df = file[key].arrays(library="pd"[[], [], [], [], [], [], [], ..., [...], [0.00287], [], [], [], [6.84e-05], []]
# )
#             print(file[key])
#             print(df, df["SignalTime"])
#             # for detector in range(len(list_detector)):
#             #
#             #     if list_detector[detector] in key:
#             #         print(key)
#             #         #we need to merge all of the events that are split into ;1, ;2, ;3, ;4
#             #         #very slow process -> discuss with Nick, could we retort to the other format?
#             #         if len(detectors_pd)<=detector:
#             detectors_pd.append(df)
#             # print(detectors_pd['triggerTime'])
#                     # else:
#                     # #     # print(detectors_pd[detector][0], df)
#                     #      detectors_pd[detector] = detectors_pd[detector].append(df, ignore_index = True)
#
#         # print( detectors_pd[1])
#
#     else: #if we have the easier format where everything is in the same branch
    for key in file.keys()[:-1]:
        df = file[key].arrays(library="pd")
        detectors_pd.append(df)
            # detectors_pd_charge.append(df)

    # tree_array = file["ACT0L;1"].arrays()
    # # print(tree_array, 'HEYY')
    # print(tree_array["SignalTime"], 'First cycle: length: ', len(tree_array["SignalTime"]))
    # tree_array = file["ACT0L;2"].arrays()
    # # print(tree_array, 'HEYY')
    # print(tree_array["SignalTime"], 'Second cycle: length: ', len(tree_array["SignalTime"]))
    # tree_array = file["ACT0L;3"].arrays()
    # # print(tree_array, 'HEYY')
    # print(tree_array["SignalTime"], 'Third cycle: length: ', len(tree_array["SignalTime"]))
    #
    # raise end


    #extract the hit times
    for pmt in range(19):
        print(pmt)
        hit_times = pd.DataFrame(detectors_pd[pmt]['SignalTime'].values.tolist())
        detectors_pd[pmt] = pd.concat([detectors_pd[pmt],hit_times], axis=1)
        print(detectors_pd[pmt]['triggerTime'])


    nPoints = 1000
    offset = 0
    plt.plot(detectors_pd[0]['triggerTime'][offset:nPoints+offset]-detectors_pd[9]['triggerTime'][offset:nPoints+offset], marker = 'x')
    # fit_params = np.polyfit(np.array(np.array(detectors_pd[-1]['triggerTime'][offset:nPoints+offset])), detectors_pd[9]['triggerTime'][offset:nPoints+offset], 1)
    # p = np.poly1d(fit_params)


    lin_min = min(min(detectors_pd[-1]['triggerTime'][offset:nPoints+offset]), min(detectors_pd[9]['triggerTime'][offset:nPoints+offset]))
    lin_max = max(max(detectors_pd[-1]['triggerTime'][offset:nPoints+offset]), max(detectors_pd[9]['triggerTime'][offset:nPoints+offset]))
    x = np.linspace(lin_min, lin_max, 100)

    # plt.plot(x, p(x), 'r--', label = 'linear fit: %.3fx + %.3e'%(fit_params[0], fit_params[1]))
    # plt.plot([lin_min, lin_max], [lin_min, lin_max], 'k--')
    plt.ylabel('Trigger time for the ACT - Trigger time for the TOF')
    plt.xlabel('Event')
    plt.grid()
    plt.legend()
    plt.show()

    nPoints = 10000
    offset = 0
    plt.scatter(np.array(np.array(detectors_pd[0]['timeStamp'][offset:nPoints+offset])), detectors_pd[9]['timeStamp'][offset:nPoints+offset], marker = 'x')
    fit_params = np.polyfit(np.array(np.array(detectors_pd[0]['timeStamp'][offset:nPoints+offset])), detectors_pd[9]['timeStamp'][offset:nPoints+offset], 1)
    p = np.poly1d(fit_params)

    lin_min = min(min(detectors_pd[0]['timeStamp'][offset:nPoints+offset]), min(detectors_pd[9]['timeStamp'][offset:nPoints+offset]))
    lin_max = max(max(detectors_pd[0]['timeStamp'][offset:nPoints+offset]), max(detectors_pd[9]['timeStamp'][offset:nPoints+offset]))
    x = np.linspace(lin_min, lin_max, 100)

    plt.plot(x, p(x), 'r--', label = 'linear fit: %.3fx + %.3e'%(fit_params[0], fit_params[1]))
    # plt.plot([lin_min, lin_max], [lin_min, lin_max], 'k--')
    plt.ylabel('time Stamp for the ACTs')
    plt.xlabel('time Stamp for the TOFs')
    plt.grid()
    plt.legend()
    plt.show()
    # print(detectors_pd[pmt]['triggerTime'])

    for pmt in range(19):
        all_hit_times.append([])
        print("Pmt:", pmt, " max number of peaks is:", max(detectors_pd[pmt]['nPeaks']))
        for hit in range(max(detectors_pd[pmt]['nPeaks'])):
            all_hit_times[pmt].append(detectors_pd[pmt][hit])


    #now peak matching:
    for window in [2]:
        matchedCounter, notMatchedCounter, certified_hit_times1, number_of_hits_per_event1 = peakMatchInTOf(detectors_pd, all_hit_times,8, 9, 10, 11, 1000, 1100, window)
        if window == 2:
            plt.plot(number_of_hits_per_event1, 'rx', label = 'TOF0')
        else:
            plt.plot(window, matchedCounter, 'rx')

        matchedCounter, notMatchedCounter, certified_hit_times, number_of_hits_per_event = peakMatchInTOf(detectors_pd, all_hit_times,12, 13, 14, 15, 1000, 1100, window)
        TOFmatch = np.where(number_of_hits_per_event==number_of_hits_per_event1, number_of_hits_per_event, -1 )

        if window == 2:
            plt.plot(number_of_hits_per_event, 'bx', label = 'TOF1')
            plt.plot(TOFmatch, 'gx', label = 'Same number in TOF0 and TOF1')

        else:
            plt.plot(window, matchedCounter, 'bx')

    plt.xlabel('Event')
    plt.ylabel('Number of peak matched in the TOFs')
    plt.title('Number of peaks matched per Event when the window is %.1fns'%window)
    plt.grid()
    plt.ylim(-0.5, 3.5)
    plt.legend()
    plt.show()
    #cher

    plt.figure(1)
    plt.plot(certified_hit_times1 - certified_hit_times, 'bx')
    plt.show()


    # plt.hist(certified_hit_times1 - certified_hit_times)



    return 0



def multipeaksQualityCheck(argv):
    root_filenames = argv[1]
    file = ur.open(root_filenames)

    detectors_pd = []
    detectors_pd_charge = []
    detectors_pd_window_charge = []

    print(file.keys())
    print(file['TOF03;1'].keys())

    # onePeakToF = True

    n = 1

    df = file['TOF00'].arrays(library="pd")
    onePeakToF = (df['nPeaks']==n)
    df = file['TOF01'].arrays(library="pd")
    onePeakToF = (df['nPeaks']==n) & onePeakToF
    df = file['TOF02'].arrays(library="pd")
    onePeakToF = (df['nPeaks']==n) & onePeakToF
    df = file['TOF03'].arrays(library="pd")
    onePeakToF = (df['nPeaks']==n) & onePeakToF
    df = file['TOF10'].arrays(library="pd")
    onePeakToF = (df['nPeaks']==n) & onePeakToF

    #need to get rid of the accidentals in TOF10 -> pretty unlikely in this case of only one
    # onePeakToF = np.where(df['IntCharge'][:][0]<=0.015, False, onePeakToF)
    df = file['TOF11'].arrays(library="pd")
    onePeakToF = (df['nPeaks']==n) & onePeakToF
    df = file['TOF12'].arrays(library="pd")
    onePeakToF = (df['nPeaks']==n) & onePeakToF
    df = file['TOF13'].arrays(library="pd")
    onePeakToF = (df['nPeaks']==n) & onePeakToF


    # df = file['ACT0L'].arrays(library="pd")
    # onePeakToF = (df['nPeaks']==1) & onePeakToF





#read in each of the branches
    for key in file.keys()[:-1]:
        df = file[key].arrays(library="pd")
        #it is useful to look at only the case where we have one peak in all of the TOF
        df = df.reset_index()
        detectors_pd.append(df[onePeakToF])
        detectors_pd_charge.append(df[onePeakToF])
        detectors_pd_window_charge.append(df[onePeakToF])


    all_hit_times = []
    all_hit_charge = []
    all_hit_window_charge = []

    #make an array of [DataframePMT1, DataFramePMT2, ...] where DataFramePMT2 has columns 0,1, 2, 3... corresponding to the hit time of each hit, 0th 1st...
    for pmt in range(19):
        print(pmt)
        hit_times = pd.DataFrame(detectors_pd[pmt]['SignalTime'].values.tolist())
        hit_charge = pd.DataFrame(detectors_pd[pmt]['IntCharge'].values.tolist())
        hit_window_charge = pd.DataFrame(detectors_pd[pmt]['WindowIntCharge'].values.tolist())
        detectors_pd[pmt] = pd.concat([detectors_pd[pmt],hit_times], axis=1)
        detectors_pd_charge[pmt] = pd.concat([detectors_pd_charge[pmt],hit_charge], axis=1)
        detectors_pd_window_charge[pmt] = pd.concat([detectors_pd_window_charge[pmt],hit_window_charge], axis=1)
        # all_second_hit_times.append(hit_times[1])

    #here, reading in the dataframes to fill an array of all of the hits for each PMT (easier than to access the dataframe (a bit convoluted, I agree,
    for pmt in range(19):
        all_hit_times.append([])
        all_hit_charge.append([])
        all_hit_window_charge.append([])
        print("Pmt:", pmt, " max number of peaks is:", max(detectors_pd[pmt]['nPeaks']))
        #need the number of hits of the corresponding entry
        for hit in range(n): # int(max(detectors_pd[pmt]['nPeaks']) #int(min(max(detectors_pd[pmt]['nPeaks']), max(detectors_pd[12]['nPeaks'])
            print(hit)
            # print(
            all_hit_times[pmt].append(detectors_pd[pmt][hit])
            all_hit_charge[pmt].append(detectors_pd_charge[pmt][hit])
            all_hit_window_charge[pmt].append(detectors_pd_window_charge[pmt][hit])


    #to look at different detectors and plot them more easily
    i_min = 0 #the range of detectors of interest
    i_max = 8
    list_markers = ['o','v', '^', 's', '<', '>', 'X', 'd', 'H', 'h', 'D', 'P', 'p','*', 'o','v', '^', 's', '>', '<']
    PMTref = 11

    mean_first_hit_TOF1 = (detectors_pd[12][0] + detectors_pd[13][0] + detectors_pd[14][0] + detectors_pd[15][0])/4
    mean_first_hit_TOF0 = (detectors_pd[8][0] + detectors_pd[9][0] + detectors_pd[10][0] + detectors_pd[11][0])/4
    print(mean_first_hit_TOF1- mean_first_hit_TOF0)
    # print(all_hit_times)

    def gaussian(x, mean, amplitude, standard_deviation):
        return amplitude * np.exp( - (x - mean)**2 / (2*standard_deviation ** 2))

    def linear(x, slope, intercept):
        return x * slope + intercept

    #now we will compare the histogram of the difference between each detector and the reference: mean TOF
    for i in range(i_min, i_max):

        particle = 0

        hist = all_hit_times[i][particle]
        #TOF1 - TOFparticle
        particleTOF =  (all_hit_times[12][particle]+all_hit_times[13][particle]+all_hit_times[14][particle]+all_hit_times[15][particle] - (all_hit_times[8][particle] + all_hit_times[9][particle] + all_hit_times[10][particle] + all_hit_times[11][particle]))/4

        print(particleTOF, all_hit_times[8][particle])


        hist_charge = all_hit_charge[i][particle]
        leadGlass_charge = all_hit_charge[18][particle]
        hist_window_charge = all_hit_window_charge[i][particle]

        #
        # print(detectors_pd_window_charge[i]["WindowMinusPeak"], 'Window Minus peak Int charge ')
        # print(detectors_pd_window_charge[i][detectors_pd_window_charge[i]["WindowMinusPeak"])
        #



        hist_window_charge = all_hit_window_charge[i][particle]


        detectors_pd_window_charge[i]["WindowMinusPeak"] =  hist_window_charge -  hist_charge
        print(detectors_pd_window_charge[i][detectors_pd_window_charge[i]["WindowMinusPeak"] > 0.025])

        counter = 1

        plt.figure(2)
        plt.hist2d(hist_charge, hist_window_charge,  bins = [100, 100], range = [[0, 0.1], [0, 0.1]],  label = file.keys()[i], norm = colors.LogNorm())

        plt.figure(2)
        plt.grid()
        plt.title("Window integrated charge vs peak integrated charge, one hit in TOFs %s"%file.keys()[i])
        plt.colorbar(label="Number of triggers/(0.001V*0.001V)")
        plt.ylabel('Window integrated hit charge (V)')
        plt.xlabel('Peak integrated hit charge (V)')
        plt.show()

        plt.figure(3)
        plt.hist2d( detectors_pd[i]['PedestalSigma'], hist_window_charge -hist_charge,  bins = [100, 400], range = [[0, 0.1], [-0.2, 0.2]], label = file.keys()[i], norm = colors.LogNorm())
        plt.grid()
        plt.title("One hit in TOFs %s"%file.keys()[i])
        plt.colorbar(label="Number of triggers/(0.001V*0.001V)")
        plt.ylabel('Window integrated - Peak integrated hit charge (V)')
        plt.xlabel('PedestalSigma (V)')
        plt.show()

        plt.figure(4)
        plt.hist2d(hist-mean_first_hit_TOF1, hist_window_charge - hist_charge,  bins = [100, 200], range = [[-80, -10],[-0.1, 0.1]],  label = file.keys()[i], norm = colors.LogNorm())
        plt.grid()
        plt.title("%i hit in TOFs  %s"%(n, file.keys()[i]))
        plt.colorbar(label="Number of triggers/(0.5ns*0.002V)")
        plt.ylabel('Window integrated - Peak integrated hit charge (V)')
        plt.xlabel('HitTime - TOF1 Mean hit time(ns)')
        plt.show()


        plt.figure(5)
        plt.hist2d(particleTOF, hist_window_charge - hist_charge,  bins = [200, 200], range = [[10,13],[-0.1, 0.1]],  label = file.keys()[i], norm = colors.LogNorm())
        plt.grid()
        plt.title("%i hit in TOFs  %s"%(n, file.keys()[i]))
        plt.colorbar(label="Number of triggers/(0.1ns*0.002V)")
        plt.ylabel('Window integrated - Peak integrated hit charge (V)')
        plt.xlabel('Particle Time of flight(ns)')
        plt.show()


        plt.figure(6)
        plt.subplot(2,1,1)
        plt.hist2d(leadGlass_charge, hist_charge,  bins = [80, 70], range = [[-0.1, 0.7],[-0.04, 0.66]],  label = file.keys()[i], norm = colors.LogNorm())
        plt.grid()
        plt.title("%i hit in TOFs  %s"%(n, file.keys()[i]))
        plt.colorbar(label="Number of triggers/(0.01V*0.01V)")
        plt.ylabel('Peak integrated hit charge (V)')
        plt.xlabel('LeadGlass peak integrated charge')

        plt.subplot(2,1,2)
        plt.hist2d(leadGlass_charge, hist_window_charge,  bins = [80, 70], range = [[-0.1, 0.7],[-0.04, 0.66]],  label = file.keys()[i], norm = colors.LogNorm())
        plt.grid()
        # plt.title("%i hit in TOFs  %s"%(n, file.keys()[i]))
        plt.colorbar(label="Number of triggers/(0.01V*0.01V)")
        plt.ylabel('Window integrated integrated hit charge (V)')
        plt.xlabel('LeadGlass peak integrated charge')

        plt.figure(7)
        plt.subplot(2,1,1)
        plt.hist2d(particleTOF, hist_charge,  bins = [60, 70], range = [[10,16],[-0.04, 0.66]],  label = file.keys()[i], norm = colors.LogNorm())
        plt.grid()
        plt.title("%i hit in TOFs  %s"%(n, file.keys()[i]))
        plt.colorbar(label="Number of triggers/(0.1ns*0.01V)")
        plt.ylabel('Peak integrated hit charge (V)')
        plt.xlabel('Particle Time of flight(ns)')

        plt.subplot(2,1,2)
        plt.hist2d(particleTOF, hist_window_charge,  bins = [60, 70], range = [[10,16],[-0.04, 0.66]],  label = file.keys()[i], norm = colors.LogNorm())
        plt.grid()
        plt.title("%i hit in TOFs  %s"%(n, file.keys()[i]))
        plt.colorbar(label="Number of triggers/(0.1ns*0.01V)")
        plt.ylabel('Window integrated integrated hit charge (V)')
        plt.xlabel('Particle Time of flight(ns)')







        plt.show()







    raise end

def nPeakInToF(file, n):
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

def noHoleHit(file, conditionInitial):
    df = file['Hole0'].arrays(library="pd")
    condition = (df['nPeaks']==0) & conditionInitial
    print(condition, sum((df['nPeaks']==0)) )
    df = file['Hole1'].arrays(library="pd")
    condition = (df['nPeaks']==0) & condition
    print(condition, (df['nPeaks']==0) )
    return condition

def meanDistancePeaks(argv):
    root_filenames = argv[1]
    file = ur.open(root_filenames)
    nPeaksInToF = nPeakInToF(file, 2)
    detectors_pd = []

    for detector in range(len(file.keys()[:-1])):
            key = file.keys()[detector]
            print(key)
            df = file[key].arrays(library="pd")
            detectors_pd.append(df[nPeaksInToF])

    for pmt in range(len(file.keys()[:-1])):
        if pmt>=8 and pmt < 16:
            #all of the TOF distributions
            hit_times = pd.DataFrame(detectors_pd[pmt]['SignalTime'].values.tolist())
            detectors_pd[pmt] = pd.concat([detectors_pd[pmt],hit_times], axis=1)

            mean_delay = detectors_pd[pmt][1] - detectors_pd[pmt][0]
            print(mean_delay)
            plt.hist(mean_delay, bins = 55, range = (0,110), label = '%s'%file.keys()[pmt][:-1], alpha = 0.4)
    plt.legend()
    plt.title('%s'%root_filenames)
    plt.xlabel('Delay between two consecutive particles [peak detect] (ns)')
    plt.ylabel('Number of triggers/2ns')
    plt.grid()
    plt.show()

    #deviation accros the first hit - we know they are small from the 1sigma variation
    #but we can also select the case when they all have only 1 signal

    #print(detectors_pd[8]['SignalTime'])


    # print(file['TOF03;1']['nPeaks'])

    # T0F03 = file.arrays(library="pd")
    # print(T0F03[T0F03['nPeaks']!=1]['nPeaks'])
    # print(T0F03)

    # for i in range(len(file['TOF03;1']['nPeaks'].array())):
    #     print(file['TOF03;1']['nPeaks'][i])

def cutTOF(argv):
    root_filenames = argv[1]
    file = ur.open(root_filenames)
    nPeaksInToF = nPeakInToF(file, False)
    detectors_pd_el = []
    detectors_pd_NONel = []
    detectors_pd_el_wind = []
    detectors_pd_NONel_wind = []
    detectors_pd = []

    for detector in range(len(file.keys()[:-1])):
            key = file.keys()[detector]
            df = file[key].arrays(library="pd")
            detectors_pd.append(df[nPeaksInToF])

    detectors_pd_safe = detectors_pd.copy()
    detectors_pd_wind = detectors_pd.copy()

    for pmt in range(len(file.keys()[:-1])):
        # if pmt>=8 and pmt < 16:
            #all of the TOF distributions
        hit_times = pd.DataFrame(detectors_pd[pmt]['SignalTime'].values.tolist())
        intCharge = pd.DataFrame(detectors_pd_safe[pmt]['IntCharge'].values.tolist())
        detectors_pd_safe[pmt] = pd.concat([detectors_pd_safe[pmt],intCharge], axis=1)
        detectors_pd[pmt] = pd.concat([detectors_pd[pmt],hit_times], axis=1)

        intCharge = pd.DataFrame(detectors_pd_safe[pmt]['WindowIntCharge'].values.tolist())
        detectors_pd_wind[pmt] = pd.concat([detectors_pd_wind[pmt],intCharge], axis=1)

    mean_first_hit_TOF1 = (detectors_pd[12][0] + detectors_pd[13][0] + detectors_pd[14][0] + detectors_pd[15][0])/4
    mean_first_hit_TOF0 = (detectors_pd[8][0] + detectors_pd[9][0] + detectors_pd[10][0] + detectors_pd[11][0])/4

    tof = np.array(np.array(mean_first_hit_TOF1)-np.array(mean_first_hit_TOF0))
    print(tof)


    limCutTOF = 12.5

    for detector in range(len(file.keys()[:-1])):
        a = pd.Series(data = tof, name = 'TOF')
        detectors_pd_safe[detector] = pd.concat([detectors_pd_safe[detector],a], axis=1)
        detectors_pd_wind[detector] = pd.concat([detectors_pd_wind[detector],a], axis=1)

        detectors_pd_el.append(detectors_pd_safe[detector][detectors_pd_safe[detector]['TOF']<limCutTOF])
        detectors_pd_NONel.append(detectors_pd_safe[detector][detectors_pd_safe[detector]['TOF']>limCutTOF])

        detectors_pd_el_wind.append(detectors_pd_wind[detector][detectors_pd_wind[detector]['TOF']<limCutTOF])
        detectors_pd_NONel_wind.append(detectors_pd_wind[detector][detectors_pd_wind[detector]['TOF']>limCutTOF])


    det = 9
    print(detectors_pd_el[det])
    # plt.hist(detectors_pd_el[det][0], bins = 75, range = (0, 0.75), histtype = 'step', label = 'TOF < %.1fns'%limCutTOF)
    plt.hist(detectors_pd_NONel[det][0], bins = 50, range = (0, 1), histtype = 'step', linestyle = '--', label = 'TOF > %.1fns PEAK integrated'%limCutTOF)
    plt.hist(detectors_pd_NONel_wind[det][0] , bins = 50, range = (0, 1), histtype = 'step', label = 'TOF > %.1fns WINDOW integrated '%limCutTOF)
    plt.hist(detectors_pd_NONel[det]["MaxVoltage"] * (0.9/2) , bins = 50, range = (0, 1), histtype = 'step', label = 'TOF > %.1fns max voltage * (0.9/2)'%limCutTOF)
    plt.title("%s - Peak int charge" %file.keys()[det])
    plt.xlabel('Charge (nC)')
    plt.ylabel('Triggers/0.01nC')
    # plt.semilogy()
    plt.grid()
    plt.legend()
    # plt.hist(detectors_pd_el[det][1], bins = 100, range = (0, 0.75), histtype = 'step')
    # plt.hist(detectors_pd_NONel[det][1], bins = 100, range = (0, 0.75), histtype = 'step')
    plt.show()

    plt.title("PEAK integrated - TOF > %.1fns"%limCutTOF)

    plt.hist2d(detectors_pd_NONel[18][0], detectors_pd_NONel_wind[det][0]+detectors_pd_NONel_wind[det+1][0],  bins = (100, 100), range = ((0,0.25), (0,1)), norm = colors.LogNorm())
    plt.colorbar(label="Number of triggers")
    plt.xlabel(file.keys()[18])
    plt.ylabel("%s - WINDOW "%file.keys()[det][:-3])
    plt.show()

    plt.title("PEAK integrated - TOF < %.1fns"%limCutTOF)

    plt.hist2d(detectors_pd_el[18][0], detectors_pd_el[det][0]+detectors_pd_el[det+1][0],  bins = (100, 100), range = ((0,0.25), (0,1)), norm = colors.LogNorm())
    plt.colorbar(label="Number of triggers")
    plt.xlabel(file.keys()[18])
    plt.ylabel(file.keys()[det][:-3])
    plt.show()

    det = 0
    a = detectors_pd_safe[det][detectors_pd_safe[det]["nPeaks"] == 1]
    #check the peak height vs the integrated charge to see if there is an easy relationship
    plt.hist2d(a["MaxVoltage"], a[0], bins = (100,100), range= ((0,2), (0,.5)), label = 'Max voltage vs peakInt %s'%file.keys()[det], norm = colors.LogNorm())
    plt.xlabel('Max Voltage')
    plt.ylabel('Peak integrated charge')
    plt.title("%s - Run 502" %file.keys()[det])
    plt.colorbar(label="Number of triggers")
    plt.grid()
    plt.show()


    for det in range(0,8):
        nbins = 450
        plt.hist(detectors_pd[det]["MaxVoltage"], bins = nbins, range = (0, 2), histtype = 'step', label = 'Max voltage %s'%file.keys()[det])
        plt.title("%s - Run 502" %file.keys()[det])
        plt.xlabel('Max voltage')
        plt.xlim((None, None))
        plt.ylabel('Tiggers/%.2fmV'%(2/nbins * 1000))
        plt.grid()
        plt.semilogy()
        plt.show()



def getHist(array_df):
    Mean_ACT23 = (array_df[4][0] + array_df[5][0] +array_df[6][0] + array_df[7][0])/4
    Mean_ACT1 = (array_df[3][0] + array_df[2][0] )/2
    tof = array_df[0]['TOF']
    leadGlass = array_df[18][0]
    return Mean_ACT23, Mean_ACT1, tof, leadGlass


def sum2Gaussians(x, a1, m1, s1, a2, m2, s2, c):
    return abs(a1) * np.exp(-(m1 - x) **2/ (2 * s1 **2) ) + abs(a2) * np.exp(-(m2 - x) **2/ (2 * s2 **2)) + abs(c)
def oneGaussian(x, a1, m1, s1):
    return a1 * np.exp(-(m1 - x) **2/ (2 * s1 **2) )


def distToLinearCut(x1, y1, a, b):
    #simple trigonometry
    yA = a*x1 + b
    xA = (y1-b)/a
    theta = np.arctan((yA-y1)/(x1-xA))
    return np.sin(theta) * (x1-xA)

def singlePE(argv):
    root_filenames = argv[1]
    file = ur.open(root_filenames)
    nPeaksInToF = nPeakInToF(file, 1)
    print(sum(nPeaksInToF))
    detectors_pd_el = []
    detectors_pd_NONel = []
    detectors_pd_el_wind = []
    detectors_pd_NONel_wind = []
    detectors_pd = []
    detectors_pd_2D_tofLG_selection = []
    detectors_pd_2D_tofLG_selection_e = []
    detectors_pd_2D_ACTs_selection = []
    detectors_pd_2D_ACTs_selection_e = []
    detectors_pd_all = []
    # nPeaksInToF = noHoleHit(file, nPeaksInToF)
    print(sum(nPeaksInToF))
    for detector in range(len(file.keys()[:-1])):
        key = file.keys()[detector]

        df = file[key].arrays(library="pd")
        df = df[nPeaksInToF]
        df = df.reset_index()
        detectors_pd.append(df)

    detectors_pd_safe = detectors_pd.copy()
    detectors_pd_wind = detectors_pd.copy()

    for pmt in range(len(file.keys()[:-1])):
        # if pmt>=8 and pmt < 16:
            #all of the TOF distributions
        hit_times = pd.DataFrame(detectors_pd[pmt]['SignalTime'].values.tolist())
        intCharge = pd.DataFrame(detectors_pd_safe[pmt]['IntCharge'].values.tolist())

        detectors_pd_safe[pmt] = pd.concat([detectors_pd_safe[pmt],intCharge], axis=1)
        detectors_pd[pmt] = pd.concat([detectors_pd[pmt],hit_times], axis=1)

        intCharge = pd.DataFrame(detectors_pd_safe[pmt]['WindowIntPE'].values.tolist())
        detectors_pd_wind[pmt] = pd.concat([detectors_pd_wind[pmt],intCharge], axis=1)

    mean_first_hit_TOF1 = (detectors_pd[12][0] + detectors_pd[13][0] + detectors_pd[14][0] + detectors_pd[15][0])/4
    mean_first_hit_TOF0 = (detectors_pd[8][0] + detectors_pd[9][0] + detectors_pd[10][0] + detectors_pd[11][0])/4

    tof = np.array(np.array(mean_first_hit_TOF1)-np.array(mean_first_hit_TOF0))
    print(tof)


    limCutTOF = 12.5

    #2d cut in LG and TOF
    linearCutA = 0.07286
    linearCutB = -0.704

    #2d cut in ACT1 ACT2
    ACTlinearA = -0.9952#-0.9243 #-0.8436#-1.983 #-1.431#-1.572 #-0.9413 #-1.572 #-1.983 #-1.9375
    ACTlinearB = 7.71#6.09 #4.806 #4.504 #10.00 #6.89 #8.3 #5.59 #6 #10.0 #19

    Pbupper = 0.18

    ACTupper = 11
    ACTlower = 3 #1

    plot_cuts = True

    ACT_selection = (detectors_pd_wind[4][0]+detectors_pd_wind[5][0] + detectors_pd_wind[6][0] + detectors_pd_wind[7][0])/4 > ACTlower


    ACT_selection = ACT_selection & ((detectors_pd_wind[4][0] + detectors_pd_wind[5][0] + detectors_pd_wind[6][0] + detectors_pd_wind[7][0])/4 < ACTlinearA * (detectors_pd_wind[2][0] + detectors_pd_wind[3][0])/2 + ACTlinearB)

    # ACT_selection = ACT_selection & ((detectors_pd_wind[4][0]+detectors_pd_wind[5][0] + detectors_pd_wind[6][0] + detectors_pd_wind[7][0])/4 < ACTupper)


    ACT_selection_e = np.where(ACT_selection == True, False, True)

    a = pd.Series(data = tof, name = 'TOF')
    TOF_selection = a>limCutTOF

    PbTOF_selection = detectors_pd_safe[18][0]<(linearCutA * a + linearCutB)
    PbTOF_selection = PbTOF_selection & (detectors_pd_safe[18][0]< Pbupper)

    PbTOF_selection_e = np.where(PbTOF_selection == True, False, True)



    for detector in range(len(file.keys()[:-1])):

        detectors_pd_safe[detector] = pd.concat([detectors_pd_safe[detector],a], axis=1)
        detectors_pd_wind[detector] = pd.concat([detectors_pd_wind[detector],a], axis=1)

        detectors_pd_el.append(detectors_pd_safe[detector][detectors_pd_safe[detector]['TOF']<limCutTOF])

        detectors_pd_NONel.append(detectors_pd_safe[detector][detectors_pd_safe[detector]['TOF']>limCutTOF])

        detectors_pd_el_wind.append(detectors_pd_wind[detector][detectors_pd_wind[detector]['TOF']<limCutTOF])





        #have peak detect for TOF and lead glass but windPE for ACTs
        if detector <= 7:
            detectors_pd_2D_tofLG_selection.append(detectors_pd_wind[detector][PbTOF_selection])
            detectors_pd_2D_tofLG_selection_e.append(detectors_pd_wind[detector][PbTOF_selection_e])

            detectors_pd_2D_ACTs_selection.append(detectors_pd_wind[detector][ACT_selection])
            detectors_pd_2D_ACTs_selection_e.append(detectors_pd_wind[detector][ACT_selection_e])

            detectors_pd_NONel_wind.append(detectors_pd_wind[detector][TOF_selection])
            detectors_pd_all.append(detectors_pd_wind[detector])

        else:
            detectors_pd_2D_tofLG_selection.append(detectors_pd_safe[detector][PbTOF_selection])
            detectors_pd_2D_tofLG_selection_e.append(detectors_pd_safe[detector][PbTOF_selection_e])

            detectors_pd_NONel_wind.append(detectors_pd_safe[detector][TOF_selection])
            detectors_pd_2D_ACTs_selection.append(detectors_pd_safe[detector][ACT_selection] )
            detectors_pd_2D_ACTs_selection_e.append(detectors_pd_safe[detector][ACT_selection_e])
            detectors_pd_all.append(detectors_pd_safe[detector])
    #
    # for det in range(0,8):
    #     nbins = 200
    #     maxRange = 50
    #     # print(len(detectors_pd_el_wind[det][0]))
    #     plt.hist(detectors_pd_NONel[det][0], bins = nbins, range = (0, maxRange), histtype = 'step', label = 'Integrated pe for non-e - %s'%file.keys()[det])
    #     plt.title("No Hits Hodoscope 1 particle - Run 502")
    #     plt.xlabel('Peak Integrated nb of PE')
    #     plt.xlim((None, None))
    #     plt.ylabel('Tiggers/%.2fpe'%(maxRange/nbins))
    #     plt.grid()
    #     plt.legend()
    #     # plt.semilogy()
    # plt.show()

    # plt.hist2d(detectors_pd_NONel_wind[4][0], detectors_pd_NONel_wind[5][0],  bins = (40, 40), range = ((0, 10), (0,10)), norm = colors.LogNorm())
    # plt.colorbar(label="Number of triggers")
    # plt.xlabel("Number of PE window integrated %s" % file.keys()[4])
    # plt.ylabel("Number of PE window integrated %s" % file.keys()[5])
    # plt.title("PMT left and right comparision - Window Integrated PE")
    # plt.show()

    ACT23_all, ACT1_all, TOF_all, LG_all = getHist(detectors_pd_all)
    ACT23_selection1, ACT1_selection1, TOF_selection1, LG_selection1 = getHist(detectors_pd_NONel_wind)
    ACT23_selection2, ACT1_selection2, TOF_selection2, LG_selection2 = getHist(detectors_pd_2D_tofLG_selection)
    ACT23_selection3, ACT1_selection3, TOF_selection3, LG_selection3 = getHist(detectors_pd_2D_ACTs_selection)
    ACT23_selection4, ACT1_selection4, TOF_selection4, LG_selection4 = getHist(detectors_pd_2D_ACTs_selection_e)
    ACT23_selection5, ACT1_selection5, TOF_selection5, LG_selection5 = getHist(detectors_pd_2D_tofLG_selection_e)


    h_ACTs_selection3 = distToLinearCut(ACT1_selection3, ACT23_selection3, ACTlinearA, ACTlinearB)
    h_ACTs_selection4 = distToLinearCut(ACT1_selection4, ACT23_selection4, ACTlinearA, ACTlinearB)

    # h_ACTs_selection3 = distToLinearCut(TOF_selection2, LG_selection2, linearCutA, linearCutB)
    # h_ACTs_selection4 = distToLinearCut(TOF_selection5, LG_selection5, linearCutA, linearCutB)


    plt.hist(h_ACTs_selection3, bins = 100, range = (-20, 10), label = 'ACT cut - non electrons')
    plt.hist(h_ACTs_selection4, bins = 100, range = (-20, 10), label = 'ACT cut - electrons')
    plt.grid()
    plt.xlabel("Perpendicular distance to the cut line")
    plt.ylabel("number of triggers")
    plt.legend()
    plt.show()




    xmin = 0
    xmax = 20

    plt.hist2d(ACT1_all, ACT23_all, bins = (200, 200), range = ((xmin,xmax), (0, 25)), norm = colors.LogNorm())
    if plot_cuts:
        plt.plot([xmin, xmax], [ACTlower, ACTlower], 'k--', label = '2D cut: lower = %.2f, upper = %.2F'%(ACTupper, ACTlower))
        # plt.plot([xmin, xmax], [ACTupper, ACTupper], 'k--', label = '2D cut: lower = %.2f, upper = %.2F'%(ACTupper, ACTlower))
        plt.plot([xmin, xmax], [xmin*ACTlinearA+ACTlinearB, xmax*ACTlinearA+ACTlinearB], 'k--')
    plt.colorbar(label="Number of triggers")
    plt.xlabel("Mean ACT 1 window PE")
    plt.ylabel("Mean ACT23 window PE")
    plt.title("Run 523")
    plt.legend()
    plt.show()

    xmin = 11
    xmax = 14

    plt.hist2d(TOF_all, LG_all, bins = (200, 200), range = ((xmin,xmax), (0, 0.5)), norm = colors.LogNorm())
    if plot_cuts:
        plt.plot([xmin, xmax], [xmin*linearCutA+linearCutB, xmax*linearCutA+linearCutB], 'k--', label = '2D cut: A = %.3f, B = %.3F'%(linearCutA, linearCutB))
        plt.plot([xmin, xmax], [Pbupper, Pbupper], 'k--')
    plt.colorbar(label="Number of triggers")
    plt.xlabel("Time of flight (ns)")
    plt.ylabel("Peak integrated lead glass")
    plt.title("Run 523")
    plt.legend()
    plt.show()

    xmax = 10
    xmin = -30
    nbins = 100


    bin_edges = np.linspace(xmin, xmax, nbins+1)
    a, _ = np.histogram(h_ACTs_selection4, bins = bin_edges)
    b, _ = np.histogram(h_ACTs_selection3, bins = bin_edges)
    x_ref = np.linspace(xmin, xmax, 1000)


    plt.hist(h_ACTs_selection4, bins = nbins, range = (xmin, xmax), label = 'ACT cut - electrons %i triggers' % len(h_ACTs_selection4), histtype = 'step')

    popt, pcov = curve_fit(oneGaussian, bin_edges[:-1] + (xmax-xmin)/(2*nbins), a, p0 = [10000, -7, 2])
    n_e = quad(oneGaussian, xmin, xmax, args = (popt[0], popt[1], popt[2]))[0] *(nbins)/ (xmax-xmin)
    n_e_inMuPi = quad(oneGaussian, 0, xmax, args = (popt[0], popt[1], popt[2]))[0] *(nbins)/ (xmax-xmin)
    n_e_inE = quad(oneGaussian, xmin, 0, args = (popt[0], popt[1], popt[2]))[0] *(nbins)/ (xmax-xmin)

    plt.plot(x_ref, oneGaussian(x_ref, *popt), 'k--', label = 'Electrons: fitted number %.1f'%(n_e))


    plt.hist(h_ACTs_selection3, bins = nbins, range = (xmin, xmax), label = 'ACT cut - non electrons %i triggers'%len(h_ACTs_selection3) , histtype = 'step')

    popt, pcov = curve_fit(sum2Gaussians, bin_edges[:-1] + (xmax-xmin)/(2*nbins), b, p0 = [200, 1.5, 2, 50, 6, 2, 0])
    n_mu_pi = quad(sum2Gaussians, xmin, xmax, args = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6]))[0] *(nbins)/ (xmax-xmin)
    n_mu_pi_inE = quad(sum2Gaussians, xmin, 0, args = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6]))[0] *(nbins)/ (xmax-xmin)
    n_mu_pi_inMuPi = quad(sum2Gaussians, 0, xmax, args = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6]))[0] *(nbins)/ (xmax-xmin)

    n_mu = quad(oneGaussian, xmin, xmax, args = (popt[0], popt[1], popt[2]))[0] *(nbins)/ (xmax-xmin)
    n_pi = quad(oneGaussian, xmin, xmax, args = (popt[3], popt[4], popt[5]))[0] *(nbins)/ (xmax-xmin)

    plt.plot(x_ref, sum2Gaussians(x_ref, *popt), 'r--', label = 'Non-electrons: fitted number %.1f'%(n_mu_pi))

    print("number of muons: %.1f"%n_mu)
    print("number of pion: %.1f"%n_pi)


    # print(n_e, len(h_ACTs_selection4))
    #
    print("Electron Purity: %.3f"%(n_e_inE/(n_mu_pi_inE + n_e_inE)))
    print("Electron Efficiency: %.3f"% (n_e_inE/n_e))
    #
    print("MuonAndPion Purity: %.3f"%(n_mu_pi_inMuPi/(n_e_inMuPi + n_mu_pi_inMuPi)))
    # print("MuonAndPion Purity: %.3f"%(1-(n_e_inMuPi/len(h_ACTs_selection3)))
    print("MuonAndPion Efficiency: %.3f"%(n_mu_pi_inMuPi/n_mu_pi))
    #
    # print(*popt)

    plt.grid()
    plt.xlabel("Perpendicular distance to the cut line")
    plt.ylabel("number of triggers per %.3f"%((xmax-xmin)/(nbins)))
    # plt.semilogy()
    plt.legend()
    plt.show()


    xmin = 11
    xmax = 13.25
    nbins = 25

    f, (a0, a1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})
    # plt.hist(detectors_pd_wind[0]['TOF'], bins = nbins, range = (xmin, xmax), histtype = 'step', label = 'All particles')
    # plt.hist(TOF_selection1, bins = nbins, range = (xmin, xmax), histtype = 'step', label = '1D TOF cut')
    bin_edges = np.linspace(xmin, xmax, nbins+1)
    print(bin_edges)
    a0.hist(TOF_selection2, bins = nbins, range = (xmin, xmax), histtype = 'step', label = '2D TOF-LG cut %i triggers' %len(TOF_selection2))
    a0.hist(TOF_selection3, bins = nbins, range = (xmin, xmax), histtype = 'step', label = '2D ACT23-ACT1 cut %i triggers' %len(TOF_selection3))



    a, _ = np.histogram(TOF_selection2, bins = bin_edges)
    b, _ = np.histogram(TOF_selection3, bins = bin_edges)

    a0.errorbar(bin_edges[:-1] + (xmax-xmin)/(2*nbins), a, yerr = 1/np.sqrt(a), marker = 'o', ls='none')
    a0.errorbar(bin_edges[:-1] + (xmax-xmin)/(2*nbins), b, yerr = 1/np.sqrt(b), marker = 'o', ls='none')

    print(a)

    a1.errorbar(bin_edges[:-1]+ (xmax-xmin)/(2 * nbins), (b-a)/a, yerr= np.sqrt(1/a + 1/b)/a , marker = 'x', color = 'r', label = 'ACTnonE - TOFLGnonE', ls='none')
    a0.set_xlabel('time of flight (ns)')
    a0.set_ylabel('Number of tiggers per %.3fns'%((xmax-xmin)/nbins))
    a0.set_title("TOF comparison")
    a0.legend()
    a0.grid()
    a1.set_xlabel('time of flight (ns)')
    a1.set_ylabel('Frac. Diff. in nb tiggers')
    a1.legend()
    a1.grid()
    plt.show()




    xmin = 0
    xmax = 0.2
    nbins = 50
    #
    # # plt.hist(LG_selection1, bins = nbins, range = (xmin, xmax), histtype = 'step', label = '1D TOF cut')
    # plt.hist(LG_selection2, bins = nbins, range = (xmin, xmax), histtype = 'step', label = '2D TOF-LG cut')
    # plt.hist(LG_selection3, bins = nbins, range = (xmin, xmax), histtype = 'step', label = '2D ACT23-ACT1 cut')
    # plt.xlabel('LeadGlass peak Int')
    # plt.ylabel('Number of tiggers per %.3fns'%((xmax-xmin)/nbins))
    # plt.grid()
    # plt.title("Lead Glass comparison")
    # plt.legend()
    # plt.show()


    f, (a0, a1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})
    # plt.hist(detectors_pd_wind[0]['TOF'], bins = nbins, range = (xmin, xmax), histtype = 'step', label = 'All particles')
    # plt.hist(TOF_selection1, bins = nbins, range = (xmin, xmax), histtype = 'step', label = '1D TOF cut')
    bin_edges = np.linspace(xmin, xmax, nbins)
    print(bin_edges)
    a0.hist(LG_selection2, bins = nbins, range = (xmin, xmax), histtype = 'step', label = '2D TOF-LG cut %i'%len(LG_selection2))
    a0.hist(LG_selection3, bins = nbins, range = (xmin, xmax), histtype = 'step', label = '2D ACT23-ACT1 cut %i triggers'%len(LG_selection3))

    a, _ = np.histogram(LG_selection2, bins = bin_edges)
    b, _ = np.histogram(LG_selection3, bins = bin_edges)

    print(a)

    a1.errorbar(bin_edges[:-1], (b-a)/a, yerr= np.sqrt(1/a + 1/b)/a , marker = 'x', color = 'r', label = 'ACTnonE - TOFLGnonE', ls='none')
    a0.set_xlabel('LeadGlass peak Int')
    a0.set_ylabel('Number of tiggers per %.3f'%((xmax-xmin)/nbins))
    a0.set_title("LeadGlass comparison")
    a0.legend()
    a0.grid()
    a1.set_xlabel('LeadGlass peak Int')
    a1.set_ylabel('Frac. Diff. in nb tiggers')
    a1.legend()
    a1.grid()
    plt.show()


    xmin = 0
    xmax = 6
    nbins = 100

    f, (a0, a1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1, 1]})
    # plt.hist(detectors_pd_wind[0]['TOF'], bins = nbins, range = (xmin, xmax), histtype = 'step', label = 'All particles')
    # plt.hist(TOF_selection1, bins = nbins, range = (xmin, xmax), histtype = 'step', label = '1D TOF cut')
    bin_edges = np.linspace(xmin, xmax, nbins)
    print(bin_edges)

    a0.hist(ACT23_selection2, bins = nbins, range = (xmin, xmax), histtype = 'step', label = 'Non-electrons-like %i triggers (LG-TOF selection)'%len(ACT23_selection2))
    a0.hist(ACT23_selection5, bins = nbins, range = (xmin, xmax), histtype = 'step', label = 'Electrons-like %i triggers (LG-TOF selection)'%len(ACT23_selection5))




    a0.semilogy()
    # a1.semilogy()



    print(a)
    #
    # a1.errorbar(bin_edges[:-1], (b-a)/a, yerr= np.sqrt(1/a + 1/b)/a , marker = 'x', color = 'r', label = 'ACTnonE - TOFLGnonE', ls='none')

    a1.hist(ACT23_selection3, bins = nbins, range = (xmin, xmax), histtype = 'step', label = 'Non-electrons-like %i triggers (ACTs selection)'%len(ACT23_selection3))
    a1.hist(ACT23_selection4, bins = nbins, range = (xmin, xmax), histtype = 'step', label = 'Electrons-like %i triggers (ACTs selection)'%len(ACT23_selection4))

    a, _ = np.histogram(ACT23_selection3, bins = bin_edges)
    b, _ = np.histogram(ACT23_selection4, bins = bin_edges)



    x_ref = np.linspace(xmin, xmax, 1000)
    #here the fit to the ACT23 information
    popt, pcov = curve_fit(sum2Gaussians, bin_edges[:-1] + (xmax-xmin)/(2*nbins), a, p0 = [332.15273834,   0.91885663,   0.51667649, 35, 4, 6, 0])
    # popt = [343, 1.15, 0.1, 376, 1.45, 0.4, 75]
    # a1.plot(x_ref, sum2Gaussians(x_ref, *popt), 'k-')

    #check the purity
    window_low_tot = 0
    window_high_tot = 25

    n_muon_tot = quad(oneGaussian, window_low_tot, window_high_tot, args = (popt[3], popt[4], popt[5]))[0]*(nbins)/ (xmax-xmin)
    n_pion_tot = quad(oneGaussian, window_low_tot, window_high_tot, args = ( popt[0], popt[1], popt[2]))[0]*(nbins)/ (xmax-xmin)

    # print(n_muon_tot, n_pion_tot)

    print("Total number of muons: %.1f" % n_muon_tot)
    print("Total number of pions: %.1f" % n_pion_tot)



    a1.plot(x_ref, oneGaussian(x_ref, *popt[0:3]), 'g--', label = "Total fitted number of muons: %.1f" % n_muon_tot)
    a1.plot(x_ref, oneGaussian(x_ref, *popt[3:6]), 'r--', label = "Total fitted number of pions: %.1f" % n_pion_tot)

    print("Best fit parameters:", popt)


    a0.set_xlabel('ACT23 mean Number of PE')
    a0.set_ylabel('Number of tiggers per %.3f'%((xmax-xmin)/nbins))
    a0.set_title("Light yield")
    a0.legend()
    a0.grid()

    a1.set_xlabel('ACT23 mean Number of PE')
    a1.set_ylabel('Number of tiggers per %.3f'%((xmax-xmin)/nbins))
    a1.set_title("Light yield")
    a1.legend()
    a1.grid()
    # a1.set_xlabel('LeadGlass peak Int')
    # a1.set_ylabel('Frac. Diff. in nb tiggers')
    # a1.legend()
    # a1.grid()
    plt.show()

    #Plot the purities
    window_low_mu = 11
    window_high_mu = 25

    muon_efficiency =[]
    muon_purity =[]
    window_start =[]

    pion_efficiency=[]
    pion_purity =[]
    window_end =[]


    for window_low_mu in np.arange(window_low_tot+1, window_high_tot-1, 1):
        n_muon_windowM = quad(oneGaussian, window_low_mu, window_high_mu, args = (popt[3], popt[4], popt[5]))[0]*(nbins)/ (xmax-xmin)
        n_pion_windowM = quad(oneGaussian, window_low_mu, window_high_mu, args = (popt[0], popt[1], popt[2]))[0]*(nbins)/ (xmax-xmin)

        muon_efficiency.append(n_muon_windowM/n_muon_tot)
        muon_purity.append((n_muon_windowM)/(n_pion_windowM+n_muon_windowM))

        if window_low_mu == 9:
            print("Efficiciency mu %.3f"%(n_muon_windowM/n_muon_tot))
            print("purity mu %.3f"% ((n_muon_windowM)/(n_pion_windowM+n_muon_windowM)))
        window_start.append(window_low_mu)

    window_low_pi = 0
    window_high_pi = 11


    for window_high_pi in np.arange(window_low_tot+1, window_high_tot-1, 1):
        n_muon_windowP = quad(oneGaussian, window_low_pi, window_high_pi, args = (popt[3], popt[4], popt[5]))[0]*(nbins)/ (xmax-xmin)
        n_pion_windowP = quad(oneGaussian, window_low_pi, window_high_pi, args = ( popt[0], popt[1], popt[2]))[0]*(nbins)/ (xmax-xmin)

        if window_high_pi == 9:
            print("Efficiciency pi %.3f"% (n_pion_windowP/n_pion_tot))
            print("purity pi %.3f"% ((n_pion_windowP)/(n_pion_windowP+n_muon_windowP)))

        pion_efficiency.append(n_pion_windowP/n_pion_tot)
        pion_purity.append((n_pion_windowP)/(n_pion_windowP+n_muon_windowP))
        window_end.append(window_high_pi)




    fig, ax1 = plt.subplots()
    ax1.set_title("Muon Purity and Efficiency - Run 490")
    ax1.set_xlabel('Lower ACT23 cut (PE)')
    ax1.set_ylabel('Purity', color = 'green')
    ax1.plot(window_start, muon_purity, color = 'green',  marker = 'o')
    ax1.tick_params(axis ='y', labelcolor = 'green')

    ax2 = ax1.twinx()
    ax2.set_ylabel('Efficiency', color = 'red')
    ax2.plot(window_start, muon_efficiency, color = 'red', marker = 'o')
    ax2.tick_params(axis ='y', labelcolor = 'red')
    plt.show()

    fig, ax1 = plt.subplots()
    ax1.set_title("Pion Purity and Efficiency - Run 490")
    ax1.set_xlabel('Upper ACT23 cut (PE)')
    ax1.set_ylabel('Purity', color = 'green')
    ax1.plot(window_end, pion_purity, color = 'green',  marker = 'o')
    ax1.tick_params(axis ='y', labelcolor = 'green')

    ax2 = ax1.twinx()
    ax2.set_ylabel('Efficiency', color = 'red')
    ax2.plot(window_end, pion_efficiency, color = 'red', marker = 'o')
    ax2.tick_params(axis ='y', labelcolor = 'red')
    plt.show()







    xmin = 0
    xmax = 0.5
    nbins = 50

    plt.hist(LG_all, bins = nbins, range = (xmin, xmax), histtype = 'step', label = 'All')
    # plt.hist(LG_selection2, bins = nbins, range = (xmin, xmax), histtype = 'step', label = '2D TOF-LG cut')
    plt.xlabel('LeadGlass peak Int')
    plt.ylabel('Number of tiggers per %.3fns'%((xmax-xmin)/nbins))
    plt.grid()
    plt.title("Lead Glass all")
    plt.legend()
    plt.show()

    xmin = 11
    xmax = 14

    plt.hist2d(TOF_selection3, ACT23_selection3, bins = (200, 200), range = ((xmin,xmax), (0, 20)), norm = colors.LogNorm())
    plt.colorbar(label="Number of triggers")
    plt.xlabel("Time of flight (ns)")
    plt.ylabel("Mean ACT23 window PE")
    plt.show()


    raise end








if __name__ == "__main__":
    # execute only if run as a script"
    # lookingAtEventsAndPlotting(sys.argv)
    # multipeaksQualityCheck(sys.argv)
    # meanDistancePeaks(sys.argv)
    # cutTOF(sys.argv)
    singlePE(sys.argv)
    # peakMatching(sys.argv)
    # print(root_filenames)




