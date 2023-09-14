import uproot as ur

import matplotlib.pyplot as plt

import numpy as np

import pandas as pd
import sys

import matplotlib.colors as colors

from scipy.optimize import curve_fit

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

    df = file['TOF00'].arrays(library="pd")
    onePeakToF = (df['nPeaks']!=100)
    df = file['TOF01'].arrays(library="pd")
    onePeakToF = (df['nPeaks']==1) & onePeakToF
    df = file['TOF02'].arrays(library="pd")
    onePeakToF = (df['nPeaks']==1) & onePeakToF
    df = file['TOF03'].arrays(library="pd")
    onePeakToF = (df['nPeaks']==1) & onePeakToF
    df = file['TOF10'].arrays(library="pd")
    onePeakToF = (df['nPeaks']==1) & onePeakToF
    print(df['IntCharge'], np.where(df['IntCharge'][:][0]<=0.015, False, onePeakToF))
    #need to get rid of the accidentals in TOF10 -> pretty unlikely in this case of only one
    onePeakToF = np.where(df['IntCharge'][:][0]<=0.015, False, onePeakToF)
    df = file['TOF11'].arrays(library="pd")
    onePeakToF = (df['nPeaks']==1) & onePeakToF
    df = file['TOF12'].arrays(library="pd")
    onePeakToF = (df['nPeaks']==1) & onePeakToF
    df = file['TOF13'].arrays(library="pd")
    onePeakToF = (df['nPeaks']==1) & onePeakToF


    # df = file['ACT0L'].arrays(library="pd")
    # onePeakToF = (df['nPeaks']==1) & onePeakToF





#read in each of the branches
    for key in file.keys()[:-1]:
        df = file[key].arrays(library="pd")
        #it is useful to look at only the case where we have one peak in all of the TOF
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
        for hit in range(int(min(max(detectors_pd[pmt]['nPeaks']), max(detectors_pd[12]['nPeaks'])))): # int(max(detectors_pd[pmt]['nPeaks'])
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

        hist = all_hit_times[i][0]
        #TOF1 - TOF0
        particleTOF =  (all_hit_times[12][0]+all_hit_times[13][0]+all_hit_times[14][0]+all_hit_times[15][0] - (all_hit_times[8][0] + all_hit_times[9][0] + all_hit_times[10][0] + all_hit_times[11][0]))/4

        print(particleTOF, all_hit_times[8][0])


        hist_charge = all_hit_charge[i][0]
        leadGlass_charge = all_hit_charge[18][0]
        hist_window_charge = all_hit_window_charge[i][0]

        #
        # print(detectors_pd_window_charge[i]["WindowMinusPeak"], 'Window Minus peak Int charge ')
        # print(detectors_pd_window_charge[i][detectors_pd_window_charge[i]["WindowMinusPeak"])
        #



        hist_window_charge = all_hit_window_charge[i][0]


        detectors_pd_window_charge[i]["WindowMinusPeak"] =  hist_window_charge -  hist_charge
        print(detectors_pd_window_charge[i][detectors_pd_window_charge[i]["WindowMinusPeak"] > 0.025])

        counter =0

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
        plt.title("One hit in TOFs  %s"%file.keys()[i])
        plt.colorbar(label="Number of triggers/(0.5ns*0.002V)")
        plt.ylabel('Window integrated - Peak integrated hit charge (V)')
        plt.xlabel('HitTime - TOF1 Mean hit time(ns)')
        plt.show()


        plt.figure(5)
        plt.hist2d(particleTOF, hist_window_charge - hist_charge,  bins = [200, 200], range = [[10,13],[-0.1, 0.1]],  label = file.keys()[i], norm = colors.LogNorm())
        plt.grid()
        plt.title("One hit in TOFs  %s"%file.keys()[i])
        plt.colorbar(label="Number of triggers/(0.1ns*0.002V)")
        plt.ylabel('Window integrated - Peak integrated hit charge (V)')
        plt.xlabel('Particle Time of flight(ns)')
        plt.show()


        plt.figure(6)
        plt.subplot(2,1,1)
        plt.hist2d(leadGlass_charge, hist_charge,  bins = [60, 300], range = [[0, 0.6],[-0.04, 0.28*2]],  label = file.keys()[i], norm = colors.LogNorm())
        plt.grid()
        plt.title("One hit in TOFs  %s"%file.keys()[i])
        plt.colorbar(label="Number of triggers/(0.01V*0.001V)")
        plt.ylabel('Peak integrated hit charge (V)')
        plt.xlabel('LeadGlass peak integrated charge')

        plt.subplot(2,1,2)
        plt.hist2d(leadGlass_charge, hist_window_charge,  bins = [60, 300], range = [[0, 0.6],[-0.04, 0.28*2]],  label = file.keys()[i], norm = colors.LogNorm())
        plt.grid()
        # plt.title("One hit in TOFs  %s"%file.keys()[i])
        plt.colorbar(label="Number of triggers/(0.01V*0.001V)")
        plt.ylabel('Window integrated integrated hit charge (V)')
        plt.xlabel('LeadGlass peak integrated charge')

        plt.figure(7)
        plt.subplot(2,1,1)
        plt.hist2d(particleTOF, hist_charge,  bins = [200, 300], range = [[10,13],[-0.04, 0.28*2]],  label = file.keys()[i], norm = colors.LogNorm())
        plt.grid()
        plt.title("One hit in TOFs  %s"%file.keys()[i])
        plt.colorbar(label="Number of triggers/(0.1ns*0.001V)")
        plt.ylabel('Peak integrated hit charge (V)')
        plt.xlabel('Particle Time of flight(ns)')

        plt.subplot(2,1,2)
        plt.hist2d(particleTOF, hist_window_charge,  bins = [200, 300], range = [[10,13],[-0.04, 0.28*2]],  label = file.keys()[i], norm = colors.LogNorm())
        plt.grid()
        plt.title("One hit in TOFs  %s"%file.keys()[i])
        plt.colorbar(label="Number of triggers/(0.1ns*0.001V)")
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

if __name__ == "__main__":
    # execute only if run as a script"
    # lookingAtEventsAndPlotting(sys.argv)
    # multipeaksQualityCheck(sys.argv)
    meanDistancePeaks(sys.argv)
    # peakMatching(sys.argv)
    # print(root_filenames)




