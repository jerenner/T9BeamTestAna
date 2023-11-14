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
        for hit in range(int(max(detectors_pd[pmt]['nPeaks']))):
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



def lookingAtEventsAndPlotting(argv):
    root_filenames = argv[1]
    file = ur.open(root_filenames)
    detectors_pd = []
    detectors_pd_charge = []

    print(file.keys())
    print(file['TOF03;1'].keys())

    #it is useful to look at only the case where we have one peak in all of the TOF PMTs - a bit disgusting
    df = file['TOF00'].arrays(library="pd")
    onePeakToF = (df['nPeaks']==1)
    df = file['TOF01'].arrays(library="pd")
    onePeakToF = (df['nPeaks']==1) & onePeakToF
    df = file['TOF02'].arrays(library="pd")
    onePeakToF = (df['nPeaks']==1) & onePeakToF
    df = file['TOF03'].arrays(library="pd")
    onePeakToF = (df['nPeaks']==1) & onePeakToF
    df = file['TOF10'].arrays(library="pd")
    onePeakToF = (df['nPeaks']==1) & onePeakToF
    df = file['TOF11'].arrays(library="pd")
    onePeakToF = (df['nPeaks']==1) & onePeakToF
    df = file['TOF12'].arrays(library="pd")
    onePeakToF = (df['nPeaks']==1) & onePeakToF
    df = file['TOF13'].arrays(library="pd")
    onePeakToF = (df['nPeaks']==1) & onePeakToF
    #
    df = file['%s'%file.keys()[6]].arrays(library="pd")
    for i in range(len(df['nPeaks'])):
        # print(i/len(df['nPeaks']),  np.all(df['IntCharge'][i]>=0.02) & onePeakToF[i])
        onePeakToF[i] = np.any(df['IntCharge'][i]>=0.015) & onePeakToF[i]


    # print(onePeakToF, df['nPeaks'])

#read in each of the branches
    for key in file.keys()[:-1]:
        df = file[key].arrays(library="pd")
        #it is useful to look at only the case where we have one peak in all of the TOF
        detectors_pd.append(df[onePeakToF])
        detectors_pd_charge.append(df[onePeakToF])
        print(df[onePeakToF])

    all_first_hit_times = []
    all_second_hit_times = []
    all_hit_times = []
    all_hit_charge = []

    #make an array of [DataframePMT1, DataFramePMT2, ...] where DataFramePMT2 has columns 0,1, 2, 3... corresponding to the hit time of each hit, 0th 1st...
    for pmt in range(19):
        print(pmt)
        hit_times = pd.DataFrame(detectors_pd[pmt]['SignalTime'].values.tolist())
        hit_charge = pd.DataFrame(detectors_pd[pmt]['IntCharge'].values.tolist())
        detectors_pd[pmt] = pd.concat([detectors_pd[pmt],hit_times], axis=1)
        detectors_pd_charge[pmt] = pd.concat([detectors_pd_charge[pmt],hit_charge], axis=1)
        all_first_hit_times.append(hit_times[0])
        # all_second_hit_times.append(hit_times[1])

    #Calculate the mean first hit time in each of the TOF detectors - no peak matching at this point
    mean_first_hit_TOF0 = (detectors_pd[11][0] + detectors_pd[10][0] + detectors_pd[9][0] + detectors_pd[8][0])/4
    mean_first_hit_TOF1 = (detectors_pd[12][0] + detectors_pd[13][0] + detectors_pd[14][0] + detectors_pd[15][0])/4




    #here, reading in the dataframes to fill an array of all of the hits for each PMT (easier than to access the dataframe (a bit convoluted, I agree,
    for pmt in range(19):
        all_hit_times.append([])
        all_hit_charge.append([])
        print("Pmt:", pmt, " max number of peaks is:", max(detectors_pd[pmt]['nPeaks']))
        if pmt == 6:
            for hit in range(int(max(detectors_pd[pmt]['nPeaks']))): # int(max(detectors_pd[pmt]['nPeaks'])
                all_hit_times[pmt].append(detectors_pd[pmt][hit])
                all_hit_charge[pmt].append(detectors_pd_charge[pmt][hit])


#to look at different detectors and plot them more easily

    i_min = 6 #the range of detectors of interest
    i_max = 7
    list_markers = ['o','v', '^', 's', '<', '>', 'X', 'd', 'H', 'h', 'D', 'P', 'p','*', 'o','v', '^', 's', '>', '<']
    PMTref = 11
    # print(all_hit_times)

    def gaussian(x, mean, amplitude, standard_deviation):
        return amplitude * np.exp( - (x - mean)**2 / (2*standard_deviation ** 2))

    def linear(x, slope, intercept):
        return x * slope + intercept

    #now we will compare the histogram of the difference between each detector and the reference: mean TOF
    for i in range(i_min, i_max):
        plt.figure(1)

        hist = all_hit_times[i][0]
        hist_charge = all_hit_charge[i][0]
        # for k in range(end):
        #     if np.isnan(hist.loc[k]):
        #         print('yes', k)
        #         hist = hist.drop(index = k)
        #         hist_charge = hist_charge.drop(index = k)
        #         mean_first_hit_TOF1 = mean_first_hit_TOF1.drop(index = k)
        #         counter +=1

        # .dropna(how='any') #here read the first entry of the given PMT
        # hist = hist.fillna(value=100) #hist.mean()
        # mean_first_hit_TOF1 = mean_first_hit_TOF1.dropna(how='any')
        # hist = list(hist)
        # condition = (
        # hist = hist[np.isnan(hist)==False]
        # mean_first_hit_TOF1 = mean_first_hit_TOF1[not(np.isnan(hist))]

        print('i:', i, hist)
        #Making a histogram of the distrubution of the delay between each detector and the TOF1 mean hit time and fitting it
        bin_heights, bin_borders, _ = plt.hist(hist - mean_first_hit_TOF1,  bins = 100, range = (-100, 100) , alpha = 0.3)
        bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
        popt, _ = curve_fit(gaussian, bin_centers, bin_heights, p0=[1., 0., 1.])

        x_interval_for_fit = np.linspace(bin_borders[0], bin_borders[-1], 10000)
        plt.plot(x_interval_for_fit, gaussian(x_interval_for_fit, *popt), label='%s: mean = %.2fns, std = %.2fns '%(file.keys()[i],popt[0],popt[2]) )

        plt.figure(2)
        plt.hist2d(hist-mean_first_hit_TOF1, hist_charge, bins = [200, 100], range = [[-100, 100], [0, 0.5]], label = file.keys()[i], norm = colors.LogNorm())

        plt.figure(4)
        plt.plot(hist-mean_first_hit_TOF1, 'x', label = file.keys()[i])
        x = np.linspace(0, len(hist), len(hist))
        # fit_params = np.polyfit(x[:100], hist[:100]-mean_first_hit_TOF1[:100], 2)
        # p = np.poly1d(fit_params)
        # print(hist[10:20], mean_first_hit_TOF1[10:20], hist, mean_first_hit_TOF1)
        # plt.plot(x, p(x), 'r--', label = 'linear fit: %.3fx + %.3e'%(fit_params[0], fit_params[1]))
        # popt, _ = curve_fit(linear, x[:end-counter], hist[:end-counter]-mean_first_hit_TOF1[:end-counter], p0=[0., -40.])
        # plt.plot(x, linear(x, *popt), label='%s: slope = %.2fns, intercept = %.2fns '%(file.keys()[i],popt[0], popt[1] ))
        # H, xedges, yedges = np.histogram2d(hist_charge, hist-mean_first_hit_TOF1)
        # H = H.T
        # # fig = plt.figure(figsize=(7, 3))
        # fig = plt.figure(figsize=(7, 3))
        # ax = fig.add_subplot(131, title='imshow: square bins')
        # plt.imshow(H, interpolation='nearest', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])



    plt.figure(1)
    plt.grid()
    plt.legend()
    # plt.semilogy()
    plt.xlabel('EARLIEST HIT MeanTOF1hitTime - EARLIEST PMT hit time (ns)')
    plt.ylabel('Triggers/2ns')
    # plt.show()
    plt.figure(2)
    plt.grid()
    plt.title("Hit delay compared to mean tof1 vs hit charge %s\n One Hit in all TOFs"%file.keys()[i])
    plt.colorbar(label="Number of triggers/(2ns,0.005V)")
    plt.xlabel('Other PMT EARLIEST hitTime - Mean TOF1 EARLIEST(ns)')
    plt.ylabel('Other PMT EARLIEST integrated hit charge (V)')

    plt.figure(4)

    plt.grid()
    plt.ylabel('Other PMT EARLIEST hitTime - Mean TOF1 EARLIEST(ns)')
    plt.xlabel('Event number')
    plt.legend()
    # plt.legend()
    plt.show()


    #Now look at all of the hit times for a fw events, see if we can match the signal seen in the four TOF modules, useful to find events of interest, e.g when more than one particle crosses the TOF within the window.
    i = 0
    plt.figure(3)
    for event in range(25):
        for pmt in [8, 9, 10, 11]:
            if event == 0:
                print(i, pmt)
                plt.scatter(event*20+3, detectors_pd[pmt][0][event*20+3], marker = '%s'%list_markers[pmt], c = 0, vmin=0., vmax=6.,cmap = 'Set1', label = file.keys()[pmt])
                i = i+1

            for hit in range(int(max(detectors_pd[pmt]['nPeaks']))):
                plt.scatter(event*20+3, detectors_pd[pmt][hit][event*20+3], marker = '%s'%list_markers[pmt], c = hit, vmin=0., vmax=6.,cmap = 'Set1')

    plt.grid()
    plt.legend()
    plt.xlabel('Event')
    plt.colorbar(label="Hit Index")
    plt.ylabel('Hit time')
    # plt.show()


    #For a given event look at all of the hit times
    i = 0
    x_ticks_labels = []
    ax = plt.figure(4)
    hit_times = []
    event_of_interest = 181
    for event in [event_of_interest]:
        print('TOFTriggerTime: ', detectors_pd[9]['triggerTime'][event])
        print('ACTTriggerTime: ', detectors_pd[0]['triggerTime'][event])
        # print('ACTTimeStamp: ', detectors_pd[0]['triggerTime'][event])
        print('ACTTriggerTime - TOFTriggerTime: ', detectors_pd[0]['triggerTime'][event]-detectors_pd[9]['triggerTime'][event])
        print('ACTtimeStamp - TOFtimeStamp: ', detectors_pd[0]['timeStamp'][event]-detectors_pd[9]['timeStamp'][event])

        print('ACT3L hit 1 Time', detectors_pd[6][0][event])
        print('ACT3R hit 1 Time', detectors_pd[7][0][event])

        print('Mean TOF hit 1 Time', (detectors_pd[12][0][event]+detectors_pd[13][0][event]+detectors_pd[14][0][event]+detectors_pd[15][0][event])/4)

        # print('ACT3L hit 2 Time', detectors_pd[6][1][event])
        # print('ACT3R hit 2 Time', detectors_pd[7][1][event])
        #
        # print('Mean TOF hit 2 Time', (detectors_pd[12][1][event]+detectors_pd[13][1][event]+detectors_pd[14][1][event]+detectors_pd[15][1][event])/4)

        for pmt in range(19):
            hit_times.append([])

            if pmt != 16 and pmt != 17:
                x_ticks_labels.append(file.keys()[pmt])

                if pmt == 0:

                    # print(i, pmt)
                    plt.scatter(file.keys()[pmt], detectors_pd[pmt][0][event], marker = '%s'%list_markers[event-event_of_interest], c = 0, vmin=0., vmax=6.,cmap = 'Set1', label = 'Event %i'%(event))
                    i = i+1

                for hit in range(int(max(detectors_pd[pmt]['nPeaks']))):

                    if detectors_pd[pmt][hit][event] > 0:
                        hit_times[pmt].append(detectors_pd[pmt][hit][event])

                    if hit == 0 :
                        print(detectors_pd[pmt][0][event], "," , detectors_pd[pmt][1][event], ",",detectors_pd[pmt][2][event])



                    plt.scatter(file.keys()[pmt], detectors_pd[pmt][hit][event], marker = '%s'%list_markers[event-event_of_interest], c = hit, vmin=0., vmax=6.,cmap = 'Set1')

    # raise end

    print("Event: ", event_of_interest, ' hit times: ', hit_times)
    # raise end
    plt.grid()
    plt.legend()
    plt.xlabel('PMT')
    plt.colorbar(label="Hit Index")
    plt.ylabel('Hit time')
    plt.show()






    raise end

    nsperpoint = 2
    print(mean_first_hit_TOF1[0])
    mean_first_hit_TOF1 = mean_first_hit_TOF1/(4*nsperpoint)
    print(mean_first_hit_TOF1[0])
    mean_first_hit_TOF0 = mean_first_hit_TOF0/(4*nsperpoint)
    # print(mean_first_hit_TOF0)
    plt.xlabel('Mean time of the first SignalTime in TOF1(ns)')
    plt.ylabel('Events')
    plt.grid()
    plt.hist(mean_first_hit_TOF1, bins = 100)
    plt.show()


    plt.scatter(mean_first_hit_TOF0[:200], mean_first_hit_TOF1[:200])
    plt.xlabel('Mean time of the first hit in the TOF0')
    plt.ylabel('Mean time of the first hit in the TOF1')
    plt.grid()
    plt.plot([min(mean_first_hit_TOF0[:200]), max(mean_first_hit_TOF0[:200])], [min(mean_first_hit_TOF0[:200])+11.64, max(mean_first_hit_TOF0[:200])+11.64], 'k--', label = 'particle travelling at c')
    plt.legend()
    plt.show()


    plt.hist(mean_first_hit_TOF1-mean_first_hit_TOF0, bins = 100)
    plt.grid()
    plt.xlabel("TOF1-TOF0")
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
    lookingAtEventsAndPlotting(sys.argv)
    # peakMatching(sys.argv)
    # print(root_filenames)




