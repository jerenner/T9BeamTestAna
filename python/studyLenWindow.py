#this is a small plotting code to figure out the optimal size of the integration window.
#For each detector we want the integrated charge as a function of the integration window
import uproot as ur

import matplotlib.pyplot as plt

import numpy as np

import pandas as pd
import sys

import matplotlib.colors as colors

from scipy.optimize import curve_fit

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




def studyLenWindow(argv):
    root_filenames = [f for f in argv[1:]]
    window_integrated = [] #to host the integrated window charge
    peak_integrated = [] #useful to have the peak integrated charge too to compare with
    # (eventtually we might want to do that for a subset of events, or actually for each event and for each window size setup?
    #we look at the files ones after the next
    fig = plt.figure(1,figsize=(25, 15))
    fig.tight_layout()

    event = 10000
    event2 = 20000
    event3 = 30000
    lenLongAverage = 5000

    TwoHits = False

    if TwoHits:
        fig.suptitle("Run 523 - study of window size\nTwo particle in each TOF", fontsize=16, fontweight='bold')
    else:
        fig.suptitle("Run 523 - study of window size\nOne particle in each TOF", fontsize=16, fontweight='bold')

    y_array_1 = [[],[],[],[]]
    y_array_2 = [[],[],[],[]]
    y_array_3 = [[],[],[],[]]
    y_array_4 = [[],[],[],[]]
    y_array_5 = [[],[],[],[]]

    if TwoHits:
        y_array_12 = [[],[],[],[]]
        y_array_22 = [[],[],[],[]]
        y_array_32 = [[],[],[],[]]
        y_array_42 = [[],[],[],[]]
        y_array_52 = [[],[],[],[]]


    list_xlabels = []
    #look at TOF00, TOF10, ACT0L, ACT2L
    detectorOfInterest = [8, 12, 5, 7]
    for file_name in root_filenames:

        print(file_name)

        lower_bound = float(file_name.split("_")[1][:-2])
        upper_bound = float(file_name.split("_")[-1][:-7])
        list_xlabels.append(-lower_bound + upper_bound)
        # print(lower_bound
        file = ur.open(file_name)
        if TwoHits:
            nPeaksInToF = nPeakInToF(file, 2) #we only want one hit
        else:
            nPeaksInToF = nPeakInToF(file, 1)

        file_window_sum = []
        file_int_sum = []

        if TwoHits:
            file_window_sum2 = []
            file_int_sum2 = []

        for detector in range(len(file.keys()[:-1])):
            key = file.keys()[detector]
            print(key)
            df = file[key].arrays(library="pd")
            df = df[nPeaksInToF]
            df = df.reset_index() #very important!



            intCharge = pd.DataFrame(df['IntCharge'].values.tolist())
            dfPeak = pd.concat([df, intCharge], axis = 1)

            windowIntCharge = pd.DataFrame(df['WindowIntCharge'].values.tolist())
            dfWindow = pd.concat([df, windowIntCharge], axis = 1)
            # detectors_pd.append(df[onePeakToF])
            if (not TwoHits):
                file_int_sum.append(list(map(sum, df['IntCharge'].values.tolist())))
                file_window_sum.append(list(map(sum, df['WindowIntCharge'].values.tolist())))


            if TwoHits:
                dfWindow[0] = dfWindow[0].fillna(0)
                dfPeak[0] = dfPeak[0].fillna(0)

                dfWindow[1] = dfWindow[1].fillna(0)
                dfPeak[1] = dfPeak[1].fillna(0)
                file_int_sum2.append(list(dfPeak[1]))
                file_window_sum2.append(dfWindow[1])
                 # print(dfPeak[0])
                file_int_sum.append(list(dfPeak[0]))
                file_window_sum.append(list(dfWindow[0]))




        for plot in range(4):
            y_array_4[plot].append(sum(file_window_sum[detectorOfInterest[plot]][:])/len(file_int_sum[detectorOfInterest[plot]]))
            y_array_5[plot].append(sum(file_int_sum[detectorOfInterest[plot]][:])/len(file_int_sum[detectorOfInterest[plot]]))


            if (not TwoHits):
                y_array_1[plot].append(sum(file_window_sum[detectorOfInterest[plot]][event-lenLongAverage:event+lenLongAverage])/(2*lenLongAverage))
                y_array_2[plot].append(sum(file_window_sum[detectorOfInterest[plot]][event2-lenLongAverage:event2+lenLongAverage])/(2*lenLongAverage))
                y_array_3[plot].append(sum(file_window_sum[detectorOfInterest[plot]][event3-lenLongAverage:event3+lenLongAverage])/(2*lenLongAverage))

            if TwoHits:
                y_array_42[plot].append(sum(file_window_sum2[detectorOfInterest[plot]][:])/len(file_int_sum2[detectorOfInterest[plot]]))
                y_array_52[plot].append(sum(file_int_sum2[detectorOfInterest[plot]][:])/len(file_int_sum2[detectorOfInterest[plot]]))




    for plot in range(4):
            plt.subplot(2,2,plot+1)
            plt.title('%s'%(file.keys()[detectorOfInterest[plot]][:-2])) #
            plt.plot(list_xlabels, y_array_4[plot], 'gx-', label = 'Mean windowIntCharge all Peak0 events')
            plt.plot(list_xlabels, y_array_5[plot], 'k-', label = 'Mean Peak0 peakIntCharge')
            if (not TwoHits):
                plt.plot(list_xlabels, y_array_1[plot], 'rx-', label = 'Mean windowIntCharge Peak0 %i evts centre: %i'%(lenLongAverage*2, event))
                plt.plot(list_xlabels, y_array_2[plot], 'x-', color = 'darkgoldenrod', label = 'Mean windowIntCharge%i Peak0 evts centre: %i'%(lenLongAverage*2, event2))
                plt.plot(list_xlabels, y_array_3[plot], 'bx-', label = 'Mean windowIntCharge %i Peak0 evts centre: %i'%(lenLongAverage*2, event3))

            if TwoHits:
                plt.plot(list_xlabels, y_array_42[plot], 'gx--', label = 'Mean windowIntCharge all Peak1 events')
                plt.plot(list_xlabels, y_array_52[plot], 'k--', label = 'Mean Peak1 peakIntCharge')


            plt.legend()
            plt.grid()
            plt.ylabel("Window integrated hit charge")
            plt.xlabel("Length of the integration window")

    plt.subplot(2,2,1)
    #
    # plt.title("Run 490 - study of window size")
    if TwoHits:
        plt.savefig("pdf_results/studyLenWindow_run490_scan_TwoHits.pdf")
    else:
        plt.savefig("pdf_results/studyLenWindow_run490_scan_oneHit.pdf")
    plt.show()





def overlayACTprofiles(argv):
    #this is to overlay the histograms of the ACT signal to get an idea of the separation
    root_filenames = [f for f in argv[1:]]
    TwoHits = False
    for file_name in root_filenames:
        lower_bound = float(file_name.split("_")[1][:-2])
        upper_bound = float(file_name.split("_")[-1][:-7])
        file = ur.open(file_name)
        if TwoHits:
            nPeaksInToF = nPeakInToF(file, 2) #we only want one hit
        else:
            nPeaksInToF = nPeakInToF(file, 1)

        for detector in range(len(file.keys()[:-1])):
            key = file.keys()[detector]
            print(key)
            df = file[key].arrays(library="pd")
            df = df[nPeaksInToF]
            df = df.reset_index() #very important!

            print(df.columns)


            intCharge = pd.DataFrame(df['IntCharge'].values.tolist())
            dfPeak = pd.concat([df, intCharge], axis = 1)

            windowIntCharge = pd.DataFrame(df['WindowIntCharge'].values.tolist())
            dfWindow = pd.concat([df, windowIntCharge], axis = 1)

            plt.figure(detector)
            plt.title('%s'%file.keys()[detector])
            plt.hist(dfWindow[0], bins = 100, range = (0, 0.75), label = 'Window Charge (%.0fns, %.0fns)'%(lower_bound, upper_bound), histtype = 'step')
            plt.hist(dfPeak[0], bins = 100, range = (0, 0.75), color = 'r', label = 'Peak Charge', histtype = 'step')
            plt.xlabel('Integrated charge [nC]')
            plt.ylabel('Triggers')

            # plt.show()
    for detector in range(len(file.keys()[:-1])):
        plt.figure(detector)
        plt.grid()
        plt.legend()
    plt.show()



if __name__ == "__main__":
    # studyLenWindow(sys.argv)
    overlayACTprofiles(sys.argv)
