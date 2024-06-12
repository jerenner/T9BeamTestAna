#this is a code to ensure coincidence between the signal accros the 4 PMTs of a given TOF detector
import sys
import pandas as pd
import numpy as np
import uproot as ur
import matplotlib.pyplot as plt
import scipy.optimize as spo


def isInCoincidence(df_substracted, expectedMean = 0, expectedStd = 10, nSigma = 3):
    mask = (df_substracted > expectedMean - nSigma * expectedStd) & (df_substracted < expectedMean + nSigma * expectedStd)
    return mask

def isInLooseCoincidence(df_substracted, expectedMean = 0, expectedStd = 10, protonMaxTOF = 40, deltaTOFelectronMuons = 40, nSigma = 5):
    '''once the matching has been done within one TOF we can roughly match them to the second TOF detector, need to accept protons'''
    if (expectedMean>0) :
        '''the order of the TOFs is TOF1 - TOF0'''
        mask = (df_substracted > expectedMean - nSigma * expectedStd - protonMaxTOF) & (df_substracted < expectedMean + nSigma * expectedStd)
    else:
        '''The order is reversed '''
        mask = (df_substracted > expectedMean - nSigma * expectedStd - protonMaxTOF) & (df_substracted < expectedMean + nSigma * expectedStd)

    return mask

def getCoincidenceAcrossTOFs(hit_timing, TOFa_base, TOFa_test, TOFb_test, runNumber, first_hodoscope_run = 579):
    matched_event_times, good_candidate_times, mean_matched_time = getMatchedEventsTiming(hit_timing, TOFa_base, TOFa_test, False)


    #correct the relevant branch
    hit_timing[TOFa_base] = mean_matched_time
    # check loose coincidence for protons, deuterium in LM setup
    if runNumber<first_hodoscope_run:
        matched_event_times_TOFb, matched_event_times_TOFa, mean_matched_time = getMatchedEventsTiming(hit_timing, TOFb_test, TOFa_base, True)
    else:#keep a loose coincidence for now (keep protons!!!)
        matched_event_times_TOFb, matched_event_times_TOFa, mean_matched_time = getMatchedEventsTiming(hit_timing, TOFb_test, TOFa_base, True)

    #in case we only want one coincidence
    # matched_event_times_TOFb, matched_event_times_TOFa = matched_event_times, good_candidate_times
    #careful, there is an intended switch in the order
    return matched_event_times_TOFa, matched_event_times_TOFb



def get2plus2CoincidenceAcrossTOFs(hit_timing, TOFa_base, TOFa_test, TOFb_test, TOFb_base, runNumber, first_hodoscope_run = 579):
    
    #TS1 1+1 coincidence
    matched_event_times, good_candidate_times, mean_matched_time = getMatchedEventsTiming(hit_timing, TOFa_base, TOFa_test, False)

    #TS0 1+1 coincidence
    matched_event_timesb, good_candidate_timesb, mean_matched_timeb = getMatchedEventsTiming(hit_timing, TOFb_base, TOFb_test, False)

    #correct the relevant TS1 branch
    hit_timing[TOFa_base] = mean_matched_time
    #correct the relevant TS0 branch
    hit_timing[TOFb_test] = mean_matched_timeb

    # check loose coincidence for protons, deuterium in LM setup
    if runNumber<first_hodoscope_run:
        matched_event_times_TOFb, matched_event_times_TOFa, mean_matched_time = getMatchedEventsTiming(hit_timing, TOFb_test, TOFa_base, True)

    else:#keep a loose coincidence for now (keep protons!!!) can decide to remove them later
        matched_event_times_TOFb, matched_event_times_TOFa, mean_matched_time = getMatchedEventsTiming(hit_timing, TOFb_test, TOFa_base, True)

    #in case we only want one coincidence
    # matched_event_times_TOFb, matched_event_times_TOFa = matched_event_times, good_candidate_times
    #careful, there is an intended switch in the order
    return matched_event_times_TOFa, matched_event_times_TOFb





def getMatchedEventsTiming(hit_timing, basePMT, testPMT, matchingTOF0andTOF1 = False):
    matched_event_times = pd.DataFrame()
    mean_matched_time = pd.DataFrame()
    good_candidate_times = pd.DataFrame()

    #calculate mean delay between the two PMTs
    # print(hit_timing)
    mean, std = float(pd.DataFrame(hit_timing[basePMT][0]-hit_timing[testPMT][0]).dropna().mean()), float(pd.DataFrame(hit_timing[basePMT][0]-hit_timing[testPMT][0]).dropna().std())

    nbins = 1000
    alpha = 0.3
    rangeLow = -100
    rangeHigh = 100
    #fit a gaussian to the data to have a smaller std
    tof0_tof2, bins = np.histogram(hit_timing[basePMT][0]-hit_timing[testPMT][0], bins = nbins, range = (rangeLow, rangeHigh))

    bins_centre = (bins[:-1]+bins[1:])/2
    popt, pcov = spo.curve_fit(Gaussian, bins_centre, tof0_tof2)
    mean = popt[1]
    std = abs(popt[2])

    # print(mean, std)
    # mean = -0.02
    # std = 6

    # print("\ninitial hit time in PMT %i: \n"%basePMT, hit_timing[basePMT], "\nintial hit times in PMT %i:\n"%testPMT, hit_timing[testPMT])

    #for now look at only the first matched hit -> could look at other ones later
    #but I am having hit multiplicities issues TODO
    largerBaseHit = max(hit_timing[basePMT].columns)

    #run through all the events in the testPMT frame and try to find matches
    for candidate_hit in hit_timing[testPMT].columns:

        #this is the canditate that we will try to match in the other PMT
        time_difference = hit_timing[basePMT].subtract(hit_timing[testPMT][candidate_hit], axis = 0)
        #this is the reference event that we will try to match with the base
        candidate_times = hit_timing[testPMT][candidate_hit]

        #this makes a mask of the base events saying if they are in coincidence or not
        if not(matchingTOF0andTOF1):
            coincidence = isInCoincidence(time_difference, mean, std)
        else:
            coincidence = isInLooseCoincidence(time_difference, mean, std)


        #for each event we select the base hit time correcsponding to that target
        matched_timing = pd.DataFrame(hit_timing[basePMT][coincidence].reset_index(drop=True))
        # #need to reset the index to the reference, technically not essential but safer
        matched_timing.index = candidate_times.index

        # print("Sum of the events in the 3 columns of matched events:", len(matched_timing[0].dropna())+len(matched_timing[1].dropna())+len(matched_timing[2].dropna()))

        #now we are removing the nans and collapsing all of the columns into one so only the earliest matched hit remains or nan
        # print(largerBaseHit, max(matched_timing.columns))
        if (largerBaseHit==2 and max(matched_timing.columns)>=2): #candidate_hit == 0 and candidate_hit<= largerBaseHit):
            matched_timing_compressed = matched_timing[0].combine_first(matched_timing[1]).combine_first(matched_timing[2])
        elif (largerBaseHit==3 and max(matched_timing.columns)>=3): #candidate_hit == 0 and candidate_hit<= largerBaseHit):
            matched_timing_compressed = matched_timing[0].combine_first(matched_timing[1]).combine_first(matched_timing[2]).combine_first(matched_timing[3])
        elif (largerBaseHit==4 and max(matched_timing.columns)>=4): #candidate_hit == 0 and candidate_hit<= largerBaseHit):
            matched_timing_compressed = matched_timing[0].combine_first(matched_timing[1]).combine_first(matched_timing[2]).combine_first(matched_timing[3]).combine_first(matched_timing[4])
        elif (largerBaseHit==5 and max(matched_timing.columns)>=5): #candidate_hit == 0 and candidate_hit<= largerBaseHit):
            matched_timing_compressed = matched_timing[0].combine_first(matched_timing[1]).combine_first(matched_timing[2]).combine_first(matched_timing[3]).combine_first(matched_timing[4]).combine_first(matched_timing[5])
        else:
            raise Exception("DimensionDisparity")


        # elif (candidate_hit ==1 and candidate_hit<= largerBaseHit ):
        #     matched_timing_compressed = matched_timing[1].combine_first(matched_timing[2]).combine_first(matched_timing[3])
        #
        # elif (candidate_hit ==2 and candidate_hit<= largerBaseHit):
        #     matched_timing_compressed = matched_timing[2].combine_first(matched_timing[3])
        #
        # elif (candidate_hit == 3 and candidate_hit<= largerBaseHit):
        #     matched_timing_compressed = matched_timing[3]
        #
        # elif (candidate_hit == 4 and candidate_hit<= largerBaseHit):
        #     matched_timing_compressed = matched_timing[4]
        #
        # elif (candidate_hit ==5 and candidate_hit<= largerBaseHit):
        #     matched_timing_compressed = matched_timing[5]

        #for that matched event, get the mean timing for that hit
        matched_time = matched_timing_compressed

        # print("The length after compression:", len(matched_time.dropna()))

#         #a candidate is a good candidate and should be kept only if it has a match
#         if not(matchingTOF0andTOF1):
#             good_candidate_times[candidate_hit] = np.where(coincidence.any(axis = 1), candidate_times, np.nan)
#             good_candidate_times.index = candidate_times.index
#         else:
#             # good_candidate_times = hit_timing[testPMT]
#

        #We need a time of flight measurement, so we can apply the same method for this new timing and the hit time in the TOF 1


        #append the timing of this nth matched event to the new dataframe
        if (candidate_hit<= max(hit_timing[testPMT].columns) and candidate_hit<= max(hit_timing[basePMT].columns) and candidate_hit<=largerBaseHit):
            good_candidate_times[candidate_hit] = np.where(coincidence.any(axis = 1), candidate_times, np.nan)
            good_candidate_times.index = candidate_times.index

            # print(len(good_candidate_times[candidate_hit].dropna()))

            matched_event_times[candidate_hit] = matched_time
            # print(len(matched_event_times[candidate_hit].dropna()))


            mean_matched_time[candidate_hit] = (matched_time + good_candidate_times[candidate_hit])/2


    # print( "\nmatched times in PMT%i:\n"%basePMT, matched_event_times, "\n matched times in PMT %i:\n"%testPMT, good_candidate_times)

    # print(len(matched_event_times[1]), len(matched_event_times[0].dropna()), len(matched_event_times[2].dropna()), len(matched_event_times[3].dropna()))
    return matched_event_times, good_candidate_times, mean_matched_time

def performCoincidence(argv, TOFa_base = 0, TOFa_test = 1, TOFb_test = 6, runNumber = 0 , checkWarning = True, plot = False, first_hodoscope_run= 579):
    root_filenames = argv #[1]
    #read in the file
    # print("Hello from coincidence")
    file = ur.open(root_filenames)

    run = root_filenames[-8:-5]

    signalTimeTarget = 'SignalTimeCorrected'
    targetTOF = 'TOF'

    hit_timing = []
    #append the data to the initial dataframe
    for detector in range(len(file.keys()[:-1])):
        key = file.keys()[detector]
        df = file[key].arrays(library="pd")

        list_index = df.index
        df = df.reset_index()
        #only store the TOF detectors
        if key[:-4] == '%s'%targetTOF:
            digitiser0 = file[file.keys()[0]].arrays(library="pd")
            signalTimeCorrected =  pd.DataFrame(df['%s'%signalTimeTarget].values.tolist())

            hit_timing.append(signalTimeCorrected) #.iloc[200:300]#829:849

    key = file.keys()[-1]
    df_enventInfo =  file[key].arrays(library="pd")

    if plot:
        plotting(hit_timing)

    #take PMTs far away from each other to be as stable as possible


    if ((TOFa_base > 3) and (TOFa_test>3) and (TOFb_test<3)):
        print("You have paired up the PMTs correctely")
    elif (checkWarning):
        raise Exception("The PMTs aren't paired up properly, %i and %i should be on TOF1, i.e. above 3 and %i should be on the other side of 3 i.e. on the other TOF which is is not"%(TOFa_base,TOFa_test, TOFb_test))

    matched_event_times_TOFa, matched_event_times_TOFb = getCoincidenceAcrossTOFs(hit_timing, TOFa_base, TOFa_test, TOFb_test, runNumber, first_hodoscope_run)



    # print(len(matched_event_times_TOFa[0]), len(matched_event_times_TOFa[0].dropna()))

    return matched_event_times_TOFa, matched_event_times_TOFb

def perform2plus2Coincidence(argv, TOFa_base = 0, TOFa_test = 1, TOFb_test = 6, TOFb_base = 5, runNumber = 0 , checkWarning = True, plot = False, first_hodoscope_run= 579):
    """"Perform coincidence across paris of PMTs in each TS before doing it accros them"""
    root_filenames = argv #[1]
    #read in the file
    # print("Hello from coincidence")
    file = ur.open(root_filenames)

    run = root_filenames[-8:-5]

    signalTimeTarget = 'SignalTimeCorrected'
    targetTOF = 'TOF'

    hit_timing = []
    #append the data to the initial dataframe
    for detector in range(len(file.keys()[:-1])):
        key = file.keys()[detector]
        df = file[key].arrays(library="pd")

        list_index = df.index
        df = df.reset_index()
        #only store the TOF detectors
        if key[:-4] == '%s'%targetTOF:
            digitiser0 = file[file.keys()[0]].arrays(library="pd")
            signalTimeCorrected =  pd.DataFrame(df['%s'%signalTimeTarget].values.tolist())

            hit_timing.append(signalTimeCorrected) #.iloc[200:300]#829:849

    key = file.keys()[-1]
    df_enventInfo =  file[key].arrays(library="pd")

    if plot:
        plotting(hit_timing)

    #take PMTs far away from each other to be as stable as possible


    if ((TOFa_base > 3) and (TOFa_test>3) and (TOFb_test<3)):
        print("You have paired up the PMTs correctely")
    
    elif (checkWarning):
        raise Exception("The PMTs aren't paired up properly, %i and %i should be on TOF1, i.e. above 3 and %i should be on the other side of 3 i.e. on the other TOF which is is not"%(TOFa_base,TOFa_test, TOFb_test))

    matched_event_times_TOFa, matched_event_times_TOFb = get2plus2CoincidenceAcrossTOFs(hit_timing, TOFa_base, TOFa_test, TOFb_test, TOFb_base, runNumber, first_hodoscope_run)

    # print(len(matched_event_times_TOFa[0]), len(matched_event_times_TOFa[0].dropna()))

    return matched_event_times_TOFa, matched_event_times_TOFb

def Gaussian(x, A, mu, sigma):
    return A * np.exp(-(x-mu)**2/(2*sigma**2))


def plotting(hit_timing):
    df_TOF = hit_timing
    targetTOF = "TOF1"
    signalTimeTarget = 'SignalTime'
    nbins = 50
    alpha = 0.3
    rangeLow = -3
    rangeHigh = 3
    run = 393

    # np.hist(df_TOF[0][0]-df_TOF[1][0], bins = 100, range = (-2.5, 2.5), alpha = 0.5)
    #need to afit a gaussian!!
    plt.figure(figsize = (10, 15))
    for basePMT in range(2):
        for testPMT in range(basePMT+1, 2):
            print('PMT %i - PMT %i'%(basePMT, testPMT))
            mean, std = np.array(df_TOF[basePMT][0]-df_TOF[testPMT][0]).mean(), np.array(df_TOF[basePMT][0]-df_TOF[testPMT][0]).std()
            tof0_tof2, bins, _ = plt.hist(df_TOF[basePMT][0]-df_TOF[testPMT][0], bins = nbins, range = (rangeLow, rangeHigh), alpha = alpha, label = '(%s) PMT %i - PMT %i: mean %.2f ns, std: %.2f ns'%(targetTOF, basePMT, testPMT, mean, std))

            bins_centre = (bins[:-1]+bins[1:])/2

            plt.plot(bins_centre, tof0_tof2, 'rx')

            popt, pcov = spo.curve_fit(Gaussian, bins_centre, tof0_tof2)
            plt.plot(bins_centre, Gaussian(bins_centre, *popt), 'k--', label='Gaussian fit: A=%5.3f, mean=%5.3f, std=%5.3f' % tuple(popt))




    plt.grid()
    plt.legend()
    plt.xlabel("Delay between timing of earliest hit in pairs of %s PMTs (ns)"%targetTOF,  fontsize = 15)
    plt.ylabel("Occurences",  fontsize = 12)
    plt.title("Run %s: relative timing between \n pairs of PMTs in %s - %s"%(run, targetTOF, signalTimeTarget), weight = 'bold', fontsize = 12)

    # plt.savefig('/home/ac4317/Laptops/Year1/WCTE/BeamTestJuly2023/DataAnalysis/data/analysis-code/T9BeamTestAna/pdf_results/relative_Timings_%s_%s_coincidence_run%s.pdf'%(targetTOF, signalTimeTarget, run))

    plt.show()




if __name__ == "__main__":
    """If we call the code directly, for testing"""
    TOFa_base = 4
    TOFa_test = 5
    TOFb_test = 1
    TOF0, TOF1  = checkCoincidence(sys.argv[1], TOFa_base, TOFa_test, TOFb_test, False, False)
    #check the efficiecency of the selection
    accountedFor = 0
    for i in max(TOF0.columns)-TOF0.columns:
        print("The number of %i events (or more!) in coincidence in PMTs %i and %i and %i is: %.2f (%.2f percent), total length is %i"%(i+1, TOFa_base, TOFa_test, TOFb_test, len(TOF0[i].dropna()),  len(TOF0[i].dropna())/len(TOF0) * 100, len(TOF0)))

        print("Total number of exaclty %i events: %.2f (%.3f percent)"%(i+1, len(TOF0[i].dropna())-accountedFor, (len(TOF0[i].dropna())-accountedFor)/len(TOF0) * 100))
        accountedFor = len(TOF0[i].dropna())




#plotting
    # raise BelowPlottingOfTypicalDelay


