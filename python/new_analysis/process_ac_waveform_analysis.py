#!/usr/bin/env python3
#
# Process waveform analysis of root file
# Usage:
#    python3 process_waveform_analysis.py [input_files] [config_file] [output_file]
#        input_file, input_file2, ...: any number of root files from midas2root
#        config_file: config JSON file
#        output_file: output ntuple root file

import sys

import numpy as np
import uproot as ur
import json5
from waveform_analysis_ac import WaveformAnalysis
import awkward as ak
import pandas as pd
import matplotlib.pyplot as plt

import coincidence as coincidence

def compare_lengths(arr1, arr2):
            """
            Recursively prints the lengths of nested arrays within an Awkward Array.
            """
            i = 0
            while len(arr1[i]) == len(arr2[i]) and i<len(arr1)-1:
                i+=1
            if i==len(arr1)-1:
                print("same dimensions!")
                return 0
            else:
                print("Dimension of array 1:", len(arr1[i]), "Dimension of array 2:",len(arr2[i]))
                print("Failed at %.1f percent i= %i"%(i/len(arr1)*100, i))
                return i


def print_usage(argv):
    print('Usage:')
    print(f"{argv[0]} input_root_file input_ntuple_file config_file output_file")
    print(f"    input_root_file one root files from midas2root")
    print(f"    input_ntuple_file one root file which contains multi-peak info")
    print(f"    config_file: config JSON file")
    print(f"    output_file: output ntuple root file")
    print(f"    ")
    print(f"    If you are having issues, try to swith on the debug variable in the process_ac_waveform_analysis.py code")
    return


def process_waveforms(argv):

    global debug 
    debug = False #some potentially useful print outs

    if len(argv) < 4:
        print("Not enough arguments provided.")
        print_usage(argv)
        return

    root_filenames = argv[1]
    ntuple_filename = argv[2]
    config_filename = argv[-2]
    output_filename = argv[-1]



    print(f"Using config file {config_filename}")
    with open(config_filename) as file:
        config = json5.load(file)["WaveAnalysis"]

    output_file = ur.recreate(output_filename)
    print(f"{len(root_filenames)} input root files to process.")
    # for f in root_filenames:
    process_file(root_filenames, ntuple_filename, config, output_file)
    output_file.close()

#
def process_file(root_filename, ntuple_filename, config, output_file):

    
    print(f"Processing root file {root_filename}\n")
#     "open the root file"
    run_file = ur.open(root_filename)
    ref_file = ur.open(ntuple_filename)
    global checkCoincidence
    checkCoincidence = config["checkCoincidence"]
    mainCoincidencePMTa = config["mainCoincidencePMTa"]
    mainCoincidencePMTb = config["mainCoincidencePMTb"]
    looseCoincidencePMT = config["looseCoincidencePMT"]

    #TODO: make check coincidence and the PMTs used part of the config file (and merge the config files)
    print("You have chosen to set the CheckCoincidence to", checkCoincidence)
    if (checkCoincidence):
        print("The PMTs used for the first coincidence are %i and %i and the PMT used for second coincidence is %i (loose or tight depends on set-up) \n"%(mainCoincidencePMTa, mainCoincidencePMTb, looseCoincidencePMT))

    signalTimeBranch = "SignalTimeCorrected"
    run_number = int(root_filename[-11:-5])

    if not(checkCoincidence):
        #reding the timings for calibration, technically we shouldn't need following Arturo's modifications but good sanity check
        #for a narrower integration window but because the integration is made with respect to the actual
        #waveforms it doesn't work well for now, we would have to apply the timing corrections there
        #I don't think that it is a big deal because the window is set quite large so it does capture all the signal, even if it uses the old timing.

        df_TOF1_hitTimes = ref_file['TOF10'].arrays(library="pd")
        hit_times = pd.DataFrame(df_TOF1_hitTimes['%s'%signalTimeBranch].values.tolist())
        df_TOF1_hitTimes =  pd.concat([df_TOF1_hitTimes, hit_times], axis=1)

        df_TOF11 = ref_file['TOF11'].arrays(library="pd")
        hit_times = pd.DataFrame(df_TOF11['%s'%signalTimeBranch].values.tolist())
        df_TOF11 =  pd.concat([df_TOF11, hit_times], axis=1)

        df_TOF12 = ref_file['TOF12'].arrays(library="pd")
        hit_times = pd.DataFrame(df_TOF12['%s'%signalTimeBranch].values.tolist())
        df_TOF12 =  pd.concat([df_TOF12, hit_times], axis=1)

        df_TOF13 = ref_file['TOF13'].arrays(library="pd")
        hit_times = pd.DataFrame(df_TOF13['%s'%signalTimeBranch].values.tolist())
        df_TOF13 =  pd.concat([df_TOF13, hit_times], axis=1)

        # print(df_TOF1_hitTimes)

    else: #look at the coincidence
        #now instead check the coincidence
        #the order will depend on the ones you choose, be careful
        #By default we take the coincidence between two PMTs in TOF1 and loosely match it with another PMT in TOF0,
        #TODO: make this coincidence accross all PMTs (increase efficiency)
        #we have tight-loose coincidence for TOF1-TOF0 in LM setup
        #and tight-tight in TG run, can change default first_hodoscope_run
        #in coincidence.performCoincidence
        SignalTimeMatchedTOF1, SignalTimeMatchedTOF0 = coincidence.performCoincidence(ntuple_filename, mainCoincidencePMTa, mainCoincidencePMTb, looseCoincidencePMT, run_number)
        #set the reference timing for the hits as the time matched dataframe
        df_TOF1_hitTimes = SignalTimeMatchedTOF1.copy()
        max_column_labels = df_TOF1_hitTimes.idxmax(axis=1, skipna=True)
        #need to have the same format to calculate the TOF when computing the position of the window
        df_TOF0_hitTimes = SignalTimeMatchedTOF0.copy()
        #some events will not see any data
        max_column_labels = np.where(np.isnan(max_column_labels), 0, max_column_labels)

        # Create a new DataFrame with max column labels to have the number of matched hits
        df_TOF1_hitTimes["nPeaks"] = max_column_labels.astype('int32')+1

        # print("In coincidence mode here are the timing of the hits:", df_TOF1_hitTimes)



    digitizers = ["midas_data_D300", "midas_data_D301", "midas_data_D302", "midas_data_D303"]
    channels_per_digitizer = config["NumberOfChannels"] // len(digitizers)
    digitizer_channels = {d: {} for d in digitizers}
    for i in range(config["NumberOfChannels"]):
        channel = f"Channel{i % channels_per_digitizer}"
        digitizer = digitizers[i // channels_per_digitizer]
        if config["ActiveChannels"][i]:
            if config['ChannelNames'][i] is None:
                raise ValueError(f"Channel number {i}: {digitizer}/{channel} is active, but has no channel name in the config file.")
            digitizer_channels[digitizer][f"{channel}"] = i
        print(f"Channel number {i}: {digitizer}/{channel} is {'active' if config['ActiveChannels'][i] else 'not active'}")
    active_digitizers = [d for d in digitizers if len(digitizer_channels[d]) > 0]
    print("active_digitisers", active_digitizers)
    events_per_digitizer = [run_file[d].num_entries for d in active_digitizers]
    if min(events_per_digitizer) != max(events_per_digitizer):
        raise ValueError("Digitizers wth active channels have different numbers of events: " +
                         ', '.join([f'{d}:{n}' for d, n in zip(active_digitizers, events_per_digitizer)]))
    n_channels = sum(len(channels) for channels in digitizer_channels.values())
    n_events = events_per_digitizer[0]
    if n_events == 0:
        print(f"WARNING: The digitizers in {root_filename} have zero events. Skipping this file.")
        return

    print(f"All {n_channels} active channels have {n_events} events")

    #we need to perform the calibration of the timing for the window offset
    
    for d in active_digitizers:
        optional_branch_names = [b for b in ["timeStamp", "triggerTime", "spillNumber"] if b in run_file[d].keys()]
        channels = digitizer_channels[d]
        branch_names = list(channels.keys()) + optional_branch_names
        print(f"Reading digitizer {d}")
        for batch, report in run_file[d].iterate(branch_names, step_size="100 MB", report=True):
            print(f"... events {report.tree_entry_start} to {report.tree_entry_stop}")
            for c, i in channels.items():
                print(f"... ... processing {d}/{c} into {config['ChannelNames'][i]}", end="", flush=True)
                print("\nWindow integration bounds are: [%.1f, %.1f] and the timing branch considered is %s. Coincidence check is %s\n"%(config["window_lower_bound"][i], config["window_upper_bound"][i], signalTimeBranch, checkCoincidence))
                #calculate the time offset to the reference
                df = ref_file[config['ChannelNames'][i]].arrays(library="pd")
                hit_times = pd.DataFrame(df['%s'%signalTimeBranch].values.tolist())
                #absolute hit time offset
                df =  pd.concat([df, hit_times], axis=1)
                if not(checkCoincidence):
                    #the wire length calibration that need to be made is the mean of the difference between the PMTs hit time and the TOF1 
                    #mean hit time (taking the earliest hit)
                    WindowTimeOffset = (df[0]-(df_TOF1_hitTimes[0]+df_TOF11[0]+df_TOF12[0]+df_TOF13[0])/4).mean()

                else:
                    #If there is coincidence, no need to offset, aligment has been done by Arturo's timing corrections
                    #We still remove the mean of the offset
                    #with some discretisation to make sure we are where we should be 
                    
                    WindowTimeOffset = (df[0]-df_TOF1_hitTimes[0]).mean()

                    #Here, the expected time offset is proportional to the particle travel time 
                    #(df[0] - SignalTimeMatchedTOF1[0]).mean()

                config_args = {
                    "threshold": config["Thresholds"][i],
                    "analysis_window": (config["AnalysisWindowLow"][i], config["AnalysisWindowHigh"][i]),
                    "pedestal_window": (config["PedestalWindowLow"][i], config["PedestalWindowHigh"][i]),
                    "reverse_polarity": (config["Polarity"][i] == 0),
                    "voltage_scale": config["VoltageScale"],
                    "time_offset": config["TimeOffset"][i],
                    "window_time_offset": WindowTimeOffset,
                    "PMTgain": config["PMTgain"][i],
                    "window_lower_bound": config["window_lower_bound"][i],
                    "window_upper_bound": config["window_upper_bound"][i],
                    "PMTdistanceToTOF1": config["PMTLocRelativeToTOF1"][i],
                    "distanceTOF1toTOF0": config["distanceTOF1toTOF0"],
                    "checkCoincidence": config["checkCoincidence"]
                }
                waveforms = batch[c]
                optional_branches = {b: batch[b] for b in optional_branch_names}

                if not(checkCoincidence):
                    process_batch(waveforms, optional_branches, config["ChannelNames"][i], output_file, config_args, df_TOF1_hitTimes.loc[int(report.tree_entry_start):int(report.tree_entry_stop)-1], df.loc[int(report.tree_entry_start):int(report.tree_entry_stop)-1],config["window_lower_bound"][i], config["window_upper_bound"][i])

                if checkCoincidence:
                    #reshape the timing entires so we can save them with the data
                    SignalTimeMatchedTOF0_reduced = np.array(SignalTimeMatchedTOF0.loc[int(report.tree_entry_start):int(report.tree_entry_stop)-1].values.squeeze().tolist())

                    SignalTimeMatchedTOF1_reduced  = np.array(SignalTimeMatchedTOF1.loc[int(report.tree_entry_start):int(report.tree_entry_stop)-1].values.squeeze().tolist())

                    #match the size of the dataframes to the size of the batch
                    df_TOF1_hitTimes_reduced = df_TOF1_hitTimes.loc[int(report.tree_entry_start):int(report.tree_entry_stop)-1]


                    df_TOF0_hitTimes_reduced = df_TOF0_hitTimes.loc[int(report.tree_entry_start):int(report.tree_entry_stop)-1]

                    df_reduced = df.loc[int(report.tree_entry_start):int(report.tree_entry_stop)-1]


                    SignalTimeMatchedTOF1_reduced = ak.Array([[val for val in SignalTimeMatchedTOF1_reduced[index] if (not np.isnan(val))] for index in range(len(SignalTimeMatchedTOF1_reduced))])

                    SignalTimeMatchedTOF0_reduced = ak.Array([[val for val in SignalTimeMatchedTOF0_reduced[index] if (not np.isnan(val))] for index in range(len(SignalTimeMatchedTOF0_reduced))])


                    process_batch(waveforms, optional_branches, config["ChannelNames"][i], output_file, config_args, df_TOF1_hitTimes_reduced, df_reduced,  config["window_lower_bound"][i], config["window_upper_bound"][i], df_TOF0_hitTimes_reduced, SignalTimeMatchedTOF1_reduced, SignalTimeMatchedTOF0_reduced)

    print(f"Saving event info for run {run_number}")
    output_file["EventInfo"] = {
        "RunNumber": np.repeat(run_number, n_events),
        "EventNumber": run_file[digitizers[0]]["eventNumber"].array(),
        "SpillNumber": run_file[digitizers[0]]["spillNumber"].array()
    }
    run_file.close()


def read_nTuple(df, branch_name, isSqueezed = False):
    #reading the information from the pre-processed ntuple to limit the computing time and not overwrite entries branches that are 1D (i.e. not arrays) need to be squeezed
    if isSqueezed:
        branch = ak.Array(np.array(pd.DataFrame(df['%s'%branch_name].values.tolist()).values.tolist()).squeeze())
    else:
        #faster option
        branch = pd.DataFrame(df['%s'%branch_name].values.tolist()).values.tolist()
        branch = ak.Array([[value for value in row if not pd.isna(value)] for row in branch])
        #slower option
        # branch = ak.Array(df['%s'%branch_name].dropna().tolist())
    return branch


def process_batch(waveforms, optional_branches, channel, output_file, config_args, df_TOF1_hitTimes, df, lower_bound, upper_bound, df_TOF0_hitTimes=None, SignalTimeMatchedTOF1 = None, SignalTimeMatchedTOF0 = None):

    waveform_analysis = WaveformAnalysis(waveforms, **config_args)

    #instead of overwriting the variables, copy them from the root file
    #TODO: impose coincidence and only perform the integration on peaks that have been matched
    if checkCoincidence:

        c =2.99792458 #ns/m
        electronTravelTime = c * config_args["PMTdistanceToTOF1"]

        #all particles are assumed to travel at speed c 
        #but protons and others are slower, we need to take those into account
        relativeTOF_offset = (df_TOF1_hitTimes - df_TOF0_hitTimes) * (config_args["PMTdistanceToTOF1"]/config_args["distanceTOF1toTOF0"]) - electronTravelTime

        #add the TOF offset to the expected time, which will be diferent for each peak
        for i in range(max(df_TOF1_hitTimes["nPeaks"])):
            df_TOF1_hitTimes[i] = relativeTOF_offset[i] + df_TOF1_hitTimes[i]

        waveform_analysis.run_window_analysis(df_TOF1_hitTimes, df)
        ##calculate the time of the matched hits in the waveform
        #only save for the TOF the hits that 
        shifted_Hit_Time = SignalTimeMatchedTOF1 - df["DigiTimingOffset"] + config_args["window_time_offset"]

        relativeTOF_offset_forTOFcalc = (SignalTimeMatchedTOF1 - SignalTimeMatchedTOF0) * (config_args["PMTdistanceToTOF1"]/config_args["distanceTOF1toTOF0"]) - electronTravelTime

        shifted_Hit_Time_upper = shifted_Hit_Time + upper_bound + relativeTOF_offset_forTOFcalc
        shifted_Hit_Time_lower = shifted_Hit_Time + lower_bound + relativeTOF_offset_forTOFcalc
        
        EndWaveform = waveform_analysis.waveformEnd
        ns_per_sample = waveform_analysis.ns_per_sample
        #we need at least part of the window to be included in the waveform
        #We need at least one entry in the waveform to perform the  integration so we 
        #are getting rid of all other entries in the matched timing by calculating if the expected timing is in it or not
        isAtLeastPartlyWithinRange = np.where(shifted_Hit_Time_upper<ns_per_sample, False, np.where(shifted_Hit_Time_lower>(EndWaveform)*ns_per_sample, False, True))

        SignalTimeMatchedTOF1 = SignalTimeMatchedTOF1[isAtLeastPartlyWithinRange]
        SignalTimeMatchedTOF0 = SignalTimeMatchedTOF0[isAtLeastPartlyWithinRange]
    else:
        waveform_analysis.run_window_analysis(df_TOF1_hitTimes, df)
    
    
    #window integration which has to be done here as it uses the timing
    window_integrated_charges = waveform_analysis.window_int_charge
    window_integrated_pe = waveform_analysis.window_int_pe

    window_width = waveform_analysis.windowWidth
    window_lower_bound = waveform_analysis.windowIntLowerBound
    window_upper_bound = waveform_analysis.windowIntUpperBound
    window_central_time = waveform_analysis.windowIntCentralTime
    window_central_time_corrected = waveform_analysis.windowIntCentralTimeCorrected


    peaks = ak.zip({"PeakVoltage": read_nTuple(df, "PeakVoltage", isSqueezed = False),
                    "PeakTime": read_nTuple(df, "PeakTime", isSqueezed = False),
                    "SignalTime": read_nTuple(df, "SignalTime", isSqueezed = False),
                    "IntCharge": read_nTuple(df, "IntCharge", isSqueezed = False),
                    "IntPE": read_nTuple(df, "IntPE", isSqueezed = False),
                    "SignalTimeCorrected": read_nTuple(df, "SignalTimeCorrected", isSqueezed = False)})

    #needs to be handled separately because it has another 
    if not(checkCoincidence):
        window_peaks = ak.zip({"WindowIntCharge": window_integrated_charges,
                            "WindowIntPE": window_integrated_pe,
                            "WindowLowerTime": window_lower_bound,
                            "WindowUpperTime":
                            window_upper_bound,
                            "WindowCentralTime":
                            window_central_time,
                            "WindowCentralTimeCorrected":
                            window_central_time_corrected})

        #if we do not have TOF defined.
        branches = {"Pedestal": read_nTuple(df, "Pedestal", isSqueezed = True),
                    "PedestalSigma": read_nTuple(df, "PedestalSigma", isSqueezed = True),
                    "MaxVoltage": read_nTuple(df, "MaxVoltage", isSqueezed = True),
                    "WholeWaveformInt": read_nTuple(df, "WholeWaveformInt", isSqueezed = True),
                    "DigiTimingOffset": read_nTuple(df, "DigiTimingOffset", isSqueezed = True),
                    "Peaks": peaks,
                    "WindowPeaks": window_peaks,
                    **optional_branches}


    else:
        #trying to compare the shapes of the array to see where there are issues
        if debug:
            #somethimes we have some issued with the format, should not have them anymore 
            #but just in case, useful to keep
            print("Comparing the length of SignalTimeMatchedTOF0 and SignalTimeMatchedTOF1")
            compare_lengths(SignalTimeMatchedTOF0, SignalTimeMatchedTOF1)
            print("Comparing the length of SignalTimeMatchedTOF0 and window_width")
            k = compare_lengths(SignalTimeMatchedTOF0, window_width)

            print(shifted_Hit_Time_lower[k], shifted_Hit_Time_upper[k])

            print(SignalTimeMatchedTOF1[k], window_width[k])
            print("Comparing the length of window_width and window_integrated_charges")
            compare_lengths(window_width, window_integrated_charges)

            
            print("Comparing the length of SignalTimeMatchedTOF1 and window_integrated_charges")
            compare_lengths(SignalTimeMatchedTOF1, window_integrated_charges)

        window_peaks = ak.zip({
                            "WindowIntCharge": window_integrated_charges,
                            "WindowIntPE": window_integrated_pe,
                            # "df_TOF1_hitTimes":df_TOF1_hitTimes,
                            #"nPeaksMatched": df_TOF1_hitTimes["nPeaks"],
                            "WindowWidth": window_width,
                            "SignalTimeMatchedTOF1": SignalTimeMatchedTOF1,
                            "SignalTimeMatchedTOF0": SignalTimeMatchedTOF0,
                            "WindowLowerTime": window_lower_bound,
                            "WindowUpperTime":
                            window_upper_bound,
                            "WindowCentralTime":
                            window_central_time,
                            "WindowCentralTimeCorrected":
                            window_central_time_corrected
                            })


        branches = {"Pedestal": read_nTuple(df, "Pedestal", isSqueezed = True),
                    "PedestalSigma": read_nTuple(df, "PedestalSigma", isSqueezed = True),
                    "MaxVoltage": read_nTuple(df, "MaxVoltage", isSqueezed = True),
                    "WholeWaveformInt": read_nTuple(df, "WholeWaveformInt", isSqueezed = True),
                    "WholeWaveformIntPE": read_nTuple(df, "WholeWaveformIntPE", isSqueezed = True),
                    "DigiTimingOffset": read_nTuple(df, "DigiTimingOffset", isSqueezed = True),
                    "Peaks": peaks,
                    "WindowPeaks": window_peaks,
                    **optional_branches}

    print(f" ... writing {channel} to {output_file.file_path}")

    try:
        output_file[channel].extend(branches)
    except KeyError:
        branch_types = {name: branch.type for name, branch in branches.items()}
        output_file.mktree(channel, branch_types, field_name=(lambda s, f: f))
        output_file[channel].extend(branches)



if __name__ == "__main__":
    # execute only if run as a script"
    process_waveforms(sys.argv)
