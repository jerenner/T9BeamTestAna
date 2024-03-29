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



def print_usage(argv):
    print('Usage:')
    print(f"{argv[0]} input_root_file input_ntuple_file config_file output_file")
    print(f"    input_root_file one root files from midas2root")
    print(f"    input_ntuple_file one root file which contains multi-peak info")
    print(f"    config_file: config JSON file")
    print(f"    output_file: output ntuple root file")
    return


def process_waveforms(argv):
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
        print("The PMTs used for the main coincidence are %i and %i and the PMT used for loose coincidence is %i \n"%(mainCoincidencePMTa, mainCoincidencePMTb, looseCoincidencePMT))

    signalTimeBranch = "SignalTimeCorrected" #eventually we should move to SignalTimeCorrected

    if not(checkCoincidence):
        #reding the timings for calibration, technically we shouldn't need following Arturo's
        #modifications but good sanity check


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

        print(df_TOF1_hitTimes)

    else: #look at the coincidence
        #now instead check the coincidence
        #the order will depend on the ones you choose, be careful
        #By default we take the coincidence between two PMTs in TOF1 and loosely match it with another PMT in TOF0,
        #TODO: make this coincidence accross all PMTs (increase efficiency)
        SignalTimeMatchedTOF1, SignalTimeMatchedTOF0 = coincidence.checkCoincidence(ntuple_filename, mainCoincidencePMTa, mainCoincidencePMTb, looseCoincidencePMT)
        #set the reference timing for the hits as the time matched dataframe
        df_TOF1_hitTimes = SignalTimeMatchedTOF1.copy()
        max_column_labels = df_TOF1_hitTimes.idxmax(axis=1, skipna=True)
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
    run_number = int(root_filename[-11:-5])
    for d in active_digitizers:
        optional_branch_names = [b for b in ["timeStamp", "triggerTime", "spillNumber"] if b in run_file[d].keys()]
        channels = digitizer_channels[d]
        branch_names = list(channels.keys()) + optional_branch_names
        print(f"Reading digitizer {d}")
        for batch, report in run_file[d].iterate(branch_names, step_size="100 MB", report=True):
            print(f"... events {report.tree_entry_start} to {report.tree_entry_stop}")
            for c, i in channels.items():
                print(f"... ... processing {d}/{c} into {config['ChannelNames'][i]}", end="", flush=True)
                #calculate the time offset to the reference
                df = ref_file[config['ChannelNames'][i]].arrays(library="pd")
                hit_times = pd.DataFrame(df['%s'%signalTimeBranch].values.tolist())
                #absolute hit time offset
                df =  pd.concat([df, hit_times], axis=1)
                if not(checkCoincidence):
                    #the wire length calibration that need to be made is the mean of the difference between the PMTs hit time and the TOF1 mean hit time (taking the earliest hit)
                    WindowTimeOffset = (df[0]-(df_TOF1_hitTimes[0]+df_TOF11[0]+df_TOF12[0]+df_TOF13[0])/4).mean()

                else:
                    #If there is coincidence, no need to offset, aligment has been done by Arturo's timing corrections
                    WindowTimeOffset = 0 #(df[0] - SignalTimeMatchedTOF1[0]).mean()

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
                }
                waveforms = batch[c]
                optional_branches = {b: batch[b] for b in optional_branch_names}

                if not(checkCoincidence):
                    process_batch(waveforms, optional_branches, config["ChannelNames"][i], output_file, config_args, df_TOF1_hitTimes.loc[int(report.tree_entry_start):int(report.tree_entry_stop)-1], df.loc[int(report.tree_entry_start):int(report.tree_entry_stop)-1])

                if checkCoincidence:
                    #reshape the timing entires so we can save them with the data
                    SignalTimeMatchedTOF0_reduced = np.array(SignalTimeMatchedTOF0.loc[int(report.tree_entry_start):int(report.tree_entry_stop)-1].values.squeeze().tolist())
                    #TODO get more than 1 hit !!
                    SignalTimeMatchedTOF0_reduced  = ak.Array([[val for val in values if (not np.isnan(val))] for values in SignalTimeMatchedTOF0_reduced])
                    # SignalTimeMatchedTOF0 = ak.Array([[values] for values in SignalTimeMatchedTOF0])

                    SignalTimeMatchedTOF1_reduced  = np.array(SignalTimeMatchedTOF1.loc[int(report.tree_entry_start):int(report.tree_entry_stop)-1].values.squeeze().tolist())
                    SignalTimeMatchedTOF1_reduced

                    SignalTimeMatchedTOF1_reduced  = ak.Array([[val for val in values if (not np.isnan(val))] for values in SignalTimeMatchedTOF1_reduced])
                    # SignalTimeMatchedTOF1 = ak.Array([[values] for values in SignalTimeMatchedTOF1])


                    process_batch(waveforms, optional_branches, config["ChannelNames"][i], output_file, config_args, df_TOF1_hitTimes.loc[int(report.tree_entry_start):int(report.tree_entry_stop)-1], df.loc[int(report.tree_entry_start):int(report.tree_entry_stop)-1], SignalTimeMatchedTOF1_reduced, SignalTimeMatchedTOF0_reduced)



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


def process_batch(waveforms, optional_branches, channel, output_file, config_args, df_TOF1_hitTimes, df, SignalTimeMatchedTOF1 =None, SignalTimeMatchedTOF0 = None):

    waveform_analysis = WaveformAnalysis(waveforms, **config_args)

    #instead of overwriting the variables, copy them from the root file
    #TODO: impose coincidence and only perform the integration on peaks that have been matched
    waveform_analysis.run_window_analysis(df_TOF1_hitTimes, df)
    #window integration which has to be done here as it uses the timing
    window_integrated_charges = waveform_analysis.window_int_charge
    window_integrated_pe = waveform_analysis.window_int_pe


    peaks = ak.zip({"PeakVoltage": read_nTuple(df, "PeakVoltage", isSqueezed = False),
                    "PeakTime": read_nTuple(df, "PeakTime", isSqueezed = False),
                    "SignalTime": read_nTuple(df, "SignalTime", isSqueezed = False),
                    "IntCharge": read_nTuple(df, "IntCharge", isSqueezed = False),
                    "IntPE": read_nTuple(df, "IntPE", isSqueezed = False),
                    "SignalTimeCorrected": read_nTuple(df, "SignalTimeCorrected", isSqueezed = False)})

    #needs to be handled separately because it has another format
    if not(checkCoincidence):
        window_peaks = ak.zip({"WindowIntCharge": window_integrated_charges,
                            "WindowIntPE": window_integrated_pe})


    else:
        # print(len(window_integrated_charges[:]), len(SignalTimeMatchedTOF1[:]))
        window_peaks = ak.zip({"WindowIntCharge": window_integrated_charges,
                            "WindowIntPE": window_integrated_pe,
                            "SignalTimeMatchedTOF1": SignalTimeMatchedTOF1,
                            "SignalTimeMatchedTOF0": SignalTimeMatchedTOF0,
                            # "nPeaksMatched": df_TOF1_hitTimes["nPeaks"],
                            })


    branches = {"Pedestal": read_nTuple(df, "Pedestal", isSqueezed = True),
                "PedestalSigma": read_nTuple(df, "PedestalSigma", isSqueezed = True),
                "MaxVoltage": read_nTuple(df, "MaxVoltage", isSqueezed = True),
                "WholeWaveformInt": read_nTuple(df, "WholeWaveformInt", isSqueezed = True),
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



    # for key in ref_file.keys()[:-1]:
    #     df = ref_file[key].arrays(library="pd")
    #     #it is useful to look at only the case where we have one peak in all of the TOF
    #     detectors_pd.append(df)
    #             # print(df[onePeakToF])


#     detectors = run_file.keys()
#     branch_names = run_file['TOF03'].keys()
#     digitizer_channels = {d: {} for d in detectors}
#
#     for d in detectors:
#         channels = digitizer_channels[d]
#         print(channels,digitizer_channels)
#         for batch, report in run_file[d].iterate(branch_names, step_size="1000 MB", report=True):
#             print(f"... events {report.tree_entry_start} to {report.tree_entry_stop}")
#             for c, i in channels.items():
#                 # print(f"... ... processing {d}/{c} into {config['ChannelNames'][i]}", end="", flush=True)
#                 config_args = {
#                     "threshold": config["Thresholds"][i],
#                     "analysis_window": (config["AnalysisWindowLow"][i], config["AnalysisWindowHigh"][i]),
#                     "pedestal_window": (config["PedestalWindowLow"][i], config["PedestalWindowHigh"][i]),
#                     "reverse_polarity": (config["Polarity"][i] == 0),
#                     "voltage_scale": config["VoltageScale"],
#                     "time_offset": config["TimeOffset"][i],
#                 }
#                 waveforms = batch[c]
#                 #all the branches are optional and we want to keep them stored
#                 optional_branches = {b: batch[b] for b in branch_names}
#                 process_batch(waveforms, optional_branches, config["ChannelNames"][i], output_file, config_args)
#
#     run_number = int(root_filename[-11:-5])
#     print(f"Saving event info for run {run_number}")
#     output_file["EventInfo"] = {
#         "RunNumber": np.repeat(run_number, n_events),
#         "EventNumber": run_file[digitizers[0]]["eventNumber"].array(),
#         "SpillNumber": run_file[digitizers[0]]["spillNumber"].array()
#     }
#     run_file.close()4
#
#
# def process_batch(waveforms, optional_branches, channel, output_file, config_args):
#     waveform_analysis = WaveformAnalysis(waveforms, **config_args)
#     waveform_analysis.run_analysis()
#     integrated_charges = waveform_analysis.window_int_charge
#     #
#     # pedestals = ak.Array(waveform_analysis.pedestals.squeeze())
#     # pedestal_sigmas = ak.Array(waveform_analysis.pedestal_sigmas.squeeze())
#     # peak_voltages = waveform_analysis.pulse_peak_voltages
#     # peak_times = waveform_analysis.pulse_peak_times
#     # signal_times = waveform_analysis.pulse_signal_times
#     # integrated_charges = waveform_analysis.pulse_charges
#
#     peaks = ak.zip({"PeakVoltage": peak_voltages,
#                     "PeakTime": peak_times,
#                     "SignalTime": signal_times,
#                     "IntCharge": integrated_charges})
#
#     branches = {**optional_branches}
#
#     print(f" ... writing {channel} to {output_file.file_path}")
#     if channel not in output_file.keys():
#         branch_types = {name: branch.type for name, branch in branches.items()}
#         output_file.mktree(channel, branch_types, counter_name=(lambda s: "nPeaks"), field_name=(lambda s, f: f))
#     output_file[channel].extend(branches)

if __name__ == "__main__":
    # execute only if run as a script"
    process_waveforms(sys.argv)
