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
    window_lower_bound = -15 #ns - the default value
    window_upper_bound = 35

    if len(argv) > 5:
        print("Using lower and upper bounds for the integration window.")
        print("Careful, lower bound need to be negative number")
        window_lower_bound = float(argv[3])
        window_upper_bound = float(argv[4])


    print("Using integration window: (%.2f,%.2f)"%(window_lower_bound, window_upper_bound))

    print(f"Using config file {config_filename}")
    with open(config_filename) as file:
        config = json5.load(file)["WaveAnalysis"]

    output_file = ur.recreate(output_filename)
    print(f"{len(root_filenames)} input root files to process.")
    # for f in root_filenames:
    process_file(root_filenames, ntuple_filename, config, output_file, window_lower_bound, window_upper_bound)
    output_file.close()

#
def process_file(root_filename, ntuple_filename, config, output_file, window_lower_bound, window_upper_bound):
    print(f"Processing root file {root_filename}")
#     "open the root file"
    run_file = ur.open(root_filename)


    ref_file = ur.open(ntuple_filename)

    #reding the timings for calibration, technically we shouldn't need following Arturo's
    #modifications but good sanity check

    signalTimeBranch = "SignalTime" #eventually we should move to SignalTimeCorrected
    #for a narrower integration window but because the integration is made with respect to the actual
    #waveforms it doesn't work well for now, we would have to apply the timing corrections there
    #I don't think that it is a big deal because the window is set quite large so it does capture all the signal, even if it uses the old timing.

    df_TOF10 = ref_file['TOF10'].arrays(library="pd")
    hit_times = pd.DataFrame(df_TOF10['%s'%signalTimeBranch].values.tolist())
    df_TOF10 =  pd.concat([df_TOF10, hit_times], axis=1)

    df_TOF11 = ref_file['TOF11'].arrays(library="pd")
    hit_times = pd.DataFrame(df_TOF11['%s'%signalTimeBranch].values.tolist())
    df_TOF11 =  pd.concat([df_TOF11, hit_times], axis=1)

    df_TOF12 = ref_file['TOF12'].arrays(library="pd")
    hit_times = pd.DataFrame(df_TOF12['%s'%signalTimeBranch].values.tolist())
    df_TOF12 =  pd.concat([df_TOF12, hit_times], axis=1)

    df_TOF13 = ref_file['TOF13'].arrays(library="pd")
    hit_times = pd.DataFrame(df_TOF13['%s'%signalTimeBranch].values.tolist())
    df_TOF13 =  pd.concat([df_TOF13, hit_times], axis=1)


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
                df =  pd.concat([df, hit_times], axis=1)
                #absolute hit time offset
                WindowTimeOffset = (df[0]-(df_TOF10[0]+df_TOF11[0]+df_TOF12[0]+df_TOF13[0])/4).mean()

                #still removing the mean delay between the detector and the reference detector
                print("\nThe required timing calibration is: ", WindowTimeOffset, '\n')

                config_args = {
                    "threshold": config["Thresholds"][i],
                    "analysis_window": (config["AnalysisWindowLow"][i], config["AnalysisWindowHigh"][i]),
                    "pedestal_window": (config["PedestalWindowLow"][i], config["PedestalWindowHigh"][i]),
                    "reverse_polarity": (config["Polarity"][i] == 0),
                    "voltage_scale": config["VoltageScale"],
                    "time_offset": config["TimeOffset"][i],
                    "window_time_offset": WindowTimeOffset,
                    "PMTgain": config["PMTgain"][i],
                }
                waveforms = batch[c]
                optional_branches = {b: batch[b] for b in optional_branch_names}

                process_batch(waveforms, optional_branches, config["ChannelNames"][i], output_file, config_args, df_TOF10.loc[int(report.tree_entry_start):int(report.tree_entry_stop)-1], window_lower_bound, window_upper_bound, df.loc[int(report.tree_entry_start):int(report.tree_entry_stop)-1])


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


def process_batch(waveforms, optional_branches, channel, output_file, config_args, df_TOF10, window_lower_bound, window_upper_bound, df):

    waveform_analysis = WaveformAnalysis(waveforms, window_lower_bound, window_upper_bound, **config_args)

    #instead of overwriting the variables, copy them from the root file
    #TODO: impose coincidence and only perform the integration on peaks that have been matched
    waveform_analysis.run_window_analysis(df_TOF10)
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
    window_peaks = ak.zip({"WindowIntCharge": window_integrated_charges,
                           "WindowIntPE": window_integrated_pe})

    branches = {"Pedestal": read_nTuple(df, "Pedestal", isSqueezed = True),
                "PedestalSigma": read_nTuple(df, "PedestalSigma", isSqueezed = True),
                "MaxVoltage": read_nTuple(df, "MaxVoltage", isSqueezed = True),
                "WholeWaveformInt": read_nTuple(df, "WholeWaveformInt", isSqueezed = True),
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
