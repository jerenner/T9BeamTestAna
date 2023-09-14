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

    df_TOF10 = ref_file['TOF10'].arrays(library="pd")
    #add the hit times in single entries - easier, we are so far only looking at one PMT but could do the
    #average of them in the TOF1 module
    hit_times = pd.DataFrame(df_TOF10['SignalTime'].values.tolist())
    df_TOF10 =  pd.concat([df_TOF10, hit_times], axis=1)
    # #this is the first hit time
    # df_TOF10["hit0"] = df_TOF10[0]


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
    run_number = int(root_filename[-11:-5])
    for d in active_digitizers:
        optional_branch_names = [b for b in ["timeStamp", "triggerTime"] if b in run_file[d].keys()]
        channels = digitizer_channels[d]
        branch_names = list(channels.keys()) + optional_branch_names
        print(f"Reading digitizer {d}")
        for batch, report in run_file[d].iterate(branch_names, step_size="1000 MB", report=True):
            print(f"... events {report.tree_entry_start} to {report.tree_entry_stop}")
            for c, i in channels.items():
                print(f"... ... processing {d}/{c} into {config['ChannelNames'][i]}", end="", flush=True)
                config_args = {
                    "threshold": config["Thresholds"][i],
                    "analysis_window": (config["AnalysisWindowLow"][i], config["AnalysisWindowHigh"][i]),
                    "pedestal_window": (config["PedestalWindowLow"][i], config["PedestalWindowHigh"][i]),
                    "reverse_polarity": (config["Polarity"][i] == 0),
                    "voltage_scale": config["VoltageScale"],
                    "time_offset": config["TimeOffset"][i],
                    "window_time_offset": config["WindowTimeOffset"][i],
                }
                waveforms = batch[c]
                optional_branches = {b: batch[b] for b in optional_branch_names}
                print(df_TOF10, df_TOF10[0], c, i)
                process_batch(waveforms, optional_branches, config["ChannelNames"][i], output_file, config_args, df_TOF10, window_lower_bound, window_upper_bound)
                #save the waveforms pictures
                # plt.figure(2)
                plt.title(f"{config['ChannelNames'][i]} - Window ({window_lower_bound}_{window_upper_bound})")
                plt.savefig(f"pdf_results/sample_waveforms_40nsWindow502_{config['ChannelNames'][i]}_run{run_number}_{window_lower_bound}_{window_upper_bound}.pdf")


    print(f"Saving event info for run {run_number}")
    output_file["EventInfo"] = {
        "RunNumber": np.repeat(run_number, n_events),
        "EventNumber": run_file[digitizers[0]]["eventNumber"].array(),
        "SpillNumber": run_file[digitizers[0]]["spillNumber"].array()
    }
    run_file.close()


def process_batch(waveforms, optional_branches, channel, output_file, config_args, df_TOF10, window_lower_bound, window_upper_bound):
    waveform_analysis = WaveformAnalysis(waveforms, df_TOF10, window_lower_bound, window_upper_bound, **config_args)
    waveform_analysis.run_analysis()
    pedestals = ak.Array(waveform_analysis.pedestals.squeeze())
    pedestal_sigmas = ak.Array(waveform_analysis.pedestal_sigmas.squeeze())
    peak_voltages = waveform_analysis.pulse_peak_voltages
    peak_times = waveform_analysis.pulse_peak_times
    signal_times = waveform_analysis.pulse_signal_times
    integrated_charges = waveform_analysis.pulse_charges
    window_integrated_charges = waveform_analysis.window_int_charge

    print(type(integrated_charges), type(window_integrated_charges))

    peaks = ak.zip({"PeakVoltage": peak_voltages,
                    "PeakTime": peak_times,
                    "SignalTime": signal_times,
                    "IntCharge": integrated_charges})

    #needs to be handled separately because it has another format
    window_peaks = ak.zip({"WindowIntCharge": window_integrated_charges})

    branches = {"Pedestal": pedestals,
                "PedestalSigma": pedestal_sigmas,
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
#     run_file.close()
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
