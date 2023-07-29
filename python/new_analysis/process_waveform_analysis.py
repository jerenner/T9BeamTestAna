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
from waveform_analysis import WaveformAnalysis
import channel_map as cm
import awkward as ak


def print_usage(argv):
    print('Usage:')
    print(f"{argv[0]} input_file [input_file2 ...] config_file output_file")
    print(f"    input_file, input_file2, ...: any number of root files from midas2root")
    print(f"    config_file: config JSON file")
    print(f"    output_file: output ntuple root file")
    return


def process_waveforms(argv):
    if len(argv) < 4:
        print("Not enough arguments provided.")
        print_usage(argv)
        return

    root_filenames = argv[1:-2]
    config_filename = argv[-2]
    output_filename = argv[-1]

    print(f"Using config file {config_filename}")
    with open(config_filename) as file:
        config = json5.load(file)["WaveAnalysis"]

    output_file = ur.recreate(output_filename)
    print(f"{len(root_filenames)} input root files to process.")
    for f in root_filenames:
        process_file(f, config, output_file)
    output_file.close()


def process_file(root_filename, config, output_file):
    print(f"Processing root file {root_filename}")
    run_file = ur.open(root_filename)
    digitizers = ["midas_data_D300", "midas_data_D301", "midas_data_D302", "midas_data_D303"]
    channels_per_digitizer = config["NumberOfChannels"] // len(digitizers)
    digitizer_channels = {d: [] for d in digitizers}
    for c in range(config["NumberOfChannels"]):
        channel = f"Channel{c % channels_per_digitizer}"
        digitizer = digitizers[c // channels_per_digitizer]
        if config["ActiveChannels"][c]:
            digitizer_channels[digitizer].append(f"{channel}")
        print(f"{digitizer}/{channel} is {'active' if config['ActiveChannels'][c] else 'not active'}")
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
    i = 0
    optional_branches = ["timeStamp", "triggerTime"]
    for i, (d, c) in enumerate((d, c) for d, channels in digitizer_channels.items() for c in channels):
        channel_name = cm.channel_names[i]
        print(f"Processing {d}/{c} into {channel_name}")
        branches = [c] + [b for b in optional_branches if b in run_file[d].keys()]
        for batch, report in run_file[d].iterate(branches, step_size="100 MB", report=True):
            print(f"... events {report.tree_entry_start} to {report.tree_entry_stop}")
            config_args = {
                "threshold": config["Thresholds"][i],
                "analysis_window": (config["AnalysisWindowLow"][i], config["AnalysisWindowHigh"][i]),
                "pedestal_window": (config["PedestalWindowLow"][i], config["PedestalWindowHigh"][i]),
                "reverse_polarity": (config["Polarity"][i] == 0),
                "voltage_scale": config["VoltageScale"],
                "time_offset": config["TimeOffset"][i],
            }
            waveforms = batch[c]
            extra_branches = {b: batch[b] for b in branches if b != c}
            process_batch(waveforms, extra_branches, channel_name, output_file, config_args)

    run_number = int(root_filename[-11:-5])
    print(f"Saving event info for run {run_number}")
    output_file["EventInfo"] = {
        "RunNumber": np.repeat(run_number, n_events),
        "EventNumber": run_file[digitizers[0]]["eventNumber"].array(),
        "SpillNumber": run_file[digitizers[0]]["spillNumber"].array()
    }
    run_file.close()


def process_batch(waveforms, extra_branches, channel, output_file, config_args):
    waveform_analysis = WaveformAnalysis(waveforms, **config_args)
    waveform_analysis.run_analysis()
    pedestals = ak.Array(waveform_analysis.pedestals.squeeze())
    pedestal_sigmas = ak.Array(waveform_analysis.pedestal_sigmas.squeeze())
    peak_voltages = waveform_analysis.pulse_peak_voltages
    peak_times = waveform_analysis.pulse_peak_times
    signal_times = waveform_analysis.pulse_signal_times
    integrated_charges = waveform_analysis.pulse_charges

    peaks = ak.zip({"PeakVoltage": peak_voltages,
                    "PeakTime": peak_times,
                    "SignalTime": signal_times,
                    "IntCharge": integrated_charges})

    branches = {"Pedestal": pedestals,
                "PedestalSigma": pedestal_sigmas,
                "Peaks": peaks,
                **extra_branches}

    print(f"... writing {channel} to {output_file.file_path}")
    if channel not in output_file.keys():
        branch_types = {name: branch.type for name, branch in branches.items()}
        output_file.mktree(channel, branch_types, counter_name=(lambda s: "nPeaks"), field_name=(lambda s, f: f))
    output_file[channel].extend(branches)


if __name__ == "__main__":
    # execute only if run as a script"
    process_waveforms(sys.argv)
