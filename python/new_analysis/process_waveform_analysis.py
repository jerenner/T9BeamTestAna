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
    #monitor pedestal fluctuations
    runNb = int(root_filename[14:20])
    pedestal_list = [runNb]
    print(f"Processing root file {root_filename}")
    run_file = ur.open(root_filename)
    digitizers = ["midas_data_D300", "midas_data_D301", "midas_data_D302", "midas_data_D303"]
    channels_per_digitizer = config["NumberOfChannels"] // len(digitizers)
    channels = []
    for c in range(config["NumberOfChannels"]):
        channel = f"{digitizers[c // channels_per_digitizer]}/Channel{c % channels_per_digitizer}"
        if config["ActiveChannels"][c]:
            channels.append(f"{channel}")
        print(f"{channel} is {'active' if config['ActiveChannels'][c] else 'not active'}")
    events_per_channel = [run_file[c].num_entries for c in channels]
    if min(events_per_channel) != max(events_per_channel):
        raise ValueError(f"Channels have different numbers of events: {','.join([f'{c}:{n}' for c, n in zip(channels, events_per_channel)])}")
    n_channels = len(channels)
    n_events = events_per_channel[0]
    if n_events == 0:
        print(f"WARNING: The channels in {root_filename} have zero events. Skipping this file.")
        return
    print(f"All {n_channels} active channels have {n_events} events")
    for i, c in enumerate(channels):
        print(f"Processing {c}")
        waveforms = run_file[c].array()
        waveform_analysis = WaveformAnalysis(waveforms,
                                             threshold=config["Thresholds"][i],
                                             analysis_window=(config["AnalysisWindowLow"][i], config["AnalysisWindowHigh"][i]),
                                             pedestal_window=(config["PedestalWindowLow"][i], config["PedestalWindowHigh"][i]),
                                             reverse_polarity=(config["Polarity"][i] == 0),
                                             voltage_scale=config["VoltageScale"],
                                             time_offset=config["TimeOffset"][i])
        waveform_analysis.run_analysis()
        pedestals[:, i] = waveform_analysis.pedestals.squeeze()
        pedestal_sigmas[:, i] = waveform_analysis.pedestal_sigmas.squeeze()
        peak_voltages[:, i] = waveform_analysis.peak_voltages
        peak_times[:, i] = waveform_analysis.peak_times
        signal_times[:, i] = waveform_analysis.signal_times
        integrated_charges[:, i] = waveform_analysis.integrated_charges
        n_peaks[:, i] = waveform_analysis.n_peaks
        #print("pedestals", waveform_analysis.pedestals.squeeze(), waveform_analysis.my_pedestals)
        pedestal_list.append(waveform_analysis.my_pedestals)
        pedestal_list.append(waveform_analysis.my_pedestal_sigmas)
        print("Pedestals: ", pedestal_list)
        channel_name = cm.channel_names[i]
        print(f"Processing {c} into {channel_name}")
        for batch, report in run_file[c].iterate(step_size="100 MB", report=True):
            print(f"... events {report.tree_entry_start} to {report.tree_entry_stop}")
            config_args = {
                "threshold": config["Thresholds"][i],
                "analysis_window": (config["AnalysisWindowLow"][i], config["AnalysisWindowHigh"][i]),
                "pedestal_window": (config["PedestalWindowLow"][i], config["PedestalWindowHigh"][i]),
                "reverse_polarity": (config["Polarity"][i] == 0),
                "voltage_scale": config["VoltageScale"],
                "time_offset": config["TimeOffset"][i],
            }
            process_batch(batch[batch.fields[0]], channel_name, output_file, config_args)

    run_number = int(root_filename[-11:-5])
    print(f"Saving event info for run {run_number}")
    output_file["EventInfo"] = {
        "RunNumber": np.repeat(run_number, n_events),
        "EventNumber": run_file[digitizers[0]]["eventNumber"].array(),
        "SpillNumber": run_file[digitizers[0]]["spillNumber"].array()
    }


def process_batch(waveforms, channel, output_file, config_args):
    waveform_analysis = WaveformAnalysis(waveforms, **config_args)
    waveform_analysis.run_analysis()
    pedestals = waveform_analysis.pedestals.squeeze()
    pedestal_sigmas = waveform_analysis.pedestal_sigmas.squeeze()
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
                "Peaks": peaks}

    print(f"... writing {channel} to {output_file.file_path}")
    if channel not in output_file.keys():
        branch_types = {"Pedestal": "float64", "PedestalSigma": "float64", "Peaks": peaks.type}
        output_file.mktree(channel, branch_types, counter_name=(lambda s: "nPeaks"), field_name=(lambda s, f: f))
    output_file[channel].extend(branches)

    if True:
        #chekcing the stability of pedestal
        import csv
        with open(r'pedestalStabilityChecks.csv', 'a') as f:
            writer = csv.writer(f)
            writer.writerow(pedestal_list)

if __name__ == "__main__":
    # execute only if run as a script"
    process_waveforms(sys.argv)
