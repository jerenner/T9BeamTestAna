#!/usr/bin/python
#
# Process waveform analysis of root file
#

import sys, warnings

import numpy as np
import uproot as ur
import json5
from waveform_analysis import WaveformAnalysis


def PrintUsage(argv):
    print('Usage:')
    print(f"{argv[0]} input_file config_file output_file")
    return


def process_waveforms(argv):
    if len(argv) < 4:
        PrintUsage(argv)
        return

    root_filename = argv[1]
    config_filename = argv[2]
    output_filename = argv[3]
    print(f"Processing root file {root_filename} using config file {config_filename}")
    with open(config_filename) as file:
        config = json5.load(file)["WaveAnalysis"]
    run_file = ur.open(root_filename)
    digitizers = ["midas_data_D300", "midas_data_D301", "midas_data_D302", "midas_data_D303"]
    channels_per_digitizer = config["NumberOfChannels"]//len(digitizers)
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
    print(f"All {n_channels} active channels have {n_events} events")
    pedestals = np.zeros((n_events, n_channels))
    pedestal_sigmas = np.zeros((n_events, n_channels))
    n_peaks = 1
    peak_voltages = np.zeros((n_events, n_channels, n_peaks))
    peak_times = np.zeros((n_events, n_channels, n_peaks))
    signal_times = np.zeros((n_events, n_channels, n_peaks))
    integrated_charges = np.zeros((n_events, n_channels, n_peaks))

    for i, c in enumerate(channels):
        print(f"Processing {c}")
        waveforms = run_file[c].array()
        waveform_analysis = WaveformAnalysis(waveforms,
                                             threshold=config["Thresholds"][i],
                                             analysis_window=(config["AnalysisWindowLow"][i], config["AnalysisWindowHigh"][i]),
                                             pedestal_window=(config["PedestalWindowLow"][i], config["PedestalWindowHigh"][i]),
                                             reverse_polarity=(config["Polarity"][i] == 0),
                                             voltage_scale=config["VoltageScale"])
        waveform_analysis.run_analysis()
        pedestals[:, i] = waveform_analysis.pedestals.squeeze()
        pedestal_sigmas[:, i] = waveform_analysis.pedestal_sigmas.squeeze()
        peak_voltages[:, i] = waveform_analysis.peak_voltages
        peak_times[:, i] = waveform_analysis.peak_times
        signal_times[:, i] = waveform_analysis.signal_times
        integrated_charges[:, i] = waveform_analysis.integrated_charges

    print(f"Writing to {output_filename}")
    output_file = ur.recreate(output_filename)
    output_file["anaTree"] = {
        "nChannels": np.repeat(n_channels, n_events),
        "Pedestal": pedestals,
        "PedestalSigma": pedestal_sigmas,
        "PeakVoltage": peak_voltages,
        "PeakTime": peak_times,
        "SignalTime": signal_times,
        "IntCharge": integrated_charges
    }
    output_file.close()


if __name__ == "__main__":
    # execute only if run as a script"
    process_waveforms(sys.argv)
