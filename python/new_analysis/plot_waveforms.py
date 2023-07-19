#!/usr/bin/env python3
import uproot as ur
import numpy as np
import matplotlib.pyplot as plt
import json5
from waveform_analysis import WaveformAnalysis

run_file = ur.open("run_file.root")
config = json5.load(open("../../config/config.json"))['WaveAnalysis']
waveforms = run_file['midas_data_D300']['Channel0'].array().to_numpy()
analysis = WaveformAnalysis(waveforms,
                            threshold=config["Thresholds"][0],
                            analysis_window=(config["AnalysisWindowLow"][0], config["AnalysisWindowHigh"][0]),
                            pedestal_window=(config["PedestalWindowLow"][0], config["PedestalWindowHigh"][0]),
                            reverse_polarity=(config["Polarity"][0]==0),
                            voltage_scale=config["VoltageScale"],
                            time_offset=config["TimeOffset"][0])
n_peaks = analysis.count_peaks()

has_1_peak = np.where(n_peaks == 1)
has_2_peaks = np.where(n_peaks == 2)
has_more_peaks = np.where(n_peaks > 2)

rows = 5
fig, axs = plt.subplots(rows, 3, figsize=(12,2*rows))
fig.tight_layout()
axs[0][0].set_title("1 peak")
axs[0][1].set_title("2 peaks")
axs[0][2].set_title("3+ peaks")
for i in range(rows):
    axs[i][0].plot(waveforms[has_1_peak][i, analysis.analysis_bins])
    axs[i][1].plot(waveforms[has_2_peaks][i, analysis.analysis_bins])
    axs[i][2].plot(waveforms[has_more_peaks][i, analysis.analysis_bins])
plt.show()
