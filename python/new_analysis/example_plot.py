#!/usr/bin/env python3
import uproot as ur
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# read channels from file
file = ur.open("ntuple_file.root")
tof_0 = [file[c].arrays() for c in ["TOF00", "TOF01", "TOF02", "TOF03"]]
tof_1 = [file[c].arrays() for c in ["TOF10", "TOF11", "TOF12", "TOF13"]]
act_2_and_3 = [file[c].arrays() for c in ["ACT2L", "ACT2R", "ACT3L", "ACT3R"]]
all_channels = tof_0+tof_1+act_2_and_3

# find events where all ACT2+3 and all TOF have signal over threshold, and number of peaks is 1 or more in all these waveforms
single_pulses = np.where(np.all([c["nPeaks"] == 1 for c in all_channels], axis=0))[0]
is_over_threshold = np.all([c["PeakVoltage"][single_pulses, 0] > 0.2 for c in all_channels], axis=0)
good_events = single_pulses[is_over_threshold]

# calculate the mean TOF times for each event
tof0_times = np.mean([c["SignalTime"][good_events, 0] for c in tof_0], axis=0)
tof1_times = np.mean([c["SignalTime"][good_events, 0] for c in tof_1], axis=0)
# calculate the sum of ACT2+ACT3 amplitudes for each event
act23_sum = np.sum([c["PeakVoltage"][good_events, 0] for c in act_2_and_3], axis=0)

# plot 2d and 1d histograms
fig, axs = plt.subplots(2,2, constrained_layout=True, figsize=(12, 10))
axs[0, 1].axis('off')
tof_lim = (0, 30)
# 1d histogram of TOF
axs[0, 0].hist(tof1_times-tof0_times, bins=200, range=tof_lim, histtype="step")
axs[0, 0].set_xlim(*tof_lim)
axs[0, 0].set_yscale('log')
# 1d histogram of ACT2+ACT3 amplitude
act_lim = (0, 25)
axs[1, 1].hist(act23_sum, bins=200, range=act_lim, histtype="step", orientation='horizontal')
axs[1, 1].set_ylim(*act_lim)
# 2d histogram
h = axs[1, 0].hist2d(tof1_times-tof0_times, act23_sum, bins=200, range=(tof_lim, act_lim), norm=mcolors.LogNorm())
axs[1, 0].set_xlabel("TOF")
axs[1, 0].set_ylabel("ACT2+ACT3")
fig.colorbar(h[3], ax=axs[1, 0], pad=0)
plt.show()
