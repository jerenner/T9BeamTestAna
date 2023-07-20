#!/usr/bin/env python3
import uproot as ur
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import channel_map as cm

# Read signal times and amplitudes from the run ntuple file
file = ur.open("ntuple_file.root")
signal_times = file["anaTree"]["SignalTime"].array().to_numpy()
signal_charges = file["anaTree"]["PeakVoltage"].array().to_numpy()
n_peaks = file["anaTree"]["nPeaks"].array().to_numpy()

# find events where all ACT2+3 and all TOF have signal over threshold, and number of peaks is 1 in all these waveforms
over_threshold = signal_charges[:, :, 0] > 0.2
good_events = (np.all(over_threshold[:, cm.ACT2and3], axis=1) &
               np.all(over_threshold[:, cm.TOF0and1], axis=1) &
               np.all(n_peaks[:, cm.ACT2and3]==1, axis=1) &
               np.all(n_peaks[:, cm.ACT2and3]==1, axis=1))

# calculate the mean TOF times for each event
tof0_times = signal_times[good_events, cm.TOF0, 0].mean(axis=1)
tof1_times = signal_times[good_events, cm.TOF1, 0].mean(axis=1)
# calculate the sum of ACT2+ACT3 amplitudes for each event
act23_sum = signal_charges[good_events, cm.ACT2and3, 0].sum(axis=1)

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
