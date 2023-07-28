from Run import Run
from dk_utilities import get_signal_times

import json5
import numpy as np

# Perform standard analysis on 2023 data without hodoscope

run_number = 464
n_channels = 19
redo_pedestals = False

config = json5.load(open("../../config/config.json"))['WaveAnalysis']
config["NumberOfChannels"] = n_channels

config["peak_find_thresholds"] = [500]*2 + [250]*6 + [500]*8 + [250]*3
config["peak_min_amplitudes"] = [1000]*2 + [500]*6 + [1000]*8 + [500]*3

run = Run("C:/Users/Karlen/Documents/temp/root_run_000"+str(run_number)+".root", config)

# include additional information

run.user['run_number'] = run_number
run.user['n_channels'] = n_channels

# do pedestal calc (useful for special random trigger runs: eg. 503)

if run_number in [503]:

    baseline_means = []
    baseline_stds = []
    for analysis in run.analyses:
        stds = np.std(analysis.raw_waveforms, axis=1) # std by event
        means = np.mean(analysis.raw_waveforms, axis=1) # mean by event
        low_std = np.percentile(stds,10.)
        high_std = np.percentile(stds,90.)
        condition = (stds < high_std) & (stds > low_std)
        selected_means = means[condition]
        selected_stds = stds[condition]
        baseline_means.append(np.mean(selected_means))
        baseline_stds.append(np.mean(selected_stds))

    run.user['baseline_means'] = baseline_means
    run.user['baseline_stds'] = baseline_stds

# do multiple peak finding analysis to get signal times

run.user['signal_times'] = get_signal_times(run, config)

run.save_file("C:/Users/Karlen/Documents/temp/run_000"+str(run_number)+".dk")