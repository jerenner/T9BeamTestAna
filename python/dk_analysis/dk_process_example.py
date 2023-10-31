from Run import Run
from dk_utilities import get_signal_times

import json5
import numpy as np

# Perform standard analysis on 2023 data without/with hodoscope

#folder = "C:/Users/Karlen/Documents/temp/"
folder = "C:/Users/Karlen/Documents/temp/wtime/"

#run_numbers = [628,629,630,631,632,633,634,635,636,637,638,639,640,641,642,644]
#run_numbers = [592,635]
run_numbers = [731,740,745,753]
for run_number in run_numbers:
    n_channels = 32
    if run_number <= 578:
        n_channels = 19

    #config = json5.load(open("../../config/config.json"))['WaveAnalysis']
    config = json5.load(open("../../config/config_noHodoscope.json"))['WaveAnalysis']
    config["NumberOfChannels"] = n_channels

    if n_channels == 19:
        config["peak_find_thresholds"] = [500] * 2 + [250] * 6 + [500] * 8 + [250] * 3
        config["peak_min_amplitudes"] = [1000] * 2 + [500] * 6 + [1000] * 8 + [500] * 3
    else:
        config["peak_find_thresholds"] = [500]*2 + [250] + [500] + [250]*2 + [125]*2 + [250]*8 + [500]*2 + [250]*14
        config["peak_min_amplitudes"] = [1000]*2 + [500] + [1000] + [500]*2 + [250]*2 +[500]*8 +[1000]*2 + [500]*14

    run = Run(folder+"root_run_000"+str(run_number)+".root", config)

    # include additional information

    run.user['run_number'] = run_number
    run.user['n_channels'] = n_channels

    #

    # do pedestal calc (useful for special random trigger runs: eg. 503,592)

    # NOTE: run 592: the baseline appears to be incorrect for the T2 (channel 6).
    # Need to take another random trigger run (without beam): Run 676
    # That did not solve the issue. Worse estimate for channel 6 and 4.

    if run_number in [503,592,676]:

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

        if run_number == 592:
            # for some reason this value is offset from beam data
            baseline_means[6] = 15500

        run.user['baseline_means'] = baseline_means
        run.user['baseline_stds'] = baseline_stds

    # do multiple peak finding analysis to get signal times

    run.user['signal_times'] = get_signal_times(run, config)

    run.save_file(folder+"run_000"+str(run_number)+".dk")