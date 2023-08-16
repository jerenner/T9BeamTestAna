from Run import Run
import numpy as np
from scipy.optimize import curve_fit


def get_signal_times(run, config):

    def get_time_bin_crossing(start, n_samples, waveform):
        ts = np.array([start + i for i in range(n_samples)])
        amps = np.array([waveform[its] for its in ts])
        z = np.polyfit(ts, amps, 2)
        a, b, c = z[0], z[1], z[2]

        # fine binning for plot and for linear interpolation
        fine = 20
        tsf = np.array([start + i / (1. * fine) for i in range(fine * n_samples)])
        quad_fit = a * tsf * tsf + b * tsf + c

        deviations = quad_fit - baseline_means[channel] + threshold
        crossed_thresholds = np.where(deviations < 0.)[0]
        if len(crossed_thresholds) > 0:
            crossed_threshold = crossed_thresholds[0]
            slope = (quad_fit[crossed_threshold] - quad_fit[crossed_threshold - 1]) * fine
            time_bin_crossing = start + crossed_threshold / fine - deviations[crossed_threshold] / slope
        else:
            time_bin_crossing = -999.

        return time_bin_crossing, tsf, quad_fit

    # get baseline means and standard deviations
    if run.user['n_channels'] == 19:
        run_baseline = Run.open_file("C:/Users/Karlen/Documents/temp/run_000503.dk")
    else:
        run_baseline = Run.open_file("C:/Users/Karlen/Documents/temp/run_000592.dk")
    baseline_means = run_baseline.user['baseline_means']
    baseline_stds = run_baseline.user['baseline_stds']

    signal_times = []  # channels, events, peaks
    for channel in range(run.user['n_channels']):
        threshold = config["peak_find_thresholds"][channel]
        min_amplitude = config["peak_min_amplitudes"][channel]

        analysis = run.analyses[channel]

        signal_times_for_channel = []
        for iw, waveform in enumerate(analysis.raw_waveforms):

            # loop over peaks: must return to within 3 sigma before another peak is considered
            ready = False
            exceeded_3sigma = False
            peak_ranges = []

            baseline_mean = baseline_means[channel]
            baseline_std = baseline_stds[channel]
            for t_bin in range(len(waveform)):
                if ready:
                    exceeded_3sigma = (waveform[t_bin] - baseline_mean) < -3 * baseline_std
                    if exceeded_3sigma:
                        time_bin_exceed_3sigma = t_bin - 1
                        ready = False
                if exceeded_3sigma:
                    if (waveform[t_bin] - baseline_mean) < -1 * min_amplitude:
                        peak_ranges.append([time_bin_exceed_3sigma, t_bin])
                        exceeded_3sigma = False

                ready = (waveform[t_bin] - baseline_mean) > -3 * baseline_std

            time_bin_crossings = []
            for peak_range in peak_ranges:
                start = peak_range[0]
                end = peak_range[1]
                n_samples = max(4, end - start + 1)
                if start+n_samples < len(waveform):
                    time_bin_crossing, tsf, quad_fit = get_time_bin_crossing(start, n_samples, waveform)
                    if time_bin_crossing > 0:
                        time_bin_crossings.append(time_bin_crossing)

            signal_times_for_channel.append(time_bin_crossings)

        signal_times.append(signal_times_for_channel)

    return signal_times
