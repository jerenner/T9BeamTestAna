import numpy as np
from scipy.signal import find_peaks


def find_mode(values, n_bins, range_width):
    """Finds the mode of a distribution using a histogram centered around the median"""
    hist, bin_edges = np.histogram(values, bins=n_bins, range=np.median(values) + [-range_width, range_width])
    max_bin = np.argmax(hist)
    return np.mean(bin_edges[max_bin:max_bin + 2])

class WaveformAnalysis:
    """
    This class holds an array of waveforms, for quickly doing waveform analysis on all waveforms
    """
    peak_rise_time_fraction = 0.4  # signal time is when before the peak the waveform passes this fraction of the peak amplitude
    impedance = 50  # for calculating integrated charge
    use_global_pedestal = True  # assume pedestal is stable and use one value for all waveforms for the channel

    def __init__(self, waveforms, threshold=0.01, analysis_window=(0, 200), pedestal_window=(200, 420),
                 reverse_polarity=True, ns_per_sample=2, voltage_scale=0.000610351, time_offset=0):

        self.waveforms = np.array(waveforms)
        self.threshold = threshold
        self.reverse_polarity = reverse_polarity
        self.ns_per_sample = ns_per_sample
        self.voltage_scale = voltage_scale
        self.analysis_window = analysis_window
        self.pedestal_window = pedestal_window
        self.analysis_bins = range(analysis_window[0]//ns_per_sample, analysis_window[1]//ns_per_sample)
        self.pedestal_bins = range(pedestal_window[0]//ns_per_sample, pedestal_window[1]//ns_per_sample)
        self.time_offset = time_offset
        self.peak_locations = None
        self.amplitudes = None

    def find_pedestals(self):
        """Finds the pedestal of each waveform by taking the mean in the pedestal window. Also finds the standard deviation"""
        self.amplitudes = self.waveforms * self.voltage_scale
        self.pedestals = np.mean(self.amplitudes[:, self.pedestal_bins], axis=1, keepdims=True)
        self.pedestal_sigmas = np.std(self.amplitudes[:, self.pedestal_bins], axis=1, keepdims=True)
        if self.use_global_pedestal:  # use the pedestal distribution's mode to ignore outliers from unusual waveforms
            self.my_pedestals = find_mode(self.pedestals, n_bins=100, range_width=100*self.voltage_scale)
            self.my_pedestal_sigmas = find_mode(self.pedestal_sigmas, n_bins=100, range_width=100*self.voltage_scale)
        else:
            self.my_pedestals = self.pedestals
            self.my_pedestal_sigmas = self.pedestal_sigmas
        self.amplitudes -= self.my_pedestals
        if self.reverse_polarity:
            self.amplitudes *= -1

    def find_peaks(self):
        """Finds the peak voltage in the analysis window of each waveform"""
        if self.amplitudes is None:
            self.find_pedestals()
        self.peak_locations = np.argmax(self.amplitudes[:, self.analysis_bins], axis=1, keepdims=True) + self.analysis_bins[0]
        self.peak_times = (self.peak_locations + 0.5)*self.ns_per_sample
        self.peak_voltages = np.take_along_axis(self.amplitudes, self.peak_locations, axis=1).reshape(self.peak_locations.shape)

    def calculate_signal_times(self):
        """Finds the signal time of each waveform as the interpolated time before the peak when the voltage passes 0.4*[peak voltage]"""
        if self.peak_locations is None:
            self.find_peaks()
        thresholds = self.peak_rise_time_fraction*self.peak_voltages
        below_threshold = self.amplitudes < thresholds
        indices = np.arange(self.amplitudes.shape[1])
        after_peak_rise = np.cumsum(((indices < self.peak_locations) & below_threshold)[:, ::-1], axis=1)[:, ::-1] == 0
        peak_rise_location = np.sum(after_peak_rise == 0, axis=1, keepdims=True)  # count how many samples are before the peak rise
        amp_range = np.take_along_axis(self.amplitudes, np.column_stack([peak_rise_location-1, peak_rise_location]), axis=1)
        self.signal_times = peak_rise_location - (amp_range[:, [1]] - thresholds)/np.diff(amp_range)  # interpolate
        self.signal_times = self.signal_times*self.ns_per_sample + self.time_offset

    def integrate_charges(self):
        """Calculates the integrated charge of each waveform's range around the peak that is more than 3x the standard deviation of the pedestal"""
        if self.peak_locations is None:
            self.find_peaks()
        below_threshold = self.amplitudes < (3*self.my_pedestal_sigmas)
        indices = np.arange(self.amplitudes.shape[1])
        before_peak_fall = np.cumsum((indices > self.peak_locations) & below_threshold, axis=1) == 0
        after_peak_rise = np.cumsum(((indices < self.peak_locations) & below_threshold)[:, ::-1], axis=1)[:, ::-1] == 0
        self.integrated_charges = np.sum(self.amplitudes*(before_peak_fall & after_peak_rise), axis=1, keepdims=True)
        self.integrated_charges *= self.ns_per_sample / self.impedance

    def run_analysis(self):
        """Finds each waveform peak, calculates the time, integrates the charge, and checks if it is over threshold"""
        self.find_peaks()
        self.calculate_signal_times()
        self.integrate_charges()
        self.is_over_threshold = self.peak_voltages > self.threshold
