import numpy as np
import scipy.signal as signal
import awkward as ak
import matplotlib.pyplot as plt


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
                 reverse_polarity=True, ns_per_sample=2, voltage_scale=0.000610351, time_offset=0, window_time_offset=0, PMTgain = 1, window_lower_bound = -15, window_upper_bound = 35, PMTdistanceToTOF1 = None, distanceTOF1toTOF0 = None, checkCoincidence = False):
        self.raw_waveforms = np.array(waveforms)
        self.threshold = threshold
        self.reverse_polarity = reverse_polarity
        self.ns_per_sample = ns_per_sample
        self.voltage_scale = voltage_scale
        self.analysis_window = analysis_window
        self.pedestal_window = pedestal_window
        self.PMTdistanceToTOF1 = PMTdistanceToTOF1
        self.distanceTOF1toTOF0 = distanceTOF1toTOF0
        self.checkCoincidence = checkCoincidence
        self.analysis_bins = range(analysis_window[0]//ns_per_sample, analysis_window[1]//ns_per_sample)
        self.pedestal_bins = range(pedestal_window[0]//ns_per_sample, pedestal_window[1]//ns_per_sample)
        self.time_offset = time_offset
        self.window_lower_bound = window_lower_bound
        self.window_upper_bound = window_upper_bound
        #to aligne correctly the T0F1 and the other detectors
        self.window_time_offset = window_time_offset
        self.waveforms = None
        self.pedestals = None
        self.pedestal_sigmas = None
        self.my_pedestals = None
        self.my_pedestal_sigmas = None
        self.peak_locations = None
        self.peak_times = None
        self.peak_voltages = None
        self.smoothed_waveforms = None
        self.n_peaks = None
        self.signal_times = None
        self.integrated_charges = None
        self.is_over_threshold = None
        self.pulse_peak_voltages = None
        self.pulse_peak_times = None
        self.pulse_signal_times = None
        self.pulse_charges = None
        self.pulse_pe = None
        self.windowWidth = None
        #new - window integration
        self.window_int_charge = None
        self.window_int_pe = None
        #new save bounds of the integration window and ccentral value
        self.windowIntCentralTime = None
        self.windowIntCentralTimeCorrected = None
        self.windowIntUpperTime = None
        self.windowIntLowerTime = None
        #new - complete waveform integration for two particle ID
        self.whole_waveform_int = None
        self.whole_waveform_int_pe = None
        self.max_voltage = None
        self.PMTgain = PMTgain
        self.waveformEnd = None

        

    def find_pedestals(self):
        """Finds the pedestal of each waveform by taking the mean in the pedestal window. Also finds the standard deviation"""
        self.waveforms = self.raw_waveforms * self.voltage_scale
        self.pedestals = np.mean(self.waveforms[:, self.pedestal_bins], axis=1, keepdims=True)
        self.pedestal_sigmas = np.std(self.waveforms[:, self.pedestal_bins], axis=1, keepdims=True)
        if self.use_global_pedestal:  # use the pedestal distribution's mode to ignore outliers from unusual waveforms
            self.my_pedestals = find_mode(self.pedestals, n_bins=100, range_width=100*self.voltage_scale)
            self.my_pedestal_sigmas = find_mode(self.pedestal_sigmas, n_bins=100, range_width=100*self.voltage_scale)
        else:
            self.my_pedestals = self.pedestals
            self.my_pedestal_sigmas = self.pedestal_sigmas
        self.waveforms -= self.my_pedestals
        if self.reverse_polarity:
            self.waveforms *= -1
        return self.my_pedestals, self.my_pedestal_sigmas

    def find_peaks(self):
        self.find_primary_peaks()
        self.count_peaks()
        return self.n_peaks, self.peak_times, self.peak_voltages

    def find_primary_peaks(self):
        """Finds the peak voltage in the analysis window of each waveform"""
        if self.waveforms is None:
            self.find_pedestals()
        self.peak_locations = np.argmax(self.waveforms[:, self.analysis_bins], axis=1, keepdims=True) + self.analysis_bins[0]
        self.peak_times = (self.peak_locations + 0.5)*self.ns_per_sample + self.time_offset
        self.peak_voltages = np.take_along_axis(self.waveforms, self.peak_locations, axis=1).reshape(self.peak_locations.shape)

        self.max_voltage = np.take_along_axis(self.waveforms, self.peak_locations, axis=1).reshape(self.peak_locations.shape)
        return self.peak_times, self.peak_voltages

    def count_peaks(self):
        """Counts how many peaks passing the threshold fraction of the primary peak, and not overlapping in charge integration window"""
        if self.peak_voltages is None:
            self.find_primary_peaks()
        self.smoothed_waveforms = signal.savgol_filter(self.waveforms, 5, 2, mode='nearest')
        analysis_waveform = self.smoothed_waveforms[:, self.analysis_bins]
        over_integration_threshold = analysis_waveform > np.maximum(0.05*self.peak_voltages, 3*self.my_pedestal_sigmas)
        over_pulse_threshold = analysis_waveform > np.maximum(0.2*self.peak_voltages, 3*self.my_pedestal_sigmas)
        pass_integration_threshold = np.diff(over_integration_threshold, axis=1, prepend=0) == 1
        pass_pulse_threshold = np.diff(over_pulse_threshold, axis=1, prepend=0) == 1
        passed_integration_threshold = False
        self.n_peaks = np.zeros(self.waveforms.shape[0], dtype=int)
        for i in range(0, analysis_waveform.shape[1]):  # scan over waveform samples, quicker than iterating over waveforms
            passed_integration_threshold = passed_integration_threshold | pass_integration_threshold[:, i]  # check if passed integration threshold
            self.n_peaks += pass_pulse_threshold[:, i] & passed_integration_threshold  # there's a new peak when it passes pulse threshold and has passed integration threshold
            passed_integration_threshold &= ~pass_pulse_threshold[:, i]  # after passing pulse threshold, don't consider it passed integration threshold until it passes again
        return self.n_peaks

    def find_all_peak_voltages(self):
        if self.n_peaks is None:
            self.count_peaks()
        analysis_waveform = self.smoothed_waveforms[:, self.analysis_bins]
        over_integration_threshold = analysis_waveform > np.maximum(0.05*self.peak_voltages, 3*self.my_pedestal_sigmas)
        over_pulse_threshold = analysis_waveform > np.maximum(0.2*self.peak_voltages, 3*self.my_pedestal_sigmas)
        pass_integration_threshold = np.diff(over_integration_threshold, axis=1, prepend=0) == 1
        pass_pulse_threshold = np.diff(over_pulse_threshold, axis=1, prepend=0) == 1
        passed_integration_threshold = False
        pulse_peak_times = -np.ones((self.waveforms.shape[0], self.n_peaks.max()))
        pulse_peak_voltages = -np.ones((self.waveforms.shape[0], self.n_peaks.max()))
        my_pulse_voltages = -np.ones(self.waveforms.shape[0])


        n_peaks = -np.ones(self.waveforms.shape[0], dtype=int)
        for i in range(0, analysis_waveform.shape[1]):  # scan over waveform samples, quicker than iterating over waveforms
            passed_integration_threshold = passed_integration_threshold | pass_integration_threshold[:, i]  # check if passed integration threshold
            n_peaks += pass_pulse_threshold[:, i] & passed_integration_threshold  # there's a new peak when it passes pulse threshold and has passed integration threshold
            passed_integration_threshold &= ~pass_pulse_threshold[:, i]  # after passing pulse threshold, don't consider it passed integration threshold until it passes again
            pulse_peak_indices = np.where(over_pulse_threshold[:, i] & (analysis_waveform[:, i] > my_pulse_voltages))[0]
            my_pulse_voltages[pulse_peak_indices] = analysis_waveform[pulse_peak_indices, i]
            pulse_peak_voltages[pulse_peak_indices, n_peaks[pulse_peak_indices]] = my_pulse_voltages[pulse_peak_indices]
            pulse_peak_times[pulse_peak_indices, n_peaks[pulse_peak_indices]] = (i+0.5)*self.ns_per_sample + self.time_offset
            my_pulse_voltages[~over_integration_threshold[:, i]] = -1
        self.pulse_peak_times = ak.drop_none(np.ma.MaskedArray(pulse_peak_times, pulse_peak_times==-1))
        self.pulse_peak_voltages = ak.drop_none(np.ma.MaskedArray(pulse_peak_voltages, pulse_peak_voltages==-1))

    def calculate_all_signal_times(self):
        if self.pulse_peak_voltages is None:
            self.find_all_peak_voltages()
        analysis_waveform = self.smoothed_waveforms[:, self.analysis_bins]
        over_integration_threshold = analysis_waveform > np.maximum(0.05*self.peak_voltages, 3*self.my_pedestal_sigmas)
        over_pulse_threshold = analysis_waveform > np.maximum(0.2*self.peak_voltages, 3*self.my_pedestal_sigmas)
        pass_integration_threshold = np.diff(over_integration_threshold, axis=1, prepend=0) == 1
        pass_pulse_threshold = np.diff(over_pulse_threshold, axis=1, prepend=0) == 1
        passed_integration_threshold = False
        pulse_times = -np.ones((self.waveforms.shape[0], self.n_peaks.max()))
        my_pulse_times = -np.ones(self.waveforms.shape[0])
        pulse_time_thresholds = np.zeros(self.waveforms.shape[0])
        n_peaks = -np.ones(self.waveforms.shape[0], dtype=int)
        for i in range(0, analysis_waveform.shape[1]):  # scan over waveform samples, quicker than iterating over waveforms
            passed_integration_threshold = passed_integration_threshold | pass_integration_threshold[:, i]  # check if passed integration threshold
            new_peak_indices = np.where(pass_pulse_threshold[:, i] & passed_integration_threshold)[0]
            n_peaks[new_peak_indices] += 1  # there's a new peak when it passes pulse threshold and has passed integration threshold
            my_pulse_times[new_peak_indices] = -1
            passed_integration_threshold &= ~pass_pulse_threshold[:, i]  # after passing pulse threshold, don't consider it passed integration threshold until it passes again
            pulse_time_thresholds[new_peak_indices] = self.peak_rise_time_fraction*self.pulse_peak_voltages[new_peak_indices,n_peaks[new_peak_indices]]
            signal_time_indices = np.where((n_peaks>=0) & (analysis_waveform[:, i] > pulse_time_thresholds) & (my_pulse_times == -1))[0]
            my_pulse_times[signal_time_indices] = i + self.analysis_bins[0]
            if i > 0:
                signal_time_high = analysis_waveform[signal_time_indices, i]
                signal_time_low = analysis_waveform[signal_time_indices, i - 1]
                diff = signal_time_high - signal_time_low
                is_nonzero = diff != 0
                nonzero_indices = signal_time_indices[is_nonzero]
                my_pulse_times[nonzero_indices] -= (signal_time_high[is_nonzero] - pulse_time_thresholds[nonzero_indices]) / diff[is_nonzero]
            pulse_times[signal_time_indices, n_peaks[signal_time_indices]] = my_pulse_times[signal_time_indices]*self.ns_per_sample + self.time_offset
        self.pulse_signal_times = ak.drop_none(np.ma.MaskedArray(pulse_times, pulse_times == -1))

    def calculate_all_pulse_charges(self):
        if self.n_peaks is None:
            self.count_peaks()
        analysis_waveform = self.smoothed_waveforms[:, self.analysis_bins]

        over_integration_threshold = analysis_waveform > np.maximum(0.05*self.peak_voltages, 3*self.my_pedestal_sigmas)
        over_pulse_threshold = analysis_waveform > np.maximum(0.2*self.peak_voltages, 3*self.my_pedestal_sigmas)

        # print(self.peak_voltages.shape, self.my_pedestal_sigmas.shape, np.maximum(0.05*self.peak_voltages, 3*self.my_pedestal_sigmas), "The sum:", sum(np.where(np.maximum(0.05*self.peak_voltages, 3*self.my_pedestal_sigmas)!=3*self.my_pedestal_sigmas, 1, 0)), "\n \n")
        #
        # over_integration_threshold = analysis_waveform > np.maximum(0.05*self.peak_voltages, -9999 *3*self.my_pedestal_sigmas)
        # over_pulse_threshold = analysis_waveform > np.maximum(0.2*self.peak_voltages, 3*self.my_pedestal_sigmas)

        integration_threshold_diff = np.diff(over_integration_threshold, axis=1, prepend=0)
        pass_pulse_threshold = np.diff(over_pulse_threshold, axis=1, prepend=0) == 1
        passed_integration_threshold = False
        pulse_charges = np.zeros((self.waveforms.shape[0], self.n_peaks.max()))
        pulse_pe = np.zeros((self.waveforms.shape[0], self.n_peaks.max()))
        my_pulse_charges = np.zeros(self.waveforms.shape[0])
        n_peaks = -np.ones(self.waveforms.shape[0], dtype=int)
        for i in range(0, analysis_waveform.shape[1]):  # scan over waveform samples, quicker than iterating over waveforms
            passed_integration_threshold = passed_integration_threshold | (integration_threshold_diff[:, i] == 1)  # check if passed integration threshold
            n_peaks += pass_pulse_threshold[:, i] & passed_integration_threshold  # there's a new peak when it passes pulse threshold and has passed integration threshold
            passed_integration_threshold &= ~pass_pulse_threshold[:, i]  # after passing pulse threshold, don't consider it passed integration threshold until it passes again
            my_pulse_charges[integration_threshold_diff[:, i] == 1] = 0
            my_pulse_charges[over_integration_threshold[:, i]] += analysis_waveform[over_integration_threshold[:,i],i]
            end_of_pulse = np.where(((integration_threshold_diff[:, i] == -1) | (i == analysis_waveform.shape[1]-1)) & (~passed_integration_threshold))[0]
            pulse_charges[end_of_pulse, n_peaks[end_of_pulse]] = my_pulse_charges[end_of_pulse]*self.ns_per_sample/self.impedance
            pulse_pe[end_of_pulse, n_peaks[end_of_pulse]] = my_pulse_charges[end_of_pulse]/self.PMTgain
        self.pulse_charges = ak.drop_none(np.ma.MaskedArray(pulse_charges, pulse_charges==0))
        self.pulse_pe = ak.drop_none(np.ma.MaskedArray(pulse_pe, pulse_pe==0))
        # print(self.pulse_charges)


    def integrate_charge_in_window(self):
        """using a window of known position and width to integrate"""
        # self.smoothed_waveforms = signal.savgol_filter(self.waveforms, 5, 2, mode='nearest')
        #Accept the entire waveform, because of digitiser slipping  we can have spillage over into the pedestal region, but the only important thing is the reference timing
        #technically, peaks are only found in the analysis region so only a tiny number of peaks might
        #end up being outside fo that range, but it is the most robust way to do it
        analysis_waveform = self.smoothed_waveforms
        #npew the analysis is over all bins
        # print(self.analysis_bins)
        self.waveformEnd = len(self.smoothed_waveforms[0])
        allAnalysisBins = np.tile(np.arange(0, self.waveformEnd),(self.waveforms.shape[0],1))

        list_high = []
        list_low = []
        #select all of the windows of integration

        #maximum number of hits that were found
        nPeaksMax = max(self.df_TOF1_hitTimes["nPeaks"])

        pulse_charges = np.zeros((self.waveforms.shape[0],nPeaksMax))
        pulse_pe = np.zeros((self.waveforms.shape[0],nPeaksMax))
        windowWidth_array = np.zeros((self.waveforms.shape[0],nPeaksMax))

        #saving the timing characterisics of the waveform
        windowIntCentralTime = np.zeros((self.waveforms.shape[0],nPeaksMax))
        windowIntCentralTimeCorrected = np.zeros((self.waveforms.shape[0],nPeaksMax))
        windowIntLowerBound = np.zeros((self.waveforms.shape[0],nPeaksMax))
        windowIntUpperBound = np.zeros((self.waveforms.shape[0],nPeaksMax))



        #here, take into account the particle time of flight to predict where to 
        #centre the integration window, but only when we have calculated the coincidence
        # print("Particle time of flight: ",  self.particleTimeOfFlight[:10], " PMT distance to TOF1: ", self.PMTdistanceToTOF1, "previous offset: ",  self.window_time_offset)
        
        




        for i in range(nPeaksMax):
            over_integration_threshold = False
            #need to get rid of the NaN issues and look at bin id instead of ns

            #Look at the timing as signal time: need to shift the waveform and make sure we stay within the waveform boundary
            #DigiTimingOffset (DTO) is SignalTimeCorrected (STC, Arturo's corrections) - SignalTime(ST), waveforms are given in ST and df_TOF1_hitTimes given in STC
            #expectedMin/Max are given in ST so they can be applied to the waveform

            #we keep the absolute offset so we have an easy handle on that 
            # if self.checkCoincidence:
            #     self.window_time_offset = self.particleTimeOfFlight[i] * self.PMTdistanceToTOF1/self.distanceTOF1toTOF0 + self.window_time_offset 

            #     print("Sample of time of flight induced timing offset in position of window integration: ", np.array(self.window_time_offset[:5]))

            expectedMin = (self.df_TOF1_hitTimes[i]+self.window_lower_bound+self.window_time_offset - np.array(self.df_PMT["DigiTimingOffset"]))/self.ns_per_sample


            expectedMax = (self.df_TOF1_hitTimes[i]+self.window_upper_bound+self.window_time_offset - np.array(self.df_PMT["DigiTimingOffset"]))/self.ns_per_sample

            #replace the nans and stay within the boundary
            rangeLow = np.where(np.isnan(self.df_TOF1_hitTimes[i]), -9999, np.where(0>expectedMin, np.where(0<expectedMax, 0, -9999), expectedMin))

            #replace the nans and stay within the boundary
            rangeHigh = np.where(np.isnan(self.df_TOF1_hitTimes[i]), -9999, np.where(expectedMax>self.waveformEnd, np.where(expectedMin<self.waveformEnd, self.waveformEnd, -9999), expectedMax))
            # if debug:
            #     print("Check range Low", rangeLow, rangeLow[1888], rangeLow[1890], rangeLow[1891] )
            #     print("Check range High", rangeHigh, rangeHigh[1888], rangeHigh[1890], rangeHigh[1891] )

            windowWidth = rangeHigh-rangeLow 

            #need to keep this format before the reshaping to save the bound
            lowerBound = rangeLow
            upperBound = rangeHigh

            list_high.append(rangeHigh)
            list_low.append(rangeLow)

            rangeLow = np.array(rangeLow).reshape(len(rangeLow), 1)
            rangeHigh = np.array(rangeHigh).reshape(len(rangeHigh), 1)


            # print("Range Low- range high:", rangeHigh-rangeLow, rangeHigh, rangeLow, "offset:", self.df_PMT["DigiTimingOffset"])

            #creating a set of integration windows covering (Need to make this for every hit in TOF10
            #otherwise overlapping peaks will cause an issue.

            # print((allAnalysisBins > rangeLow), (allAnalysisBins < rangeHigh))
            over_integration_threshold = over_integration_threshold | ((allAnalysisBins > rangeLow) & (allAnalysisBins < rangeHigh))

            # print(over_integration_threshold)


            #purely copy-pasted, we do the same thing with a different integration window
            #we just changed over_pusle_threshold for over_integration_threshold
            integration_threshold_diff = np.diff(over_integration_threshold, axis=1, prepend=0)
            pass_pulse_threshold = np.diff(over_integration_threshold, axis=1, prepend=0) == 1
            passed_integration_threshold = False
            # n_peaks = -np.ones(self.waveforms.shape[0], dtype=int)
            my_pulse_charges = np.zeros(self.waveforms.shape[0])

            for k in range(0, analysis_waveform.shape[1]):  # scan over waveform samples, quicker than iterating over waveforms
                passed_integration_threshold = passed_integration_threshold | (integration_threshold_diff[:, k] == 1)  # check if passed integration threshold
# #                 n_peaks += pass_pulse_threshold[:, k] & passed_integration_threshold  # there's a new peak when it passes pulse threshold and has passed integration threshold
                passed_integration_threshold &= ~pass_pulse_threshold[:, k]  # after passing pulse threshold, don't consider it passed integration threshold until it passes again
                my_pulse_charges[integration_threshold_diff[:, k] == 1] = 0
                my_pulse_charges[over_integration_threshold[:, k]] += analysis_waveform[over_integration_threshold[:,k],k]
                end_of_pulse = np.where(((integration_threshold_diff[:, k] == -1) | (k == analysis_waveform.shape[1]-1)) & (~passed_integration_threshold))[0]

                pulse_charges[end_of_pulse, i] = my_pulse_charges[end_of_pulse]*self.ns_per_sample/self.impedance
                pulse_pe[end_of_pulse, i] = my_pulse_charges[end_of_pulse]/self.PMTgain
                windowWidth_array[:, i] = windowWidth*self.ns_per_sample

                #saving the edge of each of the windows
                windowIntLowerBound[:, i] = lowerBound*self.ns_per_sample

                windowIntUpperBound[:, i] = upperBound*self.ns_per_sample

                windowIntCentralTime[:, i] = self.df_TOF1_hitTimes[i]+self.window_time_offset - np.array(self.df_PMT["DigiTimingOffset"])

                windowIntCentralTimeCorrected[:,i] = self.df_TOF1_hitTimes[i]+self.window_time_offset

                # print(pulse_charges, end_of_pulse)

            self.window_int_charge = ak.drop_none(np.ma.MaskedArray(pulse_charges, pulse_charges==0))
            self.window_int_pe = ak.drop_none(np.ma.MaskedArray(pulse_pe, pulse_pe==0))
            #need to remove everything that was smaller than 1 data point
            self.windowWidth = ak.drop_none(np.ma.MaskedArray(windowWidth_array, pulse_pe==0))

            self.windowIntCentralTimeCorrected = ak.drop_none(np.ma.MaskedArray(windowIntCentralTimeCorrected, pulse_pe==0))
            self.windowIntCentralTime = ak.drop_none(np.ma.MaskedArray(windowIntCentralTime, pulse_pe==0))
            self.windowIntLowerBound = ak.drop_none(np.ma.MaskedArray(windowIntLowerBound, pulse_pe==0))
            self.windowIntUpperBound = ak.drop_none(np.ma.MaskedArray(windowIntUpperBound, pulse_pe==0))



            # print("Integrated charge:", self.window_int_charge)
        # print("At the end of the loop, the int charge is", self.window_int_charge)


    def calculate_signal_times(self):
        """Finds the signal time of each waveform as the interpolated time before the peak when the voltage passes 0.4*[peak voltage]"""
        if self.peak_locations is None:
            self.find_peaks()
        thresholds = self.peak_rise_time_fraction*self.peak_voltages
        below_threshold = self.waveforms < thresholds
        indices = np.arange(self.waveforms.shape[1])
        after_peak_rise = np.cumsum(((indices < self.peak_locations) & below_threshold)[:, ::-1], axis=1)[:, ::-1] == 0
        peak_rise_location = np.sum(after_peak_rise == 0, axis=1, keepdims=True)  # count how many samples are before the peak rise
        self.signal_times = peak_rise_location.astype(np.float64)
        amp_range = np.take_along_axis(self.waveforms, np.column_stack([peak_rise_location - 1, peak_rise_location]), axis=1)
        diff = np.diff(amp_range)
        is_nonzero = diff != 0
        self.signal_times[is_nonzero] -= (amp_range[is_nonzero.squeeze(), [1]] - thresholds[is_nonzero]) / diff[is_nonzero]  # interpolate
        self.signal_times = self.signal_times*self.ns_per_sample + self.time_offset
        return self.signal_times

    def integrate_charges(self):
        """Calculates the integrated charge of each waveform's range around the peak that is more than 3x the standard deviation of the pedestal"""
        if self.peak_locations is None:
            self.find_peaks()
        below_threshold = self.smoothed_waveforms < np.maximum(0.05*self.peak_voltages, 3*self.my_pedestal_sigmas)
        indices = np.arange(self.smoothed_waveforms.shape[1])
        before_peak_end = np.cumsum((indices > self.peak_locations) & below_threshold, axis=1) == 0
        after_peak_start = np.cumsum(((indices < self.peak_locations) & below_threshold)[:, ::-1], axis=1)[:, ::-1] == 0
        self.integrated_charges = np.sum(self.waveforms * (before_peak_end & after_peak_start), axis=1, keepdims=True)
        self.integrated_charges *= self.ns_per_sample / self.impedance
        return self.integrated_charges

    def integrate_whole_waveform(self):
        """Calculates the total charge in the analysis window, to be used for two particle event
        identification, as per Dean Karlen's study:  https://wcte.hyperk.ca/wg/beam/meetings/2023/20231211/meeting/beam-structure-and-two-particle-events/beam_structure_v5.pdf
        There is a correction of 1000, this is not a charge, no impedance correction and get back to ADC counts instead of volts"""
        analysis_waveform = self.smoothed_waveforms[:, self.analysis_bins]
        self.whole_waveform_int = np.sum(analysis_waveform, axis=1, keepdims=True)
        self.whole_waveform_int *= 1 / (self.voltage_scale * 1000)
        self.whole_waveform_int_pe = np.sum(analysis_waveform, axis=1, keepdims=True) / self.PMTgain
        return self.whole_waveform_int


    def run_analysis(self):
        """Finds each waveform peak, calculates the time, integrates the charge, and checks if it is over threshold"""
        self.find_peaks()
        self.calculate_signal_times()
        self.integrate_charges()
        self.is_over_threshold = self.peak_voltages > self.threshold
        self.count_peaks()
        self.find_all_peak_voltages()
        self.calculate_all_signal_times()
        self.calculate_all_pulse_charges()
        # self.integrate_charge_in_window()
        self.integrate_whole_waveform()

    def run_window_analysis(self, df_TOF1_hitTimes, df_PMT):
        #df is the ntuple holding the hit timings
        self.df_TOF1_hitTimes = df_TOF1_hitTimes
        self.df_PMT = df_PMT
        self.find_peaks()
        self.integrate_charge_in_window()
