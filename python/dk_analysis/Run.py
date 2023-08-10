import sys
import pickle
from pathlib import Path
import uproot as ur
import gc
import numpy as np

sys.path.insert(0, '../new_analysis')

import channel_map as cm
from waveform_analysis import WaveformAnalysis


class Run:
    """ A class that contains the standard and additional analysis objects, and pickle file open/save methods"""

    def __init__(self, filename, config):
        """Constructor"""

        self.config = config
        self.analyses = []
        self.user = {}
        self.times = {}

        # read data, one channel at a time, do standard analysis on that channel
        run_file = ur.open(filename)

        # save triggerTime and timeStamp information
        if 'triggerTime' in run_file['midas_data_D300']:
            for digi in range(4):
                self.times[digi] = {}
                for field in ['spillNumber','timeStamp','triggerTime']:
                    self.times[digi][field] = run_file['midas_data_D30' + str(digi)][field].array().to_numpy()

        n_channels = self.config["NumberOfChannels"]
        for i_c in range(n_channels):
            digitizer_id = int(i_c / 8)
            d_channel = i_c - 8 * digitizer_id
            waveforms = run_file['midas_data_D30' + str(digitizer_id)]['Channel' + str(d_channel)].array().to_numpy()

            analysis = WaveformAnalysis(waveforms,
                                        threshold=config["Thresholds"][i_c],
                                        analysis_window=(
                                            config["AnalysisWindowLow"][i_c], config["AnalysisWindowHigh"][i_c]),
                                        pedestal_window=(
                                            config["PedestalWindowLow"][i_c], config["PedestalWindowHigh"][i_c]),
                                        reverse_polarity=(config["Polarity"][i_c] == 0),
                                        voltage_scale=config["VoltageScale"],
                                        time_offset=config["TimeOffset"][i_c])

            # disable for hodoscope running - fails
            if n_channels == 19:
                analysis.find_peaks()
                analysis.calculate_signal_times()
                analysis.integrate_charges()
                analysis.is_over_threshold = analysis.peak_voltages > analysis.threshold

            # remove large copies of data - keep raw_waveforms in int32
            del analysis.smoothed_waveforms
            del analysis.waveforms
            # gc.collect()   not necessary

            # keep only compact version of the raw_waveforms
            analysis.raw_waveforms = waveforms.astype(np.int16)

            self.analyses.append(analysis)

    def save_file(self, filename):
        """
        Save a copy of the current run to a file. The run can be restored
        at a late time using Run.open_file(filename).

        If no extension is provided, the default extension, .dk is added.

        Parameters
        ----------
        filename : Path or str
            name of file to save run

        Returns
        -------
        None.

        """

        try:
            filepath = Path(filename).resolve()
        except:
            raise TypeError('Input arg could not be converted to a valid path: {}' +
                            '\n It must be a str or Path-like.'.format(filename))

        if len(filepath.suffix) < 2:
            filepath = filepath.with_suffix('.dk')

        with open(filepath, 'wb') as f:
            pickle.dump(self, f, protocol=4)

    @classmethod
    def open_file(cls, filepath):
        """
        Restore a run that was saved to a file using Run.save_file(filename)

        Parameters
        ----------
        filepath : Path or str
            name of existing run file to open

        Returns
        -------
        Run
            The run object saved in the file

        """

        try:
            filepath = Path(filepath).resolve()
        except:
            raise TypeError('Input arg could not be converted to a valid path: {}' +
                            '\n It must be a str or Path-like.'.format(filepath))

        if not filepath.exists():
            raise ValueError('Filepath does not exist: {}'.format(filepath))

        if len(filepath.suffix) < 2:
            filepath = filepath.with_suffix('.dk')

        with open(filepath, 'rb') as f:
            return pickle.load(f)
