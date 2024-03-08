#!/usr/bin/env python3
import uproot as ur
import os
import numpy as np
import matplotlib.pyplot as plt
import json5
from waveform_analysis import WaveformAnalysis

base_dir = 'plt_waveforms_787'
evts_plt_peak_good = np.array([  1,   4,   9,  11,  13,  27,  39,  42,  46,  52,  53,  62,  64,
        74,  84,  93,  98, 106, 110, 111, 113, 116, 122, 152, 153, 164,
       166, 169, 186, 191, 198, 215, 229, 231, 236, 237, 243, 246, 249,
       252])
evts_plt_peak_bad = np.array([11109])
evts_plt_tail_bad = np.array([  165,   266,   282,   303,   475,   499,   529,   704,   761,
         863,  1058,  1149,  1345,  1351,  1796,  1818,  1853,  2120,
        2230,  2381,  2531,  2580,  2687,  2864,  3078,  3084,  3108,
        3319,  3326,  3484,  3522,  3546,  3687,  3712,  4041,  4081,
        4187,  4398,  4501,  4679,  4779,  5015,  5088,  5111,  5346,
        5386,  5434,  5553,  5788,  6196,  6206,  6323,  6417,  6631,
        6641,  6738,  6799,  7184,  7220,  7319,  7850,  8315,  8589,
        8601,  8670,  8736,  8914,  8998,  9199,  9260,  9305,  9458,
        9577,  9651,  9678,  9760,  9771,  9964, 10011, 10122, 10604,
       10661, 10886, 10895, 11053, 11066, 11078, 11203, 11424, 11456,
       11576, 11655, 11661, 11680, 11863, 12103, 12221, 12263, 12356,
       12387, 12470, 12989, 13014, 13055, 13152, 13335, 13398, 13676,
       13768, 13884, 14164, 14200, 14218, 14381, 14463, 15141, 15422,
       15454, 15658, 15743, 16635, 16684, 16764, 16896, 17036, 17209,
       17303, 17351, 17408, 17774, 17975, 18058, 18112, 18131, 18309,
       18446, 18474, 18709, 18901, 19127, 19160, 19330, 19847, 19886,
       19988, 20021, 20092, 20490, 20571, 20642, 20915, 20982, 21094,
       21169, 21361, 21380, 21480, 21725, 21853, 21920, 21942, 22084,
       22277, 22311, 22338, 22381, 22515, 22518, 22625, 22935, 23066,
       23081, 23086, 23346, 23456, 23596, 23612, 23617, 23861, 23911,
       24071, 24360, 24399, 24561, 24627, 24658, 24758, 24777])
evts_787_tail = np.array([ 25348,  76023,  81508,  87006,  87047, 110295, 115180, 119639,
       153891, 158209])
evts_plt = evts_787_tail

run_file = ur.open("hodoscope/root_run_000787.root")

waveforms_TrigScint = run_file['midas_data_D300']['Channel6'].array().to_numpy()

waveforms_TOF00 = run_file['midas_data_D301']['Channel0'].array().to_numpy()
waveforms_TOF01 = run_file['midas_data_D301']['Channel1'].array().to_numpy()
waveforms_TOF02 = run_file['midas_data_D301']['Channel2'].array().to_numpy()
waveforms_TOF03 = run_file['midas_data_D301']['Channel3'].array().to_numpy()
waveforms_TOF10 = run_file['midas_data_D301']['Channel4'].array().to_numpy()
waveforms_TOF11 = run_file['midas_data_D301']['Channel5'].array().to_numpy()
waveforms_TOF12 = run_file['midas_data_D301']['Channel6'].array().to_numpy()
waveforms_TOF13 = run_file['midas_data_D301']['Channel7'].array().to_numpy()

waveforms_PbGlass = run_file['midas_data_D302']['Channel0'].array().to_numpy()
waveforms_HD14 = run_file['midas_data_D302']['Channel7'].array().to_numpy()

events2 = run_file['midas_data_D302']['eventNumber'].array().to_numpy()
spills2 = run_file['midas_data_D302']['spillNumber'].array().to_numpy()
# print(f"{len(events)} events")
# for evt in events:
#     print(evt)
print(f"{len(waveforms_TOF00)} waveforms with shape {waveforms_TOF00.shape}")

config = json5.load(open("../../config/config_hodoscope.json"))['WaveAnalysis']
analysis = WaveformAnalysis(waveforms_HD14,
                            threshold=config["Thresholds"][0],
                            analysis_window=(config["AnalysisWindowLow"][0], config["AnalysisWindowHigh"][0]),
                            pedestal_window=(config["PedestalWindowLow"][0], config["PedestalWindowHigh"][0]),
                            reverse_polarity=(config["Polarity"][0]==0),
                            voltage_scale=config["VoltageScale"],
                            time_offset=config["TimeOffset"][0])
n_peaks, peak_times, peak_voltages = analysis.find_peaks()
awin_low  = config["AnalysisWindowLow"][0]
awin_high = config["AnalysisWindowHigh"][0]
print(f"Analysis window = ({awin_low}, {awin_high})")
print("Number of peaks: ",n_peaks.shape,"events")
print("Number of peak times: ",peak_times.shape,"events")
print("Number of peak voltages: ",peak_voltages.shape,"events")

wf_names = ['TrigScint',
            'TOF00', 'TOF01', 'TOF02', 'TOF03', 
            'TOF10', 'TOF11', 'TOF12', 'TOF13',
            'PbGlass', 'HD14']
wf_waveforms = [waveforms_TrigScint,
                waveforms_TOF00, waveforms_TOF01, waveforms_TOF02, waveforms_TOF03,
                waveforms_TOF10, waveforms_TOF11, waveforms_TOF12, waveforms_TOF13,
                waveforms_PbGlass, waveforms_HD14]
for wf_name, wf_waveform in zip(wf_names, wf_waveforms):
    print(f"[{wf_name}]: {len(wf_waveform)} waveforms")

for evt in evts_plt[0:10]:

    print(f"Plotting waveforms for event {evt} (event {events2[evt]}, spill {spills2[evt]})...")
    print(f"Found pedestal of {analysis.pedestals[evt]/analysis.voltage_scale}")
    print(f"Found pedestal sigma of {analysis.pedestal_sigmas[evt]/analysis.voltage_scale}")
    print("Found Pb peak times of:",peak_times[evt])
    print("Found Pb peak voltages of:",peak_voltages[evt]/analysis.voltage_scale)

    for wf_name, wf_waveform in zip(wf_names, wf_waveforms):

        fig, axs = plt.subplots(1, 1, figsize=(12,4))
        fig.tight_layout()
        axs.plot(wf_waveform[evt][:])
        plt.title(f"{wf_name}")
        plt.xlabel("Sample",fontsize=14)
        plt.ylabel("ADC counts",fontsize=14)

        out_dir = f"{base_dir}/{evt}"
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        plt.savefig(f"{out_dir}/wf_{wf_name}_{evt}.pdf", bbox_inches='tight')
        plt.close()

# n_peaks = analysis.count_peaks()

# has_1_peak = np.where(n_peaks == 1)
# has_2_peaks = np.where(n_peaks == 2)
# has_more_peaks = np.where(n_peaks > 2)

# rows = 5
# fig, axs = plt.subplots(rows, 3, figsize=(12,2*rows))
# fig.tight_layout()
# axs[0][0].set_title("1 peak")
# axs[0][1].set_title("2 peaks")
# axs[0][2].set_title("3+ peaks")
# for i in range(rows):
#     axs[i][0].plot(waveforms[has_1_peak][i, analysis.analysis_bins])
#     axs[i][1].plot(waveforms[has_2_peaks][i, analysis.analysis_bins])
#     axs[i][2].plot(waveforms[has_more_peaks][i, analysis.analysis_bins])
# plt.show()
