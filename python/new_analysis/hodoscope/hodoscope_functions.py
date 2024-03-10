import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D

from scipy.optimize import curve_fit
from scipy.stats import norm, moyal

# Define the custom order for the keys
custom_order = [
    'ACT0L', 'ACT0R', 'ACT1L', 'ACT1R', 'ACT3L', 'ACT3R',
    'TOF00', 'TOF01', 'TOF02', 'TOF03', 'TOF10', 'TOF11', 'TOF12', 'TOF13',
    'TriggerScint',
    'PbGlass',
    'HD0', 'HD1', 'HD2', 'HD3', 'HD4', 'HD5', 'HD6', 'HD7', 'HD8', 'HD9', 'HD10', 'HD11', 'HD12', 'HD13', 'HD14'
]

# Define the computed electron energy for hitting each hodoscope element.
elec_hit_momenta_MC_fieldmap = {0: 0.09783230585487263,
                    1: 0.11089383662305125,
                    2: 0.11718751862628744,
                    3: 0.1357407927028879,
                    4: 0.13997438543024,
                    5: 0.1471055825691591,
                    6: 0.15399394364580632,
                    7: 0.17768442931536754,
                    8: 0.1954001338291036,
                    9: 0.21160634230382658,
                    10: 0.2352525614724067,
                    11: 0.2544546042946958,
                    12: 0.2914536355311684,
                    13: 0.32090114397989367,
                    14: 0.3746842330143378}

elec_hit_momenta_MC_basicfield = {0: 0.15582261018870644,
                                1: 0.16651505160688465,
                                2: 0.17289996548544717,
                                3: 0.17841390695053236,
                                4: 0.18295397261505822,
                                5: 0.1978669429349461,
                                6: 0.20530606915601174,
                                7: 0.22249044868699686,
                                8: 0.23256616677748593,
                                9: 0.24792568763571388,
                                10: 0.26542390993691517,
                                11: 0.2862366273068532,
                                12: 0.3116743632237849,
                                13: 0.33559257093775985,
                                14: 0.38199036544635756}

elec_hit_momenta_calc = {0: 0.15058755351819803,
                    1: 0.1632022173955016,
                    2: 0.16454507453432862,
                    3: 0.17837784731694356,
                    4: 0.18160773326593496,
                    5: 0.1968708412789219,
                    6: 0.20286783238767242,
                    7: 0.21983444684589068,
                    8: 0.22999315900645517,
                    9: 0.24902380087408993,
                    10: 0.26566623596454925,
                    11: 0.28724821930912336,
                    12: 0.31449359776978913,
                    13: 0.3393035715640493,
                    14: 0.38511992183555954}

elec_hit_momenta = elec_hit_momenta_calc

# Define the computed momentum range for each hodoscope element.
momentum_range = {0: 0.0063047132725456145,
                    1: 0.006932943602899561,
                    2: 0.007724812094032252,
                    3: 0.00845980150260986,
                    4: 0.009627122746240724,
                    5: 0.010498277293921526,
                    6: 0.01225191220751115,
                    7: 0.013300445377746856,
                    8: 0.016009392120709398,
                    9: 0.017294611862594872,
                    10: 0.02164736408549034,
                    11: 0.023256557966043867,
                    12: 0.030648203744867764,
                    13: 0.03271156068397968,
                    14: 0.04630141916160452}

def ntuple_to_pd(filename):
    """
    Converts an ntuple to a single Pandas dataframe containing information about the largest peaks
    in each detector.
    """
    
    # Open the TTree anaTree and get all keys
    events = uproot.open("{}:anaTree".format(filename))
    main_keys = events.keys()
    
    # Construct the analysis dataframe
    df = pd.DataFrame()
    for key in main_keys:
        
        # Skip the nChannels key
        if(key == 'nChannels'):
            continue

        arr = events[key].array(library='np').squeeze()

        # Convert the array to a DataFrame
        df_key = pd.DataFrame(arr)

        # Rename the columns to include the key name
        df_key.columns = [f"{key}{i}" for i in df_key.columns]

        # Concatenate the new DataFrame to the existing one
        df = pd.concat([df, df_key], axis=1)
        
    return df

def ntuple_to_pd_multipeak(filename):
    """
    Converts an ntuple to multiple Pandas dataframes containing per-peak information.
    """
    
    # Open the TTree anaTree and get all keys
    events = uproot.open("{}".format(filename))
    main_keys = events.keys()
    
    # Construct the analysis dataframes for each detector element.
    df_dict = {}
    for key in main_keys:
        
        key_noversion = key.split(';')[0]
        
        # Skip the EventInfo key
        if(key_noversion == 'EventInfo'):
            continue
        
        print("Processing dataframe for",key,"...")
        
        # Get the number of peaks and timestamps
        nPeaks        = events[key]['nPeaks'].array()
        timeStamp     = events[key]['timeStamp'].array()
        triggerTime   = events[key]['triggerTime'].array()
        Pedestal      = events[key]['Pedestal'].array()
        PedestalSigma = events[key]['PedestalSigma'].array()
        PeakVoltage   = events[key]['PeakVoltage'].array()
        PeakTime      = events[key]['PeakTime'].array()
        SignalTime    = events[key]['SignalTime'].array()
        IntCharge     = events[key]['IntCharge'].array()
        
        # Iterate through the array elements and save information for each peak.
        l_evt, l_ipk, l_nPeaks, l_timeStamp, l_triggerTime = [], [], [], [], []
        l_Pedestal, l_PedestalSigma = [], []
        l_PeakVoltage, l_PeakTime, l_SignalTime, l_IntCharge = [], [], [], []
        for evt, (npk, tstamp, ttime, pedestal, spedestal, pkv, pkt, sigt, chg) in enumerate(zip(nPeaks,timeStamp,triggerTime,Pedestal,PedestalSigma,PeakVoltage,PeakTime,SignalTime,IntCharge)):
            
            # Update the lists for each peak.
            if(npk == 0):
                l_evt.append(evt)
                l_ipk.append(-1)
                l_nPeaks.append(npk)
                l_timeStamp.append(tstamp)
                l_triggerTime.append(ttime)
                l_Pedestal.append(pedestal)
                l_PedestalSigma.append(spedestal)
                l_PeakVoltage.append(-1)
                l_PeakTime.append(-1)
                l_SignalTime.append(-1)
                l_IntCharge.append(-1)
            else:
                for ipk in range(npk):
                    l_evt.append(evt)
                    l_ipk.append(ipk)
                    l_nPeaks.append(npk)
                    l_timeStamp.append(tstamp)
                    l_triggerTime.append(ttime)
                    l_Pedestal.append(pedestal)
                    l_PedestalSigma.append(spedestal)
                    l_PeakVoltage.append(pkv[ipk])
                    l_PeakTime.append(pkt[ipk])
                    l_SignalTime.append(sigt[ipk])
                    l_IntCharge.append(chg[ipk])
        
        # Create a new dataframe.
        df = pd.DataFrame({'event':  l_evt,
                           'iPeak': l_ipk,
                           'nPeaks': l_nPeaks,
                           'timeStamp': l_timeStamp,
                           'triggerTime': l_triggerTime,
                           'Pedestal': l_Pedestal,
                           'PedestalSigma': l_PedestalSigma,
                           'PeakVoltage': l_PeakVoltage,
                           'PeakTime': l_PeakTime,
                           'SignalTime': l_SignalTime,
                           'IntCharge': l_IntCharge
                          })

        # Set this as the dataframe or concatenate it to the one that is already there.
        if(key_noversion in df_dict):
            last_evt = df_dict[key_noversion]['event'].values[-1]+1
            df['event'] = df['event'] + last_evt
            print("Concatenating to",key_noversion,"starting with event number",last_evt)
            df_dict[key_noversion] = pd.concat([df_dict[key_noversion], df], ignore_index=True)
        else:
            df_dict[key_noversion] = df
            
    # Add EventInfo information.
    df = pd.DataFrame( {'RunNumber':   events['EventInfo']['RunNumber'].array(library='np'),
                        'EventNumber': events['EventInfo']['EventNumber'].array(library='np'),
                        'SpillNumber': events['EventInfo']['SpillNumber'].array(library='np')})
    df_dict['EventInfo'] = df
        
    return df_dict

def plot_2D_histogram(df_dict, detector1, quantity1, detector2, quantity2, evt_list=None, base_dir=".", rnum = 0, select_nonzero_peaks=False, per_evt_quantity=False, logscale=False, nbins=60, range=[[0,100],[0,100]]):
    """
    Generate a 2D histogram for the specified signals.
    """
    
    out_dir = f"{base_dir}/{rnum}"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    print("Plotting 2D histogram for (",detector1,",",quantity1,") vs. (",detector2,",",quantity2,") for all signals in run",rnum)

    # Create a grid of 8 rows x 4 columns
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    
    #fig.title(f'Run {rnum}', fontsize=24, y=1.02)

    # Get the detector element from the dataframe dictionary.
    df1 = df_dict[detector1]
    df2 = df_dict[detector2]

    # Select the events based on the given list.
    if(not(evt_list is None)):
        print("DF1: before selection:",len(df1))
        print("DF2: before selection:",len(df2))

        df1 = df1[df1['event'].isin(evt_list)]
        df2 = df2[df2['event'].isin(evt_list)]

        print("DF1: selecting out",len(evt_list),"events to get",len(df1))
        print("DF2: selecting out",len(evt_list),"events to get",len(df2))


    # Place a cut if we're looking at an HD element and plotting PeakVoltage or IntCharge.
    if(detector1[0:2] == 'HD' and (quantity1 == 'IntCharge' or quantity1 == 'PeakVoltage')):
        df1 = df1[df1[quantity1] > 0.02]
    if(detector2[0:2] == 'HD' and (quantity2 == 'IntCharge' or quantity2 == 'PeakVoltage')):
        df2 = df2[df2[quantity2] > 0.02]

    # SINGLE PEAK ANALYSIS
    # Consider only 1-peak events.
    # df1_single_peak = df1[df1['nPeaks'] == 1]
    # df2_single_peak = df2[df2['nPeaks'] == 1]

    # Find common events
    # common_events = set(df1_single_peak['event']).intersection(set(df2_single_peak['event']))
    # df1_final = df1_single_peak[df1_single_peak['event'].isin(common_events)]
    # df2_final = df2_single_peak[df2_single_peak['event'].isin(common_events)]

    # MULTI-PEAK ANALYSIS
    df1['event_nPeaks'] = df1['event'].astype(str) + "_" + df1['nPeaks'].astype(str)
    df2['event_nPeaks'] = df2['event'].astype(str) + "_" + df2['nPeaks'].astype(str)

    common_events_nPeaks = set(df1['event_nPeaks']).intersection(set(df2['event_nPeaks']))

    df1_final = df1[df1['event_nPeaks'].isin(common_events_nPeaks)]
    df2_final = df2[df2['event_nPeaks'].isin(common_events_nPeaks)]

    print("Shape of df1[quantity1]:", df1_final[quantity1].shape)
    print("Shape of df2[quantity2]:", df2_final[quantity2].shape)
    print("Type of df1[quantity1]:", type(df1_final[quantity1]))
    print("Type of df2[quantity2]:", type(df2_final[quantity2]))
    print("NaNs in df1[quantity1]:", df1_final[quantity1].isna().sum())
    print("NaNs in df2[quantity2]:", df2_final[quantity2].isna().sum())
    print("Infs in df1[quantity1]:", np.isinf(df1_final[quantity1]).sum())
    print("Infs in df2[quantity2]:", np.isinf(df2_final[quantity2]).sum())

    h2d = ax.hist2d(df1_final[quantity1],df2_final[quantity2],bins=[nbins,nbins],norm=matplotlib.colors.LogNorm(), range=range)
    ax.set_xlabel(f"{detector1}/{quantity1}")
    ax.set_ylabel(f"{detector2}/{quantity2}")
    ax.legend()  # Add legend
    
    if(logscale):
        ax.set_yscale('log')

    # Adjust the layout so the plots do not overlap
    plt.tight_layout()
    #plt.show()
    plt.savefig(f"{out_dir}/{detector1}_{quantity1}vs{detector2}_{quantity2}.pdf", bbox_inches='tight')
    plt.close()


def plot_histograms_for_each_signal(df_dict, evt_list=None, base_dir=".", rnum = 0, quantity='nPeaks', select_nonzero_peaks=False, per_evt_quantity=False, logscale=False, nbins=60):
    """
    Generate plots of histograms for each signal.
    """
    
    out_dir = f"{base_dir}/{rnum}"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    print("Plotting all histograms for",quantity,"for all signals in run",rnum)

    # Create a grid of 8 rows x 4 columns
    fig, axes = plt.subplots(8, 4, figsize=(16, 24))
    flat_axes = axes.ravel()
    
    fig.suptitle(f'Run {rnum}', fontsize=24, y=1.02)

    # Iterate based on the custom order
    for key, ax in zip(custom_order, flat_axes):
        df = df_dict[key]

        # Cut on the selected events if an event list is provided.
        if(not(evt_list is None)):
            print(f"[{key}] before selection {len(df)}")
            df_select = df[df['event'].isin(evt_list)]
            print(f"[{key}] selecting out {len(evt_list)} events to get {len(df_select)}")
        
        # Select only non-zero peaks if specified.
        if(select_nonzero_peaks):
            df_select = df_select[df_select['nPeaks'] > 0]
        # Otherwise, this quantity is event-wide: only keep 1 entry for each event.
        else:
            df_select = df_select.drop_duplicates(subset=['event'], keep='first')
            
        # Place a cut if we're looking at an HD element and plotting PeakVoltage or IntCharge.
        if(key[0:2] == 'HD' and (quantity == 'IntCharge' or quantity == 'PeakVoltage')):
            df_select = df_select[df_select[quantity] > -2.0]

        # Plot histogram for the current signal on its corresponding axis
        n, bins, patches = ax.hist(df_select[quantity], bins=nbins, edgecolor='black', alpha=0.7, label=key)  # capture output to use in legend
        #ax.set_title(key)
        ax.set_xlabel(quantity)
        ax.set_ylabel('Counts/bin')
        ax.legend()  # Add legend
        
        if(logscale):
            ax.set_yscale('log')

    # Adjust the layout so the plots do not overlap
    plt.tight_layout()
    #plt.show()
    plt.savefig(f"{out_dir}/{quantity}.pdf", bbox_inches='tight')
    plt.close()

def read_dataframes_from_csv(directory):
    """
    Reads a set of multi-peak dataframes from a directory containing .csv files.
    The dataframes are stored in a dictionary with keys equal to the file names.
    """

    file_names = os.listdir(directory)
    dfs = {}
    for file in file_names:
        if file.endswith('.csv'):
            key = file.split('.csv')[0]
            dfs[key] = pd.read_csv(os.path.join(directory, file))
    return dfs

def compute_statistics(df):
    
    # Compute average statistics
    avg_nPeaks = df['nPeaks'].mean()
    avg_Pedestal = df['Pedestal'].mean()
    avg_PedestalSigma = df['PedestalSigma'].mean()
    
    # For histograms, bin the data and find the bin with the most counts
    def get_hist_peak(data, bins=100):
        hist, bin_edges = np.histogram(data, bins=bins)
        peak_bin = np.argmax(hist)
        peak_value = (bin_edges[peak_bin] + bin_edges[peak_bin + 1]) / 2
        return peak_value

    peak_PeakVoltage = get_hist_peak(df[df['PeakVoltage'] > 0.2]['PeakVoltage'])
    peak_PeakTime = get_hist_peak(df[df['PeakTime'] > 10]['PeakTime'])
    peak_SignalTime = get_hist_peak(df[df['SignalTime'] > 0.2]['SignalTime'])
    peak_IntCharge = get_hist_peak(df[df['IntCharge'] > 0.025]['IntCharge'])
    
    stats = {
        'avg_nPeaks': avg_nPeaks,
        'avg_Pedestal': avg_Pedestal,
        'avg_PedestalSigma': avg_PedestalSigma,
        'peak_PeakVoltage': peak_PeakVoltage,
        'peak_PeakTime': peak_PeakTime,
        'peak_SignalTime': peak_SignalTime,
        'peak_IntCharge': peak_IntCharge
    }
    
    return stats

def plot_statistics_vs_run(statistics_data, statistic_keys, signal_set, base_dir, signal_set_name):
    
    out_dir = f"{base_dir}/summary"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    n_statistics = len(statistic_keys)
    
    # Define number of rows and columns for the subplots.
    n_cols = 1
    n_rows = math.ceil(n_statistics / n_cols)
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4*n_rows))
    
    if n_statistics == 1:
        axes = [axes]
        
    # Generate a palette with as many distinct colors as there are signals.
    #palette = sns.color_palette("husl", len(signal_set))
    
    for idx, statistic_key in enumerate(statistic_keys):
        ax = axes[idx]
        #ax = axes[idx//n_cols, idx%n_cols] if n_rows > 1 else axes[idx]
        
        # Set the color cycle for this axis
        #ax.set_prop_cycle(color=palette)
        
        for idkey,key in enumerate(signal_set):
            marker = 's'
            if(idkey > 9): marker = 'o'
            runs = list(statistics_data.keys())
            values = [statistics_data[run][key][statistic_key] for run in runs]
            ax.plot(runs, values, marker=marker, label=key)
        
        #ax.set_title(f"{statistic_key} vs. Run Number")
        ax.set_xlabel("Run Number",fontsize=14)
        ax.set_ylabel(statistic_key,fontsize=14)
        ax.legend()
        ax.grid(True)
    
    # Handle case when the number of plots is odd
    #if n_statistics % 2 != 0 and n_rows > 1:
        #axes[n_rows-1, 1].axis('off')
    #    axes[n_rows-1].axis('off')
    
    plt.tight_layout()
    #plt.show()
    plt.savefig(f"{out_dir}/{signal_set_name}.pdf", bbox_inches='tight')

# -----------------------------------------------------------------------------------------
# Functions for gamma peak analysis
# -----------------------------------------------------------------------------------------
def get_num_duplicates(arr):
    u, c = np.unique(arr, return_counts=True)
    return np.sum(c[c > 1])

def filter_range(name,df,key,rng,drop=False,debug=False):
    
    # Perform the filter.
    df_filtered = df[df[key].between(*rng)]
    nduplicates = get_num_duplicates(df_filtered['event'].values)
    if(debug): print(f"{name} filter passed",len(df_filtered),"of",len(df),"with",nduplicates,"duplicates")
    
    # Handle duplicates.
    if(nduplicates > 0):
        
        # Drop all events that appear more than once if selected.
        if(drop):
            df_filtered = df_filtered[~df_filtered['event'].duplicated(keep=False)]
        # Otherwise choose the one with the maximum integrated charge.
        else:
            idx = df_filtered.groupby('event')['IntCharge'].idxmax()
            df_filtered = df_filtered.loc[idx]
        
    if(debug): print("--> returning",len(df_filtered),"events")
    return df_filtered

def timing_analysis(df_dict, pb_timing_range, tof0_timing_range, tof0_charge_range, 
                    tof1_timing_range, tof1_charge_range, t2_timing_range, t2_charge_range,
                    act0_timing_range, act0_charge_range, act1_timing_range, act1_charge_range,
                    act3_timing_range, hd_timing_ranges, hd_charge_ranges, low_radiation = False,
                    debug = False):
    
    # Extract the total number of events
    ntot_evts = -1
    ntot_spills = -1
    for det in df_dict:

        # Confirm that each detector dataframe has recorded all events.
        if(det != 'EventInfo'):
            df_det = df_dict[det]
            evts = df_det['event'].values
            if(debug): print(f"[Detector {det}] events go from {evts[0]} to {evts[-1]} for total of {len(np.unique(evts))}")
        # Get the actual total from the 'EventInfo' dataframe.
        else:
            df_evts = df_dict[det]
            evts = df_evts['EventNumber'].values
            spills = df_evts['SpillNumber'].values
            ntot_evts = len(evts)
            ntot_spills = len(np.unique(spills))
    if(debug): print(f"* Found a total of {ntot_evts} events and {ntot_spills} spills.")
    
    # Lead glass
    # --------------------------------------------------------------------------------------------
    # Note: we require 1 lead glass peak in the final analysis
    pb_1peak = filter_range("PbGlass",df_dict['PbGlass'],'nPeaks',(1,1),debug=debug)
    pb_filtered = filter_range("PbGlass",pb_1peak,'PeakTime',pb_timing_range,debug=debug)
    npeaks_pre  = len(df_dict['PbGlass'])
    npeaks_post = len(pb_filtered)
    #if(debug): print(f"-- Peaks {npeaks_post}/{npeaks_pre} = {npeaks_post/npeaks_pre}")

    # TOF elements
    # --------------------------------------------------------------------------------------------
    # TOF0
    tof00_filtered = filter_range("TOF00",df_dict['TOF00'],'PeakTime',tof0_timing_range,debug=debug)
    tof01_filtered = filter_range("TOF01",df_dict['TOF01'],'PeakTime',tof0_timing_range,debug=debug)
    tof02_filtered = filter_range("TOF02",df_dict['TOF02'],'PeakTime',tof0_timing_range,debug=debug)
    tof03_filtered = filter_range("TOF03",df_dict['TOF03'],'PeakTime',tof0_timing_range,debug=debug)

    combined_00_01 = tof00_filtered.merge(tof01_filtered, on='event', suffixes=('_00', '_01'))
    combined_00_01_02 = combined_00_01.merge(tof02_filtered, on='event')
    combined_00_01_02 = combined_00_01_02.rename(columns={'IntCharge': 'IntCharge_02'})
    tof0_combined = combined_00_01_02.merge(tof03_filtered, on='event')
    tof0_combined = tof0_combined.rename(columns={'IntCharge': 'IntCharge_03'})
    tof0_combined['combined_charge'] = (tof0_combined['IntCharge_00'] 
                                    + tof0_combined['IntCharge_01']
                                    + tof0_combined['IntCharge_02']
                                    + tof0_combined['IntCharge_03'])

    tof0_valid = filter_range("TOF0_combined",tof0_combined,'combined_charge',tof0_charge_range,debug=debug)
    tof0_valid.loc[:,'hit_TOF0'] = 1
    if(debug): print()

    # TOF1
    tof10_filtered = filter_range("TOF10",df_dict['TOF10'],'PeakTime',tof1_timing_range,debug=debug)
    tof11_filtered = filter_range("TOF11",df_dict['TOF11'],'PeakTime',tof1_timing_range,debug=debug)
    tof12_filtered = filter_range("TOF12",df_dict['TOF12'],'PeakTime',tof1_timing_range,debug=debug)
    tof13_filtered = filter_range("TOF13",df_dict['TOF13'],'PeakTime',tof1_timing_range,debug=debug)

    combined_10_11 = tof10_filtered.merge(tof11_filtered, on='event', suffixes=('_10', '_11'))
    combined_10_11_12 = combined_10_11.merge(tof12_filtered, on='event')
    combined_10_11_12 = combined_10_11_12.rename(columns={'IntCharge': 'IntCharge_12'})
    tof1_combined = combined_10_11_12.merge(tof13_filtered, on='event')
    tof1_combined = tof1_combined.rename(columns={'IntCharge': 'IntCharge_13'})
    tof1_combined['combined_charge'] = (tof1_combined['IntCharge_10'] 
                                    + tof1_combined['IntCharge_11']
                                    + tof1_combined['IntCharge_12']
                                    + tof1_combined['IntCharge_13'])

    tof1_valid = filter_range("TOF1_combined",tof1_combined,'combined_charge',tof1_charge_range,debug=debug)
    if(not low_radiation): tof1_valid.loc[:,'hit_TOF1'] = 1
    if(debug): print()

    # Filter T2
    # --------------------------------------------------------------------------------------------
    t2_filtered = filter_range("T2",df_dict['TriggerScint'],'PeakTime',t2_timing_range,debug=debug)
    t2_valid = filter_range("T2",t2_filtered,'IntCharge',t2_charge_range,debug=debug)
    t2_valid.loc[:,'hit_T2'] = 1
    if(debug): print()
    # --------------------------------------------------------------------------------------------

    # ACT elements
    # --------------------------------------------------------------------------------------------
    # ACT0
    act0l_filtered = filter_range("ACT0L",df_dict['ACT0L'],'PeakTime',act0_timing_range,debug=debug)
    act0r_filtered = filter_range("ACT0R",df_dict['ACT0R'],'PeakTime',act0_timing_range,debug=debug)

    act0_combined = act0l_filtered.merge(act0r_filtered, on='event', suffixes=('_L', '_R'))
    act0_combined['combined_charge'] = act0_combined['IntCharge_L'] + act0_combined['IntCharge_R']

    act0_valid = filter_range("ACT0_combined",act0_combined,'combined_charge',act0_charge_range,debug=debug)
    if(not low_radiation): act0_valid.loc[:,'hit_ACT0'] = 1

    # ACT1
    act1l_filtered = filter_range("ACT1L",df_dict['ACT1L'],'PeakTime',act1_timing_range,debug=debug)
    act1r_filtered = filter_range("ACT1R",df_dict['ACT1R'],'PeakTime',act1_timing_range,debug=debug)

    act1_combined = act1l_filtered.merge(act1r_filtered, on='event', suffixes=('_L', '_R'))
    act1_combined['combined_charge'] = act1_combined['IntCharge_L'] + act1_combined['IntCharge_R']

    act1_valid = filter_range("ACT1_combined",act1_combined,'combined_charge',act1_charge_range,debug=debug)
    act1_valid.loc[:,'hit_ACT1'] = 1

    # ACT3
    act3l_filtered = df_dict['ACT3L'][df_dict['ACT3L']['nPeaks'] == 0]
    act3r_filtered = df_dict['ACT3R'][df_dict['ACT3R']['nPeaks'] == 0]

    act3_combined = act3l_filtered.merge(act3r_filtered, on='event', suffixes=('_L', '_R'))
    act3_valid = act3_combined.copy()
    act3_valid.loc[:,'nohit_ACT3'] = 1
    if(debug): print("ACT3 total number of valid events:",len(act3_valid))
    # --------------------------------------------------------------------------------------------

    # Hodoscope elements
    # --------------------------------------------------------------------------------------------
    hd_dfs = {}

    # Create the filtered dataframes (containing peaks over the threshold at the correct time).
    for i in range(15):  # 0 to 14 inclusive
        hd_key = f'HD{i}'
        hit_col_name = f'hit_{hd_key}'
        hd_time_range = hd_timing_ranges[hd_key]
        hd_charge_range = hd_charge_ranges[hd_key]
        
        # Filtering
        hd_filtered = filter_range(hd_key,df_dict[hd_key],'PeakTime',hd_time_range,debug=debug)
        hd_filtered = filter_range(hd_key,hd_filtered,'IntCharge',hd_charge_range,debug=debug)
        # hd_filtered = df_dict[hd_key][(df_dict[hd_key]['PeakTime'].between(*hd_time_range)) & 
        #                         (df_dict[hd_key]['IntCharge'].between(*hd_charge_range))]
        if(debug): print(f"HD{i}: {len(hd_filtered)} of {len(df_dict[hd_key])} events after filter")
        
        # Assign a binary value indicating a hit
        if not hd_filtered.empty:
            hd_filtered.loc[:,hit_col_name] = 1
        else:
            hd_filtered.loc[:,hit_col_name] = 0

        hd_dfs[hd_key] = hd_filtered

    # Merge all HD dataframes
    combined_hd_df = hd_dfs['HD0'][['event', 'hit_HD0']]
    for i in range(1, 15):
        hd_key = f'HD{i}'
        hit_col_name = f'hit_{hd_key}'
        
        combined_hd_df = combined_hd_df.merge(
            hd_dfs[hd_key][['event', hit_col_name]],
            on='event',
            how='outer'
        ).fillna(0)

    # Calculate total hits per event
    combined_hd_df['total_hits'] = combined_hd_df.filter(like='hit_').sum(axis=1)

    # Filter events with a hit count of 1
    hd_valid_events = combined_hd_df[combined_hd_df['total_hits'] == 1]
    if(debug): print(f"Number of hd_valid_events = {len(hd_valid_events)}")

    u, c = np.unique(combined_hd_df['event'].values, return_counts=True)
    if(debug): print("Duplicates in combined HD dataframe:",np.sum(c[c > 1]))
    # --------------------------------------------------------------------------------------------

    # Merge all relevant dataframes
    pb_filtered.columns = ['LG_' + col if col != 'event' else col for col in pb_filtered.columns]
    final_df = pb_filtered.merge(hd_valid_events, on='event', how='inner')
    ntagged_evts = len(final_df[final_df.total_hits == 1])
    if(debug): print(f"Final df number of events = {len(final_df)}; {ntagged_evts} tagged events")
    if(not low_radiation): final_df = final_df.merge(act0_valid[['event', 'hit_ACT0']], on='event', how='left')
    final_df = final_df.merge(act1_valid[['event', 'hit_ACT1']], on='event', how='left')
    final_df = final_df.merge(act3_valid[['event', 'nohit_ACT3']], on='event', how='left')
    final_df = final_df.merge(tof0_valid[['event', 'hit_TOF0']], on='event', how='left')
    if(not low_radiation): final_df = final_df.merge(tof1_valid[['event', 'hit_TOF1']], on='event', how='left')
    final_df = final_df.merge(t2_valid[['event', 'hit_T2']], on='event', how='left')
    if(debug): print(f"{len(final_df[final_df.total_hits == 1])} tagged events")

    # Fill NaN values in the hit columns with 0
    if(low_radiation): final_df[['hit_ACT1', 'nohit_ACT3', 'hit_TOF0', 'hit_T2']] = final_df[['hit_ACT1', 'nohit_ACT3', 'hit_TOF0', 'hit_T2']].fillna(0)
    else: final_df[['hit_ACT0', 'hit_ACT1', 'nohit_ACT3', 'hit_TOF0', 'hit_TOF1', 'hit_T2']] = final_df[['hit_ACT0', 'hit_ACT1', 'nohit_ACT3', 'hit_TOF0', 'hit_TOF1', 'hit_T2']].fillna(0)

    return final_df, ntot_evts, ntot_spills, ntagged_evts

def line(x, m, b):
    y = m*x + b
    return y

def gaussian(x, amplitude, mean, stddev):
        return amplitude * norm.pdf(x, loc=mean, scale=stddev)

def resfxn(x, F, C):
        return F * x**0.5 + C

def gamma_peak_plots(final_df, run, run_momentum, nbins, range, base_dir=".", timing_cuts = False, low_radiation = False):
    """
    Analyze the gamma peaks.
    """

    # Set up the output directory
    out_dir = f"{base_dir}/{run}"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    ax_lbl_font_size = 20
    ax_tick_font_size = 14
    ax_legend_font_size = 20

    # Apply timing cuts if specified.
    if(timing_cuts):

        if(not low_radiation):
            cuts = (final_df.total_hits == 1) & (final_df.hit_ACT0 == 1) & \
                        (final_df.hit_ACT1 == 1) & (final_df.nohit_ACT3 == 1) & \
                        (final_df.hit_TOF0 == 1) & (final_df.hit_TOF1 == 1) & (final_df.hit_T2 == 1)
        else:
            cuts = (final_df.total_hits == 1) & \
                (final_df.hit_ACT1 == 1) & (final_df.nohit_ACT3 == 1) & \
                (final_df.hit_TOF0 == 1) & (final_df.hit_T2 == 1)
    else:
        cuts = (final_df.total_hits == 1)

    # Get the charge arrays
    chg_hd14 = final_df[(final_df.hit_HD14 == 1) & cuts]['LG_IntCharge'].values
    chg_hd13 = final_df[(final_df.hit_HD13 == 1) & cuts]['LG_IntCharge'].values
    chg_hd12 = final_df[(final_df.hit_HD12 == 1) & cuts]['LG_IntCharge'].values
    chg_hd11 = final_df[(final_df.hit_HD11 == 1) & cuts]['LG_IntCharge'].values
    chg_hd10 = final_df[(final_df.hit_HD10 == 1) & cuts]['LG_IntCharge'].values
    chg_hd9 = final_df[(final_df.hit_HD9 == 1) & cuts]['LG_IntCharge'].values
    chg_hd8 = final_df[(final_df.hit_HD8 == 1) & cuts]['LG_IntCharge'].values
    chg_hd7 = final_df[(final_df.hit_HD7 == 1) & cuts]['LG_IntCharge'].values
    chg_hd6 = final_df[(final_df.hit_HD6 == 1) & cuts]['LG_IntCharge'].values
    chg_hd5 = final_df[(final_df.hit_HD5 == 1) & cuts]['LG_IntCharge'].values
    chg_hd4 = final_df[(final_df.hit_HD4 == 1) & cuts]['LG_IntCharge'].values
    chg_hd3 = final_df[(final_df.hit_HD3 == 1) & cuts]['LG_IntCharge'].values
    chg_hd2 = final_df[(final_df.hit_HD2 == 1) & cuts]['LG_IntCharge'].values
    chg_hd1 = final_df[(final_df.hit_HD1 == 1) & cuts]['LG_IntCharge'].values
    chg_hd0 = final_df[(final_df.hit_HD0 == 1) & cuts]['LG_IntCharge'].values

    chg_arrays = [chg_hd0, chg_hd1, chg_hd2, chg_hd3, chg_hd4, chg_hd5, chg_hd6, chg_hd7, chg_hd8,
                  chg_hd9, chg_hd10, chg_hd11, chg_hd12, chg_hd13, chg_hd14]
    labels = ["H0", "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12", "H13", "H14"]

    # Create a list of colors for the histogram of each hodoscope element.
    # hues = np.linspace(0, 1, 16)[:-1]  # Exclude the last hue since it's equal to the first
    # order = np.array([0, 9, 5, 14, 3, 12, 8, 1, 10, 6, 13, 2, 11, 7, 4])
    # ordered_hues = hues[order]
    colors = sns.hls_palette(15, l=0.4)

    # Set up the figures.

    # Figure showing only peak fits.
    fig_peaks, ax_peaks = plt.subplots(1,1,figsize=(12,7))

    # Figure showing all peak fits in a single figure.
    fig_allpeaks, axes_allpeaks = plt.subplots(3, 5, figsize=(40,18))
    flat_axes_allpeaks = axes_allpeaks.ravel()

    # Figure showing all key plots.
    fig_keyplots, axes_keyplots = plt.subplots(2, 2, figsize=(28, 14))
    flat_axes_keyplots = axes_keyplots.ravel()
    fig_keyplots.suptitle("Run {}, p = +{} MeV/c".format(run,run_momentum), fontsize=32, y=0.95)

    # Fit all gamma peaks.
    initial_params = [1000, np.mean(chg_hd14), np.std(chg_hd14)]
    fit_means, fit_smeans, fit_sigmas, fit_ssigmas = [], [], [], []
    gamma_peaks_max = -1
    for ipeak,(arr,label,color) in enumerate(zip(chg_arrays[::-1],labels[::-1],colors)):

        # Make the histogram for this hodoscope element
        hist, bin_edges = np.histogram(arr, bins=nbins, range=range)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Gaussian fit
        #print("Fitting gamma peak",label,"with initial params",*initial_params)
        popt, pcov = curve_fit(gaussian, bin_centers, hist, p0=initial_params)
        perr = np.sqrt(np.diag(pcov))
        fit_curve = gaussian(bin_centers, *popt)
        #print("-- fit result",*popt,"with errors",perr)
        initial_params = popt
        fit_means.append(popt[1])
        fit_smeans.append(perr[1])
        fit_sigmas.append(popt[2])
        fit_ssigmas.append(perr[2])

        # Plot the histogram (individual peak plots).
        lbl1 = "ELEMENT {}".format(label)
        lbl2 = "$\mu$ = {:.4f} $\pm$ {:.4f}".format(popt[1],perr[1])
        lbl3 = "$\sigma$ = {:.4f} $\pm$ {:.4f}".format(popt[2],perr[2])
        legend_elements = [Line2D([0], [0], color='none', lw=0, label=lbl1),
                           Line2D([0], [0], color='none', lw=0, label=lbl2),
                           Line2D([0], [0], color='none', lw=0, label=lbl3)]

        flat_axes_allpeaks[ipeak].bar(bin_edges[:-1], hist, width=np.diff(bin_edges)[0], align='edge', color='white')
        flat_axes_allpeaks[ipeak].plot(bin_edges[:-1], hist, color='black', drawstyle='steps-post')
        flat_axes_allpeaks[ipeak].plot(bin_centers, fit_curve, '--', color='red', linewidth=1.5, alpha=1.0)
        flat_axes_allpeaks[ipeak].set_xlabel('Lead glass charge [arbitrary units]',fontsize=ax_lbl_font_size)
        flat_axes_allpeaks[ipeak].set_ylabel('Counts/bin',fontsize=ax_lbl_font_size)
        flat_axes_allpeaks[ipeak].tick_params(axis="x", labelsize=ax_tick_font_size)
        flat_axes_allpeaks[ipeak].tick_params(axis="y", labelsize=ax_tick_font_size)
        leg = flat_axes_allpeaks[ipeak].legend(handles=legend_elements, frameon=True, handlelength=0, fontsize=16)
        for i, text in enumerate(leg.get_texts()):
            if i == 0:
                text.set_weight('bold')
            text.set_horizontalalignment('right')

        # Plot the histogram (peak plot).
        lbl = '{}, $\mu$ = {:.3f}'.format(label, popt[1])
        # h_peaks, bedge_peaks = np.histogram(arr, bins=nbins)
        # ax_peaks.bar(bin_edges[:-1], hist, label=lbl, width=np.diff(bin_edges)[0], align='edge', color=color)
        # ax_peaks.plot(bin_edges[:-1], hist, color='black', drawstyle='steps-post')
        hpeaks = ax_peaks.hist(arr, bins=nbins, alpha=1.0, label=lbl, range=range)
        ax_peaks.plot(bin_centers, fit_curve, '-', color='red', linewidth=1, alpha=1.0)
        ax_peaks.legend(fontsize=12)
        ax_peaks.set_xlabel('Lead glass charge [arbitrary units]',fontsize=ax_lbl_font_size)
        ax_peaks.set_ylabel('Counts/bin',fontsize=ax_lbl_font_size)
        ax_peaks.tick_params(axis="x", labelsize=ax_tick_font_size)
        ax_peaks.tick_params(axis="y", labelsize=ax_tick_font_size)
        ax_peaks.set_yscale('log')
        gamma_peaks_max = np.max([gamma_peaks_max, np.max(hpeaks[0])*1.1])

        # Plot the histogram (key plots).
        lbl = '{} $\sigma$/$\mu$ = {:.2f}'.format(label, popt[2]/popt[1])
        flat_axes_keyplots[0].hist(arr, bins=nbins, alpha=0.6, label=lbl, range=range)
        flat_axes_keyplots[0].plot(bin_centers, fit_curve, 'r-', alpha=0.6)
        flat_axes_keyplots[0].legend(fontsize=12)
        flat_axes_keyplots[0].set_xlabel('Lead glass charge [arbitrary units]',fontsize=ax_lbl_font_size)
        flat_axes_keyplots[0].set_ylabel('Counts/bin',fontsize=ax_lbl_font_size)
        flat_axes_keyplots[0].tick_params(axis="x", labelsize=ax_tick_font_size)
        flat_axes_keyplots[0].tick_params(axis="y", labelsize=ax_tick_font_size)

    # Set the ylimit for the gamma peaks plot.
    ax_peaks.set_ylim([0.8,gamma_peaks_max])

    # Convert relevant lists to numpy arrays.        
    fit_means = np.array(fit_means)
    fit_smeans = np.array(fit_smeans)
    fit_sigmas = np.array(fit_sigmas)
    fit_ssigmas = np.array(fit_ssigmas)
    # -----------------------------------------------------------------------------------

    # -----------------------------------------------------------------------------------
    # Plot the energy resolution sigma/mean vs. mean
    err_sigma_over_mean = np.sqrt(fit_smeans**2 + fit_ssigmas**2)
    flat_axes_keyplots[1].errorbar(fit_means,fit_sigmas/fit_means,fmt='o',yerr=err_sigma_over_mean)
    flat_axes_keyplots[1].set_xlabel('Lead glass charge [arbitrary units]',fontsize=ax_lbl_font_size)
    flat_axes_keyplots[1].set_ylabel('Gamma peak sigma/mean',fontsize=ax_lbl_font_size)
    flat_axes_keyplots[1].tick_params(axis="x", labelsize=ax_tick_font_size)
    flat_axes_keyplots[1].tick_params(axis="y", labelsize=ax_tick_font_size)
    # -----------------------------------------------------------------------------------
    
    # -----------------------------------------------------------------------------------
    # Plot the lead glass charge vs. expected energy.
    elec_hit_momenta_values = [v for v in elec_hit_momenta.values()]
    e_gamma_expected = [run_momentum - mval*1000 for mval in elec_hit_momenta_values[::-1]]
    flat_axes_keyplots[2].errorbar(e_gamma_expected, fit_means, yerr=fit_smeans, fmt='o')
    flat_axes_keyplots[2].set_ylabel('Lead glass charge [arbitrary units]',fontsize=ax_lbl_font_size)
    flat_axes_keyplots[2].set_xlabel('Expected momentum [MeV/c]',fontsize=ax_lbl_font_size)
    flat_axes_keyplots[2].tick_params(axis="x", labelsize=ax_tick_font_size)
    flat_axes_keyplots[2].tick_params(axis="y", labelsize=ax_tick_font_size)
    #flat_axes[2].set_title("Run {}, p = +{} MeV/c".format(run,run_momentum),fontsize=20)

    # Fit to a line
    p0 = [(e_gamma_expected[-1] - e_gamma_expected[0])/(fit_means[-1] - fit_means[0]),e_gamma_expected[0]]
    popt, pcov = curve_fit(line, e_gamma_expected, fit_means, p0, sigma=fit_smeans)
    x = np.linspace(e_gamma_expected[0], e_gamma_expected[-1], 1000)
    y = line(x, *popt)
    perr = np.sqrt(np.diag(pcov))
    flat_axes_keyplots[2].plot(x, y, label='C = $({:.2f} \pm {:.2f}) x 10^{{-4}} \cdot $p + $({:.3f} \pm {:.3f})$'.format(popt[0]*10000,perr[0]*10000,popt[1],perr[1]), color='red', alpha=0.8, linewidth=3, linestyle=':')
    flat_axes_keyplots[2].legend(loc=2,fontsize=20)
    # -----------------------------------------------------------------------------------

    # -----------------------------------------------------------------------------------
    # Plot the lead glass resolution.

    # Convert gamma peak means and sigmas to MeV
    mconv, merr = popt[0], perr[0]
    bconv, berr = popt[1], perr[1]
    fit_means_MeV = (fit_means - bconv)/mconv
    fit_sigmas_MeV = (fit_sigmas/fit_means)*fit_means_MeV
    fit_ssigmas_MeV = (fit_ssigmas/fit_means)*fit_means_MeV

    # Subtract the momentum ranges from the fit sigmas in quadrature
    momentum_range_values = [v for v in momentum_range.values()]
    momentum_sigmas = np.array([mval for mval in momentum_range_values[::-1]]) #/12**0.5
    lg_sigmas_MeV = (fit_sigmas_MeV**2 - (momentum_sigmas*1000)**2)**0.5

    # Fit the fit_means_MeV vs. lg_sigmas_MeV to a F*sqrt(E) + C
    p0 = [lg_sigmas_MeV[1]/fit_means_MeV[1]**0.5,0.]
    # popt, pcov = curve_fit(resfxn, fit_means_MeV[1:], lg_sigmas_MeV[1:], p0=p0, sigma=fit_ssigmas_MeV[1:])
    # x_res = np.linspace(np.min(fit_means_MeV), np.max(fit_means_MeV), 1000)
    # y_res = resfxn(x_res, *popt)
    # res_err = np.sqrt(np.diag(pcov))
    fres = 0 # popt[0]
    cres = 0 # popt[1]
    ferr = 0 # res_err[0]
    cerr = 0 # res_err[1]

    # Plot the momentum range vs. energy
    #flat_axes[3].plot(fit_means_MeV,lg_sigmas,'o',label='LG resolution',color='blue')
    flat_axes_keyplots[3].plot(fit_means_MeV,momentum_sigmas*1000,'^',label='Computed HD momentum range',color='blue')
    #flat_axes_keyplots[3].errorbar(fit_means_MeV,fit_sigmas_MeV,yerr=fit_ssigmas_MeV,fmt='s',label='$(\sigma/\mu)$ x $\mu$ (calibrated)',color='green')
    #flat_axes_keyplots[3].errorbar(fit_means_MeV,lg_sigmas_MeV,yerr=fit_ssigmas_MeV,fmt='o',label='LG resolution (calibrated)',color='red')
    #flat_axes_keyplots[3].plot(x_res, y_res, label='R = $({:.1f} \pm {:.1f}) \sqrt{{p}} + ({:.1f} \pm {:.1f})$'.format(fres,ferr,cres,cerr), color='red', alpha=0.8, linewidth=3, linestyle=':')
    flat_axes_keyplots[3].set_ylabel('Resolution [MeV]',fontsize=ax_lbl_font_size)
    flat_axes_keyplots[3].set_xlabel('Expected momentum [MeV/c]',fontsize=ax_lbl_font_size)
    flat_axes_keyplots[3].tick_params(axis="x", labelsize=ax_tick_font_size)
    flat_axes_keyplots[3].tick_params(axis="y", labelsize=ax_tick_font_size)
    flat_axes_keyplots[3].legend(fontsize=16)
    # -----------------------------------------------------------------------------------

    # Save all figures.
    fig_keyplots.savefig(f"{out_dir}/key_plots.pdf", bbox_inches='tight')
    plt.close(fig_keyplots)

    fig_peaks.savefig(f"{out_dir}/gamma_peaks.pdf", bbox_inches='tight')
    plt.close(fig_peaks)

    fig_allpeaks.savefig(f"{out_dir}/gamma_peaks_separate.pdf", bbox_inches='tight')
    plt.close(fig_allpeaks)

    return [mconv, merr, bconv, berr, fres, ferr, cres, cerr], [fit_means, fit_smeans, fit_sigmas, fit_ssigmas]
