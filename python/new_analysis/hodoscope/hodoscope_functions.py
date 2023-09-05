import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

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

def plot_histograms_for_each_signal(df_dict, base_dir=".", rnum = 0, quantity='nPeaks', select_nonzero_peaks=False, logscale=False, nbins=60):
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
        
        # Select only non-zero peaks if specified.
        if(select_nonzero_peaks):
            df_select = df[df['nPeaks'] > 0]
        else:
            df_select = df
            
        # Place a cut if we're looking at an HD element and plotting PeakVoltage or IntCharge.
        if(key[0:2] == 'HD' and (quantity == 'IntCharge' or quantity == 'PeakVoltage')):
            df_select = df[df[quantity] > 0.02]

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

def filter_range(name,df,key,rng,drop=False):
    
    # Perform the filter.
    df_filtered = df[df[key].between(*rng)]
    nduplicates = get_num_duplicates(df_filtered['event'].values)
    print(f"{name} filter passed",len(df_filtered),"of",len(df),"with",nduplicates,"duplicates")
    
    # Handle duplicates.
    if(nduplicates > 0):
        
        # Drop all events that appear more than once if selected.
        if(drop):
            df_filtered = df_filtered[~df_filtered['event'].duplicated(keep=False)]
        # Otherwise choose the one with the maximum integrated charge.
        else:
            idx = df_filtered.groupby('event')['IntCharge'].idxmax()
            df_filtered = df_filtered.loc[idx]
        
    print("--> returning",len(df_filtered),"events")
    return df_filtered

def timing_analysis(df_dict, pb_timing_range, tof0_timing_range, tof0_charge_range, 
                    tof1_timing_range, tof1_charge_range, t2_timing_range, t2_charge_range,
                    act0_timing_range, act0_charge_range, act1_timing_range, act1_charge_range,
                    act3_timing_range, hd_timing_ranges, hd_charge_ranges, low_radiation = False):
    
    # Filter lead glass
    pb_filtered = filter_range("PbGlass",df_dict['PbGlass'],'PeakTime',pb_timing_range)

    # Filter TOF elements
    # --------------------------------------------------------------------------------------------
    # TOF0
    tof00_filtered = filter_range("TOF00",df_dict['TOF00'],'PeakTime',tof0_timing_range)
    tof01_filtered = filter_range("TOF01",df_dict['TOF01'],'PeakTime',tof0_timing_range)
    tof02_filtered = filter_range("TOF02",df_dict['TOF02'],'PeakTime',tof0_timing_range)
    tof03_filtered = filter_range("TOF03",df_dict['TOF03'],'PeakTime',tof0_timing_range)

    combined_00_01 = tof00_filtered.merge(tof01_filtered, on='event', suffixes=('_00', '_01'))
    combined_00_01_02 = combined_00_01.merge(tof02_filtered, on='event')
    combined_00_01_02 = combined_00_01_02.rename(columns={'IntCharge': 'IntCharge_02'})
    tof0_combined = combined_00_01_02.merge(tof03_filtered, on='event')
    tof0_combined = tof0_combined.rename(columns={'IntCharge': 'IntCharge_03'})
    tof0_combined['combined_charge'] = (tof0_combined['IntCharge_00'] 
                                    + tof0_combined['IntCharge_01']
                                    + tof0_combined['IntCharge_02']
                                    + tof0_combined['IntCharge_03'])

    tof0_valid = filter_range("TOF0_combined",tof0_combined,'combined_charge',tof0_charge_range)
    tof0_valid.loc[:,'hit_TOF0'] = 1
    print()

    # TOF1
    tof10_filtered = filter_range("TOF10",df_dict['TOF10'],'PeakTime',tof1_timing_range)
    tof11_filtered = filter_range("TOF11",df_dict['TOF11'],'PeakTime',tof1_timing_range)
    tof12_filtered = filter_range("TOF12",df_dict['TOF12'],'PeakTime',tof1_timing_range)
    tof13_filtered = filter_range("TOF13",df_dict['TOF13'],'PeakTime',tof1_timing_range)

    combined_10_11 = tof10_filtered.merge(tof11_filtered, on='event', suffixes=('_10', '_11'))
    combined_10_11_12 = combined_10_11.merge(tof12_filtered, on='event')
    combined_10_11_12 = combined_10_11_12.rename(columns={'IntCharge': 'IntCharge_12'})
    tof1_combined = combined_10_11_12.merge(tof13_filtered, on='event')
    tof1_combined = tof1_combined.rename(columns={'IntCharge': 'IntCharge_13'})
    tof1_combined['combined_charge'] = (tof1_combined['IntCharge_10'] 
                                    + tof1_combined['IntCharge_11']
                                    + tof1_combined['IntCharge_12']
                                    + tof1_combined['IntCharge_13'])

    tof1_valid = filter_range("TOF1_combined",tof1_combined,'combined_charge',tof1_charge_range)
    if(not low_radiation): tof1_valid.loc[:,'hit_TOF1'] = 1
    print()

    # T2
    t2_filtered = filter_range("T2",df_dict['TriggerScint'],'PeakTime',t2_timing_range)
    t2_valid = filter_range("T2",t2_filtered,'IntCharge',t2_charge_range)
    t2_valid.loc[:,'hit_T2'] = 1
    print()
    # --------------------------------------------------------------------------------------------

    # ACT elements
    # --------------------------------------------------------------------------------------------
    # ACT0
    act0l_filtered = filter_range("ACT0L",df_dict['ACT0L'],'PeakTime',act0_timing_range)
    act0r_filtered = filter_range("ACT0R",df_dict['ACT0R'],'PeakTime',act0_timing_range)

    act0_combined = act0l_filtered.merge(act0r_filtered, on='event', suffixes=('_L', '_R'))
    act0_combined['combined_charge'] = act0_combined['IntCharge_L'] + act0_combined['IntCharge_R']

    act0_valid = filter_range("ACT0_combined",act0_combined,'combined_charge',act0_charge_range)
    if(not low_radiation): act0_valid.loc[:,'hit_ACT0'] = 1

    # ACT1
    act1l_filtered = filter_range("ACT1L",df_dict['ACT1L'],'PeakTime',act1_timing_range)
    act1r_filtered = filter_range("ACT1R",df_dict['ACT1R'],'PeakTime',act1_timing_range)

    act1_combined = act1l_filtered.merge(act1r_filtered, on='event', suffixes=('_L', '_R'))
    act1_combined['combined_charge'] = act1_combined['IntCharge_L'] + act1_combined['IntCharge_R']

    act1_valid = filter_range("ACT1_combined",act1_combined,'combined_charge',act1_charge_range)
    act1_valid.loc[:,'hit_ACT1'] = 1

    # ACT3
    act3l_filtered = df_dict['ACT3L'][df_dict['ACT3L']['nPeaks'] == 0]
    act3r_filtered = df_dict['ACT3R'][df_dict['ACT3R']['nPeaks'] == 0]

    act3_combined = act3l_filtered.merge(act3r_filtered, on='event', suffixes=('_L', '_R'))
    act3_valid = act3_combined.copy()
    act3_valid.loc[:,'nohit_ACT3'] = 1
    print("ACT3 total number of valid events:",len(act3_valid))
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
        filtered = df_dict[hd_key][(df_dict[hd_key]['PeakTime'].between(*hd_time_range)) & 
                                (df_dict[hd_key]['IntCharge'].between(*hd_charge_range))]
        print(f"HD{i}: {len(filtered)} of {len(df_dict[hd_key])} events after filter")
        
        # Assign a binary value indicating a hit
        if not filtered.empty:
            filtered.loc[:,hit_col_name] = 1
        else:
            filtered.loc[:,hit_col_name] = 0

        hd_dfs[hd_key] = filtered

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

    u, c = np.unique(combined_hd_df['event'].values, return_counts=True)
    print("Duplicates in combined HD dataframe:",np.sum(c[c > 1]))
    # --------------------------------------------------------------------------------------------

    # Merge all relevant dataframes
    final_df = pb_filtered.merge(hd_valid_events, on='event', how='inner')
    if(not low_radiation): final_df = final_df.merge(act0_valid[['event', 'hit_ACT0']], on='event', how='left')
    final_df = final_df.merge(act1_valid[['event', 'hit_ACT1']], on='event', how='left')
    final_df = final_df.merge(act3_valid[['event', 'nohit_ACT3']], on='event', how='left')
    final_df = final_df.merge(tof0_valid[['event', 'hit_TOF0']], on='event', how='left')
    if(not low_radiation): final_df = final_df.merge(tof1_valid[['event', 'hit_TOF1']], on='event', how='left')
    final_df = final_df.merge(t2_valid[['event', 'hit_T2']], on='event', how='left')

    # Fill NaN values in the hit columns with 0
    final_df[['hit_ACT1', 'nohit_ACT3', 'hit_TOF0', 'hit_T2']] = final_df[['hit_ACT1', 'nohit_ACT3', 'hit_TOF0', 'hit_T2']].fillna(0)

    return final_df

def gaussian(x, amplitude, mean, stddev):
        return amplitude * norm.pdf(x, loc=mean, scale=stddev)

def gamma_peak_plot(final_df, run, nbins, range, timing_cuts = False, low_radiation = False):
    """
    Analyze the gamma peaks.
    """

    # Apply timing cuts if specified.
    if(timing_cuts):

        if(low_radiation):
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
    chg_hd14 = final_df[(final_df.hit_HD14 == 1) & cuts]['IntCharge'].values
    chg_hd13 = final_df[(final_df.hit_HD13 == 1) & cuts]['IntCharge'].values
    chg_hd12 = final_df[(final_df.hit_HD12 == 1) & cuts]['IntCharge'].values
    chg_hd11 = final_df[(final_df.hit_HD11 == 1) & cuts]['IntCharge'].values
    chg_hd10 = final_df[(final_df.hit_HD10 == 1) & cuts]['IntCharge'].values
    chg_hd9 = final_df[(final_df.hit_HD9 == 1) & cuts]['IntCharge'].values
    chg_hd8 = final_df[(final_df.hit_HD8 == 1) & cuts]['IntCharge'].values
    chg_hd7 = final_df[(final_df.hit_HD7 == 1) & cuts]['IntCharge'].values
    chg_hd6 = final_df[(final_df.hit_HD6 == 1) & cuts]['IntCharge'].values
    chg_hd5 = final_df[(final_df.hit_HD5 == 1) & cuts]['IntCharge'].values
    chg_hd4 = final_df[(final_df.hit_HD4 == 1) & cuts]['IntCharge'].values
    chg_hd3 = final_df[(final_df.hit_HD3 == 1) & cuts]['IntCharge'].values
    chg_hd2 = final_df[(final_df.hit_HD2 == 1) & cuts]['IntCharge'].values
    chg_hd1 = final_df[(final_df.hit_HD1 == 1) & cuts]['IntCharge'].values
    chg_hd0 = final_df[(final_df.hit_HD0 == 1) & cuts]['IntCharge'].values

    chg_arrays = [chg_hd0, chg_hd1, chg_hd2, chg_hd3, chg_hd4, chg_hd5, chg_hd6, chg_hd7, chg_hd8,
                  chg_hd9, chg_hd10, chg_hd11, chg_hd12, chg_hd13, chg_hd14]
    labels = ["H0", "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12", "H13", "H14"]

    fig = plt.figure(figsize=(14,7))
    initial_params = [1000, np.mean(chg_hd14), np.std(chg_hd14)]
    for arr,label in zip(chg_arrays[::-1],labels[::-1]):

        # Make the histogram for this hodoscope element
        hist, bin_edges = np.histogram(arr, bins=nbins, range=range)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Gaussian fit
        print("Fitting gamma peak",label,"with initial params",*initial_params)
        popt, _ = curve_fit(gaussian, bin_centers, hist, p0=initial_params)
        fit_curve = gaussian(bin_centers, *popt)
        print("-- fit result",*popt)
        initial_params = popt

        # Plot the histogram
        lbl = '{} $\sigma$/$\mu$ = {:.2f}'.format(label, popt[2]/popt[1])
        plt.hist(arr, bins=nbins, alpha=0.6, label=lbl, range=range);
        plt.plot(bin_centers, fit_curve, 'r-', alpha=0.6)

        plt.legend()
        
        plt.xlabel('Lead Glass (charge)',fontsize=20);
        plt.ylabel('Counts/bin',fontsize=20);
        plt.title("Run {}, p = +500 MeV/c".format(run),fontsize=20);


    #plt.hist(final_df[(final_df.hit_HD14 == 1) & (final_df.total_hits == 1)]['IntCharge'],bins=nbins,label="HD14, no timing cuts",range=[0,0.5],density=False)
    # plt.hist(chg_hd14,bins=nbins,label="HD14, timing cuts",density=True,range=[0,0.5])
    #print("H14 with",len(final_df[(final_df.hit_HD14 == 1) & cuts]),"events",
    #    "and",len(final_df[(final_df.hit_HD14 == 1) & cuts & (final_df.IntCharge < 0.1) & (final_df.IntCharge > 0.01)]),"tail events")
    # plt.hist(,bins=nbins,label="HD13")
    # plt.hist(final_df[(final_df.hit_HD12 == 1) & cuts]['IntCharge'],bins=nbins,label="HD12")
    # plt.hist(final_df[(final_df.hit_HD11 == 1) & cuts]['IntCharge'],bins=nbins,label="HD11")
    # plt.hist(final_df[(final_df.hit_HD10 == 1) & cuts]['IntCharge'],bins=nbins,label="HD10")
    # plt.hist(final_df[(final_df.hit_HD9 == 1) & cuts]['IntCharge'],bins=nbins,label="HD9")
    # plt.hist(final_df[(final_df.hit_HD8 == 1) & cuts]['IntCharge'],bins=nbins,label="HD8")
    # plt.hist(final_df[(final_df.hit_HD7 == 1) & cuts]['IntCharge'],bins=nbins,label="HD7")
    # plt.hist(final_df[(final_df.hit_HD6 == 1) & cuts]['IntCharge'],bins=nbins,label="HD6")
    # plt.hist(final_df[(final_df.hit_HD5 == 1) & cuts]['IntCharge'],bins=nbins,label="HD5")
    # plt.hist(final_df[(final_df.hit_HD4 == 1) & cuts]['IntCharge'],bins=nbins,label="HD4")
    # plt.hist(final_df[(final_df.hit_HD3 == 1) & cuts]['IntCharge'],bins=nbins,label="HD3")
    # plt.hist(final_df[(final_df.hit_HD2 == 1) & cuts]['IntCharge'],bins=nbins,label="HD2")
    # plt.hist(final_df[(final_df.hit_HD1 == 1) & cuts]['IntCharge'],bins=nbins,label="HD1")
    # plt.hist(final_df[(final_df.hit_HD0 == 1) & cuts]['IntCharge'],bins=nbins,label="HD0")
    # plt.xlabel("Lead glass charge",fontsize=14)
    # plt.ylabel("Counts/bin",fontsize=14)
    # plt.title("RUN 000735, p = + 0.8 GeV/c",fontsize=20)
    # plt.legend()

