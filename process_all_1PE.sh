for run in 445 438 439 473 556
do
#this is to go from a root_run file to the window analysed and converted to number of PE 

    python python/new_analysis/process_waveform_analysis.py data/root_run_000$run.root config/config_noHodoscope.json peakAnalysed_$run.root

    python python/new_analysis/process_ac_waveform_analysis.py data/root_run_000$run.root peakAnalysed_$run.root -16 45 config/config_ac_noHodoscope.json windowPE_-16ns_45ns_run$run.root

done
