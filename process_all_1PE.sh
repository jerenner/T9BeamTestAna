for run in 464 #502
do
#this is to go from a root_run file to the window analysed and converted to number of PE 

    python python/new_analysis/process_waveform_analysis.py data/root_run_000$run.root config/config_noHodoscope.json singlePEstudy/test_$run.root

    python python/new_analysis/process_ac_waveform_analysis.py data/root_run_000$run.root singlePEstudy/test_$run.root -16 45 config/config_ac_noHodoscope.json singlePEstudy/singlePE_-16ns_45ns_run$run.root

done
