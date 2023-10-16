for run in 445 438 439 473 556
do

#    python python/new_analysis/process_waveform_analysis.py data/root_run_000$run.root config/config_noHodoscope.json singlePEstudy/test_$run.root

    python python/new_analysis/process_ac_waveform_analysis.py data/root_run_000$run.root singlePEstudy/test_$run.root -16 20 config/config_ac_noHodoscope.json singlePEstudy/singlePE_-16ns_20ns_run$run.root

done
