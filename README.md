# T9BeamTestAna
Analysis package for the T9 beam test


To compile the code type make in the main directory.

## Processing data

The first analysis stage includes processing data files. The output is a root file with anaTree TTRee. It includes, trigger timestamps, pedestals, pedestals standard deviations, signal amplitudes and times.

USAGE: ./bin/waveform_analysis.app -i input_list_file -o output_root_file -c analysis_config.json

  **-i input_list_file**      --> ASCII file with run file names (one file per line). It allows users to run with multiple runs. If you put # before  the file name, it will be "commented out" (not used by the code)
  
  **-o output_root_file**     --> output root file
  
  **-c analysis_config.json** --> json configuration file (an example is located in config/config.json

## Making plots
Navigate to scripts/ directory. 

root -q -b 'MakeDataPlots.C("output_file_from_step1", momentum)'

The script will use int momentum to recognize which cuts it needs to use. The script output is a root file with the same name as the input and "_plots" appended to the end. It contains:
  - reference histograms (they have "ref" in the name and don't change when you apply cuts) 
  - TOF histograms that are used to estimate the fraction of particles.  
