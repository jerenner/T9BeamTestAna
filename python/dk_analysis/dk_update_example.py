from Run import Run
#from Event import Event
import json5
import numpy as np
import sys

run503 = Run.open_file("C:/Users/Karlen/Documents/temp/run_000503.dk")

xxx=1

# add information about the time difference for 2 peak waveforms using a simple 2 peak finder

#n_peaks = run474.analyses[2].n_peaks

#run474.user['n_peaks'] = n_peaks

run474.save_file("C:/Users/Karlen/Documents/temp/run474_update.dk")