#this is the final version of particle ID, baed on window integration and coicidence requirements. Previous iterarions included particleID and checkScintillation ID. Should be read only with the corresponding config file (no need for anything else)

#Step 0: initialise, reading the user inputs
# a config file that has all of the information
 with open(config_filename) as file:
        config = json5.load(file)["WaveAnalysis"]

particleID_analysis = LowMomentumAnalysis(config)
particleID_analysis.openDataFile()


#Step 1: Select only events with exactly one coincidence
particleID_analysis.selectEvents()

#Step 2: dE/dx in trigger scintillator-based selection of 1 particle events (also using TOF)

#Step 3: ID events using cuts on ACT signal

#Step 4: Plots

#Step 5: Outputs
