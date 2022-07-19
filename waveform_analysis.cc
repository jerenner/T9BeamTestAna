#include <string>
#include <iostream>
#include <stdlib.h>
#include "InputParser.h"
#include "AnalysisConfig.h"
#include "Utility.h"
#include <TChain.h>
#include "TCanvas.h"
#include "TH1D.h"
#include "WaveformAnalysis.h"


using namespace std;

int main(int argc, char **argv){

    InputParser input(argc, argv);

    if(input.cmdOptionExists("-h")) {
        cerr << "USAGE: " << argv[0] << "-i input_list_file -o output_root_file -c analysis_config.json" << endl;
        exit(EXIT_SUCCESS);
    }



    if(!input.cmdOptionExists("-i")) {
        cerr << "Input list file not provided!" << endl;
        cerr << "USAGE: " << argv[0] << " -i input_list_file -o output_root_file -c analysis_config.json" << endl;
        exit(EXIT_FAILURE);
    }

    if(!input.cmdOptionExists("-o")) {
        cerr << "Output root file not provided!" << endl;
        cerr << "USAGE: " << argv[0] << " -i input_list_file -o output_root_file -c analysis_config.json" << endl;
        exit(EXIT_FAILURE);
    }


    const std::string inputFile = input.getCmdOption("-i");
    if (inputFile.empty()) {
        cerr << "Input list file not provided!" << endl;
        cerr << "USAGE: " << argv[0] << " -i input_list_file -o output_root_file -c analysis_config.json" << endl;
        exit(EXIT_FAILURE);
    }

    const std::string cfgFile = input.getCmdOption("-c");
    if (cfgFile.empty()) {
        cerr << "Analysis config file not provided!" << endl;
        cerr << "USAGE: " << argv[0] << " -i input_list_file -o output_root_file -c analysis_config.json" << endl;
        exit(EXIT_FAILURE);
    }

    const std::string outFile = input.getCmdOption("-o");
    if (outFile.empty()) {
        cerr << "Analysis config file not provided!" << endl;
        cerr << "USAGE: " << argv[0] << " -i input_list_file -o output_root_file -c analysis_config.json" << endl;
        exit(EXIT_FAILURE);
    }

    AnalysisConfig cfg(cfgFile);
    cfg.Print();

    vector<TH1D*> waveforms;
    vector<WaveformAnalysis> waveAna;

    int nChannels = 0;
    for(int i = 0; i < cfg.GetNumberOfChannels(); i++){
        waveforms.push_back(NULL);

        WaveformAnalysis temp;
        if(cfg.IsActive(i)){
            nChannels++;
            temp.SetVoltageScale(cfg.GetVoltageScale());
            temp.SetThreshold(cfg.GetThreshold(i));
            temp.SetPolarity(cfg.GetPolarity(i));
            temp.SetPedestalBinWindow(cfg.GetPedestalBinLow(i), cfg.GetPedestalBinHigh(i));
            temp.SetAnalysisBinWindow(cfg.GetAnalysisBinLow(i), cfg.GetAnalysisBinHigh(i));
            temp.SetChargeMeasurement(cfg.MeasureCharge(i));
            temp.SetTimeMeasurement(cfg.MeasureTime(i));
        }
        waveAna.push_back(temp);
    }

    vector<string> fileNames;
    utl::ReadInputList(inputFile, fileNames);

    //TH1D *waveforms[8] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    int timestamp;
    uint64_t triggerTime1;
    uint64_t triggerTime2;
    int serialnumber;
    int freqsetting;

    TChain chain1("midas_data1");
    TChain chain2("midas_data2");
    for(auto i : fileNames){
        chain1.Add(i.c_str());
        chain2.Add(i.c_str());
    }
    //chain.SetBranchAddress("timestamp",&timestamp);
    //chain.SetBranchAddress("dig1Time",&dig1Time);
    //chain.SetBranchAddress("dig2Time",&dig2Time);
    //chain.SetBranchAddress("serialnumber",&serialnumber);
    chain1.SetBranchAddress("triggerTime",&triggerTime1);
    chain2.SetBranchAddress("triggerTime",&triggerTime2);
    //chain.SetBranchAddress("freqsetting",&freqsetting);


    for(int i=0; i<cfg.GetNumberOfChannels()/2; i++){
        chain1.SetBranchAddress(Form("Channel%d",i),&(waveforms.at(i)));
    }
    for(int i=0; i<cfg.GetNumberOfChannels()/2; i++){
        chain2.SetBranchAddress(Form("Channel%d",i),&(waveforms.at(i+cfg.GetNumberOfChannels()/2)));
    }





    int nent = chain1.GetEntries();
    //nent = 81;

    TFile output(outFile.c_str(), "RECREATE");
    output.cd();
    double pedestal[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    double pedestalSigma[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    //vector<double> pedestal;
    //vector<double> pedestalSigma;
    vector<vector<double>> peakVoltage;
    vector<vector<double>> peakTime;
    vector<vector<double>> signalTime;
    vector<vector<double>> intCharge;
 
    
    TTree *ana_data = new TTree("anaTree", "");

    
    ana_data->Branch("nChannels",&nChannels,"nChannels/I");
    ana_data->Branch("triggerTime1",&triggerTime1,"triggerTime1/I");
    ana_data->Branch("triggerTime2",&triggerTime2,"triggerTime2/I");
    ana_data->Branch("Pedestal",pedestal,"Pedestal[nChannels]/D");
    ana_data->Branch("PedestalSigma", pedestalSigma, "PedestalSigma[nChannels]/D"); 
    ana_data->Branch("PeakVoltage",&peakVoltage);
    ana_data->Branch("PeakTime",&peakTime);
    ana_data->Branch("SignalTime",&signalTime);
    ana_data->Branch("IntCharge",&intCharge);



    for(int i = 0; i < nent; i++){
        //cout << "++++++++++++++++++++++++" << endl;
        chain1.GetEntry(i);
        chain2.GetEntry(i);
        if(!((i+1)%1000))
            cout << "Processing event #: " << i+1 << endl;

        peakVoltage.clear();
        peakTime.clear();
        signalTime.clear();
        intCharge.clear();
        for(int j = 0; j < cfg.GetNumberOfChannels(); j++){

            if(cfg.IsActive(j)){


                waveAna.at(j).SetHistogram(waveforms.at(j));
                waveAna.at(j).RunAnalysis();

                pedestal[j] = waveAna.at(j).GetPedestal();
                pedestalSigma[j] = waveAna.at(j).GetPedestalSigma();
                peakVoltage.push_back(waveAna.at(j).GetPeakVoltage());
                peakTime.push_back(waveAna.at(j).GetPeakTime());
                signalTime.push_back(waveAna.at(j).GetSignalTime());
                intCharge.push_back(waveAna.at(j).GetIntegratedCharge());
                

            }
        }


        ana_data->Fill();
    }
    output.cd();
    ana_data->Write();
    output.Close();
}
