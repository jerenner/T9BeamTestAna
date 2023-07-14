#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSystem.h"

#include <string>
#include <vector>
#include <iostream>

using namespace std;

// Matej Pavin 2022
// modified by Jiri Kvita 2022
// modified for 2023 32 channels July 2023

// ______________________________________________________________
// ______________________________________________________________
// ______________________________________________________________

double GetBeta(double mass, double momentum) {
    double bg = momentum/mass;
    double beta = sqrt(bg*bg/(1+bg*bg));
    return beta;
}


// ______________________________________________________________

void MakeDataPlots(string fileName, int momentum) {
    vector<vector<double> > *peakVoltage = NULL;
    vector<vector<double> > *peakTime = NULL;
    vector<vector<double> > *signalTime = NULL;
    vector<vector<double> > *intCharge = NULL;

    gSystem->Exec("mkdir -p histos/");
    
    const int nChannels = 32;
    double pedestal[nChannels];
    double pedestalSigma[nChannels];

    
    TFile inFile(fileName.c_str(), "READ");
    TTree *tree = (TTree*) inFile.Get("anaTree");

    tree->SetBranchAddress("PeakVoltage",&peakVoltage);
    tree->SetBranchAddress("PeakTime",&peakTime);
    tree->SetBranchAddress("SignalTime",&signalTime);
    tree->SetBranchAddress("IntCharge",&intCharge);
    tree->SetBranchAddress("Pedestal",&pedestal);
    tree->SetBranchAddress("PedestalSigma",&pedestalSigma);
    //tree->SetBranchAddress("PassThreshold",&passThreshold);
    
    int ent = tree->GetEntries();

    double tofmin = 10.;
    double tofmax = 30.;
    
    TH1D hTOF("hRef_TOFAll", "", 250, -100, 150);
    TH1D hTOFAll("hTOFAll", "", 120, tofmin, tofmax);
    TH1D hTOFAllWide("hTOFAllWide", "", 720, tofmin, 2*tofmax);

    TH1D hTOFEl("hTOFEl", "", 120, tofmin, tofmax);
    TH1D hTOFOther("hTOFOther", "", 120, tofmin, tofmax);
    TH1D hT0("hRef_T0", "", 270, 50, 320);
    TH1D hT1("hRef_T1", "", 270, 50, 320);

    // jiri 
    TH1D hTimeReso0("hTimeReso0", "", 200, -100, 100);
    TH1D hTimeReso1("hTimeReso1", "", 200, -100, 100);
    TH1D hTimeReso0_zoom("hTimeReso0_zoom", "", 160, 20, 30);
    TH1D hTimeReso1_zoom("hTimeReso1_zoom", "", 160, -5, 5);
    
    vector<TH1D> hCharge;
    vector<TH1D> hVoltage;
    vector<TH1D> hHit;
    vector<TH1D> hPedestalSigma;
    vector<TH1D> hTime;
    
    // no cuts
    TH2D hTOFACT1V("hRef_TOFACT1V", "; t_{1}-t_{0} [ns]; ACT1 Amplitude", 200, tofmin, tofmax, 200, 0., 1.6);
    TH2D hTOFACT2V("hRef_TOFACT2V", "; t_{1}-t_{0} [ns]; ACT2 Amplitude", 200, tofmin, tofmax, 200, 0., 0.1);
    TH2D hTOFACT3V("hRef_TOFACT3V", "; t_{1}-t_{0} [ns]; ACT3 Amplitude", 200, 37.5, 42.5, 200, 0., 0.1);
        
    TH2D hTOFACT1C("hRef_TOFACT1C", "; t_{1}-t_{0} [ns]; ACT1 Charge", 200, tofmin, tofmax, 200, 0., 0.016);
    TH2D hTOFACT2C("hRef_TOFACT2C", "; t_{1}-t_{0} [ns]; ACT2 Charge", 200, tofmin, tofmax, 200, 0., 0.002);
    TH2D hTOFACT3C("hRef_TOFACT3C", "; t_{1}-t_{0} [ns]; ACT3 Charge", 200, tofmin, tofmax, 200, 0., 0.002);

    // electrons
    TH2D hTOFACT1V_el("hRef_TOFACT1V_el", "; t_{1}-t_{0} [ns]; ACT1 Amplitude", 200, tofmin, tofmax, 200, 0., 1.6);
    TH2D hTOFACT2V_el("hRef_TOFACT2V_el", "; t_{1}-t_{0} [ns]; ACT2 Amplitude", 200, tofmin, tofmax, 200, 0., 0.1);
    TH2D hTOFACT3V_el("hRef_TOFACT3V_el", "; t_{1}-t_{0} [ns]; ACT3 Amplitude", 200, 37.5, 42.5, 200, 0., 0.1);
    
    TH2D hTOFACT1C_el("hRef_TOFACT1C_el", "; t_{1}-t_{0} [ns]; ACT1 Charge", 200, tofmin, tofmax, 200, 0., 0.016);
    TH2D hTOFACT2C_el("hRef_TOFACT2C_el", "; t_{1}-t_{0} [ns]; ACT2 Charge", 200, tofmin, tofmax, 200, 0., 0.002);
    TH2D hTOFACT3C_el("hRef_TOFACT3C_el", "; t_{1}-t_{0} [ns]; ACT3 Charge", 200, tofmin, tofmax, 200, 0., 0.002);

    // ACT2cut
    TH2D hTOFACT1V_act2cut("hRef_TOFACT1V_act2cut", "; t_{1}-t_{0} [ns]; ACT1 Amplitude non-ele", 200, tofmin, tofmax, 200, 0., 1.6);
    TH2D hTOFACT2V_act2cut("hRef_TOFACT2V_act2cut", "; t_{1}-t_{0} [ns]; ACT2 Amplitude", 200, tofmin, tofmax, 200, 0., 0.1);
    TH2D hTOFACT3V_act2cut("hRef_TOFACT3V_act2cut", "; t_{1}-t_{0} [ns]; ACT3 Amplitude", 200, 37.5, 42.5, 200, 0., 0.1);
    
    TH2D hTOFACT1C_act2cut("hRef_TOFACT1C_act2cut", "; t_{1}-t_{0} [ns]; ACT1 Charge", 200, tofmin, tofmax, 200, 0., 0.016);
    TH2D hTOFACT2C_act2cut("hRef_TOFACT2C_act2cut", "; t_{1}-t_{0} [ns]; ACT2 Charge", 200, tofmin, tofmax, 200, 0., 0.002);
    TH2D hTOFACT3C_act2cut("hRef_TOFACT3C_act2cut", "; t_{1}-t_{0} [ns]; ACT3 Charge", 200, tofmin, tofmax, 200, 0., 0.002);

    TH1D hTOF_act2cut("hTOF_act2cut", "; t_{1}-t_{0} [ns];", 120, tofmin, tofmax);

    // ACT3cut
    TH2D hTOFACT1V_act3cut("hRef_TOFACT1V_act3cut", "; t_{1}-t_{0} [ns]; ACT1 Amplitude non-ele", 200, tofmin, tofmax, 200, 0., 1.6);
    TH2D hTOFACT2V_act3cut("hRef_TOFACT2V_act3cut", "; t_{1}-t_{0} [ns]; ACT2 Amplitude", 200, tofmin, tofmax, 200, 0., 0.1);
    TH2D hTOFACT3V_act3cut("hRef_TOFACT3V_act3cut", "; t_{1}-t_{0} [ns]; ACT3 Amplitude", 200, 37.5, 42.5, 200, 0., 0.1);
    
    TH2D hTOFACT1C_act3cut("hRef_TOFACT1C_act3cut", "; t_{1}-t_{0} [ns]; ACT1 Charge", 200, tofmin, tofmax, 200, 0., 0.016);
    TH2D hTOFACT2C_act3cut("hRef_TOFACT2C_act3cut", "; t_{1}-t_{0} [ns]; ACT2 Charge", 200, tofmin, tofmax, 200, 0., 0.002);
    TH2D hTOFACT3C_act3cut("hRef_TOFACT3C_act3cut", "; t_{1}-t_{0} [ns]; ACT3 Charge", 200, tofmin, tofmax, 200, 0., 0.002);

    TH1D hTOF_act3cut("hTOF_act3cut", "; t_{1}-t_{0} [ns];", 120, tofmin, tofmax);

    // 2D ACT charges
    TH2D hACT2CACT1C("hRef_ACT2CACT1C", "; ACT2 Charge; ACT1 Charge", 200, 0., 0.002, 200, 0., 0.016);
    TH2D hACT3CACT2C("hRef_ACT3CACT2C", "; ACT3 Charge; ACT2 Charge", 200, 0., 0.002, 200, 0., 0.002);
    TH2D hACT1CACT3C("hRef_ACT1CACT3C", "; ACT1 Charge; ACT3 Charge", 200, 0., 0.016, 200, 0., 0.002);

    for(int i = 0; i < nChannels; i++) {
        string name1 = "hRef_Charge" + to_string(i);
        string name2 = "hRef_Voltage" + to_string(i);
        string name3 = "hRef_Hits" + to_string(i);
        string name4 = "hRef_PedestalSigma" + to_string(i);
        string name5 = "hRef_Time" + to_string(i);
        
        string title1 = "Channel " + to_string(i) + "; Charge [nC]; Triggers";
        string title2 = "Channel " + to_string(i) + "; Total Amplitude [V]; Triggers";
        string title3 = "Channel " + to_string(i) + "; Hits per trigger; Triggers";
        string title4 = "Channel " + to_string(i) + "; #sigma_{ped} [V]; Triggers";
        string title5 = "Channel " + to_string(i) + "; Time [ns]; Triggers";
        TH1D temp1(name1.c_str(), title1.c_str(), 200, 0., 5*0.08);
        TH1D temp2(name2.c_str(), title2.c_str(), 200, 0., 13*0.8);
        TH1D temp3(name3.c_str(), title3.c_str(), 5, -0.5, 4.5);
        TH1D temp4(name4.c_str(), title4.c_str(), 200, 0., 0.01);
        TH1D temp5(name5.c_str(), title5.c_str(), 200, 250., 650);
        hCharge.push_back(temp1);
        hVoltage.push_back(temp2);
        hHit.push_back(temp3);
        hPedestalSigma.push_back(temp4);
        hTime.push_back(temp5);
    }

    // +-------------------------------+
    // |         event loop            |
    // +-------------------------------+

    cout << "Event loop!" << endl;
    for(int i = 0; i < ent; i++) {

        tree->GetEntry(i);
        vector<int> indices(nChannels, 0);;

        for(int j = 0; j < nChannels; j++) {

            int ind = std::max_element(peakVoltage->at(j).begin(),peakVoltage->at(j).end()) - peakVoltage->at(j).begin();
            indices.at(j) = ind;
            
            hCharge.at(j).Fill(intCharge->at(j).at(indices.at(j)));
            hVoltage.at(j).Fill(peakVoltage->at(j).at(indices.at(j)));
            hTime.at(j).Fill(signalTime->at(j).at(indices.at(j)));
            hHit.at(j).Fill(peakVoltage->at(j).size());
            hPedestalSigma.at(j).Fill(pedestalSigma[j]);   
        }


	// JK's time resolution
        double t0a = (signalTime->at(8).at(indices.at(8))   + signalTime->at(11).at(indices.at(11)))/2.;
        double t0b = (signalTime->at(9).at(indices.at(9)) + signalTime->at(10).at(indices.at(10)))/2.;
        double t1a = (signalTime->at(12).at(indices.at(12)) + signalTime->at(15).at(indices.at(15)))/2.;
        double t1b = (signalTime->at(13).at(indices.at(13)) + signalTime->at(14).at(indices.at(14)))/2.;
	
	// time diffs for time resolution histogramme:
	double t0diff = t0a - t0b;
	double t1diff = t1a - t1b;

	// jiri
	// Fill resolution histograms:
        hTimeReso0.Fill(t0diff);
	hTimeReso1.Fill(t1diff);
	hTimeReso0_zoom.Fill(t0diff);
	hTimeReso1_zoom.Fill(t1diff);
	  
        double t0 = (signalTime->at(8).at(indices.at(8)) + signalTime->at(9).at(indices.at(9)) + signalTime->at(10).at(indices.at(10)) + signalTime->at(11).at(indices.at(11)))/4.;
        double t1 = (signalTime->at(12).at(indices.at(12)) + signalTime->at(13).at(indices.at(13)) + signalTime->at(14).at(indices.at(14)) + signalTime->at(15).at(indices.at(15)))/4.;
	double tof = t1-t0;

	double act1v = peakVoltage->at(0).at(indices.at(0)) + peakVoltage->at(1).at(indices.at(1));
	double act2v = peakVoltage->at(2).at(indices.at(2)) + peakVoltage->at(3).at(indices.at(3));
	double act3v = peakVoltage->at(4).at(indices.at(4)) + peakVoltage->at(5).at(indices.at(5));
        hTOFACT1V.Fill(tof, act1v);
        hTOFACT2V.Fill(tof, act2v);
        hTOFACT3V.Fill(tof, act3v);

	double act1c = intCharge->at(0).at(indices.at(0)) + intCharge->at(1).at(indices.at(1));
	double act2c = intCharge->at(2).at(indices.at(2)) + intCharge->at(3).at(indices.at(3));
	double act3c = intCharge->at(4).at(indices.at(4)) + intCharge->at(5).at(indices.at(5));
        hTOFACT1C.Fill(tof, act1c);
        hTOFACT2C.Fill(tof, act2c);
        hTOFACT3C.Fill(tof, act3c);

        hACT1CACT3C.Fill(act1c, act3c);
        hACT3CACT2C.Fill(act3c, act2c);
        hACT2CACT1C.Fill(act2c, act1c); 
	
        hTOF.Fill(tof);
        hT0.Fill(t0);
        hT1.Fill(t1);
	
        bool pass = true;
        bool isEl = false;
	
        switch(momentum)	  {

	  // just a code example from 2022
	  case 320: { //; jk
	    if (act1v > 0.18)
	      isEl = true;
	    if (act2v > 0.025)
	      isEl = true; 
	    if (act3v > 0.022)
	      isEl = true;
	    break;
	  } // 320

	 
	    /*
	  case 400: {
	    double voltageCut[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025, 0.025, 0.025, 0.025, 0.025, 0.040, 0.025, 0.025};
	    for(int j = 0; j < 16; j++) {
	      if (peakVoltage->at(j).at(indices.at(j)) < voltageCut[j]) {
		pass = false;
		break;
	      }
	    }
	    break;
	  } // 400
	    */
	        
	  default: {
	    if (i < 10)
	      cout << "WARNING: Using default settings for the " << momentum << " MeV/c beam" << endl;
	    if (act1v > 0.18)
	      isEl = true;
	    if (act2v > 0.025)
	      isEl = true; 
	    if (act3v > 0.025)
	      isEl = true;
	    break;
	  }
	  
	  } // case

        if (!pass) continue;
	
      
        if (isEl) {
	  // electrons
	  hTOFEl.Fill(tof);

	  hTOFACT1V_el.Fill(tof, act1v);
	  hTOFACT2V_el.Fill(tof, act2v);
	  hTOFACT3V_el.Fill(tof, act3v);
	      
	  hTOFACT1C_el.Fill(tof, act1c);
	  hTOFACT2C_el.Fill(tof, act2c);
	  hTOFACT3C_el.Fill(tof, act3c);
        }
        else {

	  // non-electrons
	  hTOFOther.Fill(tof);
	
	 
        } // non-electrons

	if (act1v < 0.32) { // custom electron removal cut
	
	  hTOFAll.Fill(tof);
	  hTOFAllWide.Fill(tof);

	
	  if (act2v > 0.0080) { // was: 0.0057, 0.0060, tried also 0.0055
	    hTOFACT1V_act2cut.Fill(tof, act1v);
	    hTOFACT2V_act2cut.Fill(tof, act2v);
	    hTOFACT3V_act2cut.Fill(tof, act3v);
	      
	    hTOFACT1C_act2cut.Fill(tof, act1c);
	    hTOFACT2C_act2cut.Fill(tof, act2c);
	    hTOFACT3C_act2cut.Fill(tof, act3c);

	    hTOF_act2cut.Fill(tof);
	  } // act2 cuts

	  // TO VALIDATE the cut values!!!
	  if (act3v > 0.0080) { // was: 0.0060
	    hTOFACT1V_act3cut.Fill(tof, act1v);
	    hTOFACT2V_act3cut.Fill(tof, act2v);
	    hTOFACT3V_act3cut.Fill(tof, act3v);
	      
	    hTOFACT1C_act3cut.Fill(tof, act1c);
	    hTOFACT2C_act3cut.Fill(tof, act2c);
	    hTOFACT3C_act3cut.Fill(tof, act3c);

	    hTOF_act3cut.Fill(tof);
	  } // act3 cuts

	} // custom electron removal cut


	
    } // entries
    cout << "End of event loop!" << endl;

    
    TString outFileName = fileName.substr(0, fileName.size()-5) + "_plots.root";
    outFileName = outFileName.ReplaceAll("output/", "histos/");
    
    TFile outFile(outFileName.Data(), "RECREATE");
    outFile.cd();
    

    hTOFACT1V.Write();
    hTOFACT2V.Write();
    hTOFACT3V.Write();
    hTOFACT1C.Write();
    hTOFACT2C.Write();
    hTOFACT3C.Write();

    hTOFACT1V_el.Write();
    hTOFACT2V_el.Write();
    hTOFACT3V_el.Write();
    hTOFACT1C_el.Write();
    hTOFACT2C_el.Write();
    hTOFACT3C_el.Write();
    
    hTOFACT1V_act2cut.Write();
    hTOFACT2V_act2cut.Write();
    hTOFACT3V_act2cut.Write();
    hTOFACT1C_act2cut.Write();
    hTOFACT2C_act2cut.Write();
    hTOFACT3C_act2cut.Write();


    hTOFACT1V_act3cut.Write();
    hTOFACT2V_act3cut.Write();
    hTOFACT3V_act3cut.Write();
    hTOFACT1C_act3cut.Write();
    hTOFACT2C_act3cut.Write();
    hTOFACT3C_act3cut.Write();
    
    
    hTOF.Write();
    hT0.Write();
    hT1.Write();

    hTimeReso0.Write();
    hTimeReso1.Write();
    hTimeReso0_zoom.Write();
    hTimeReso1_zoom.Write();

    for (auto hist: hVoltage) hist.Write();
    for (auto hist: hCharge) hist.Write();
    for (auto hist: hHit) hist.Write();
    for (auto hist: hPedestalSigma) hist.Write();
    for (auto hist: hTime) hist.Write();

    hACT1CACT3C.Write();
    hACT3CACT2C.Write();
    hACT2CACT1C.Write();

    hTOFAll.Write();
    hTOFAllWide.Write();
    hTOFEl.Write();
    hTOFOther.Write();

    hTOF_act2cut.Write();
    hTOF_act3cut.Write();


}


