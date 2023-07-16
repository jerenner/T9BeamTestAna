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
    
    const int nChannels = 19; // UPDATE THIS FOR HODOSCOPE PMTs! to e.g. 32!
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


    TString outFileName = fileName.substr(0, fileName.size()-5) + "_plots.root";
    outFileName = outFileName.ReplaceAll("output/", "histos/");
    TFile outFile(outFileName.Data(), "RECREATE");
    outFile.cd();
    
    double tofmin = 8.;
    double tofmax = 30.;
    double tofmaxlow = 20.;
    int ntofbins = 110;
    int ntofbins2d = 120;

    double actChargeMax = 0.40;
    double actAmplitudeMax =  10.;

    TH1D hTOFAll("hTOFAll", ";t_{TOF}^{All} [ns]", 120, tofmin, tofmax);
    TH1D hTOFAllWide("hTOFAllWide", ";t_{TOF}^{All} [ns]", 2*ntofbins, tofmin, 2*tofmax);
    TH1D hTOFEl("hTOFEl", ";t_{TOF}^{e} [ns]", ntofbins, tofmin, tofmax);
    TH1D hTOFOther("hTOFOther", ";t_{TOF}^{non-e} [ns]", ntofbins, tofmin, tofmax);

    TH1D hTOFAllLow("hTOFAllLow", ";t_{TOF}^{All} [ns]", 120, tofmin, tofmaxlow);
    TH1D hTOFElLow("hTOFElLow", ";t_{TOF}^{e} [ns]", ntofbins, tofmin, tofmax);
    TH1D hTOFOtherLow("hTOFOtherLow", ";t_{TOF}^{non-e} [ns]", ntofbins, tofmin, tofmax);


    TH1D hT0("hRef_T0", "", 270, 50, 320);
    TH1D hT1("hRef_T1", "", 270, 50, 320);

    // jiri 
    TH1D hTimeReso0("hTimeReso0", "", 200, -100, 100);
    TH1D hTimeReso1("hTimeReso1", "", 200, -100, 100);
    TH1D hTimeReso0_zoom("hTimeReso0_zoom", "", 160, 20, 30);
    TH1D hTimeReso1_zoom("hTimeReso1_zoom", "", 160, -5, 5);

    // 2023 time offset analysis
    TH1D hTimeDiffTOF01("hTimeDiffTOF01", "hTimeDiffTOF01", 100, -12.,12.);
    TH1D hTimeDiffTOF02("hTimeDiffTOF02", "hTimeDiffTOF02", 100, -12.,12.);
    TH1D hTimeDiffTOF03("hTimeDiffTOF03", "hTimeDiffTOF03", 100, -12.,12.);
    
    TH1D hTimeDiffTOF11("hTimeDiffTOF11", "hTimeDiffTOF11", 100, -12.,12.);
    TH1D hTimeDiffTOF12("hTimeDiffTOF12", "hTimeDiffTOF12", 100, -12.,12.);
    TH1D hTimeDiffTOF13("hTimeDiffTOF13", "hTimeDiffTOF13", 100, -12.,12.);
    

    //acraplet TOF analysis
    TH1D hTimeTOF0("hTimeTOF0", "hTimeTOF0", 100, 0.,50.);
    TH1D hTimeTOF1("hTimeTOF1", "hTimeTOF1", 100, 0.,50.);
    TH1D hTimeTOF2("hTimeTOF2", "hTimeTOF2", 100, 0.,50.);
    TH1D hTimeTOF3("hTimeTOF3", "hTimeTOF3", 100, 0.,50.);
    


    // standard
    vector<TH1D> hCharge;
    vector<TH1D> hVoltage;
    vector<TH1D> hHit;
    vector<TH1D> hPedestalSigma;
    vector<TH1D> hTime;
    
    // no cuts
    TH2D hTOFACT0A("hRef_TOFACT0A", "; t_{1}-t_{0} [ns]; ACT0 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
    TH2D hTOFACT1A("hRef_TOFACT1A", "; t_{1}-t_{0} [ns]; ACT1 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
    TH2D hTOFACT2A("hRef_TOFACT2A", "; t_{1}-t_{0} [ns]; ACT2 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
    TH2D hTOFACT3A("hRef_TOFACT3A", "; t_{1}-t_{0} [ns]; ACT3 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
        
    TH2D hTOFACT0C("hRef_TOFACT0C", "; t_{1}-t_{0} [ns]; ACT0 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);
    TH2D hTOFACT1C("hRef_TOFACT1C", "; t_{1}-t_{0} [ns]; ACT1 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);
    TH2D hTOFACT2C("hRef_TOFACT2C", "; t_{1}-t_{0} [ns]; ACT2 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);
    TH2D hTOFACT3C("hRef_TOFACT3C", "; t_{1}-t_{0} [ns]; ACT3 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);

    // electrons
    TH2D hTOFACT0A_el("hRef_TOFACT0A_el", "; t_{1}-t_{0} [ns]; ACT0 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
    TH2D hTOFACT1A_el("hRef_TOFACT1A_el", "; t_{1}-t_{0} [ns]; ACT1 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
    TH2D hTOFACT2A_el("hRef_TOFACT2A_el", "; t_{1}-t_{0} [ns]; ACT2 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
    TH2D hTOFACT3A_el("hRef_TOFACT3A_el", "; t_{1}-t_{0} [ns]; ACT3 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);

    TH2D hTOFACT0C_el("hRef_TOFACT0C_el", "; t_{1}-t_{0} [ns]; ACT0 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);
    TH2D hTOFACT1C_el("hRef_TOFACT1C_el", "; t_{1}-t_{0} [ns]; ACT1 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);
    TH2D hTOFACT2C_el("hRef_TOFACT2C_el", "; t_{1}-t_{0} [ns]; ACT2 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);
    TH2D hTOFACT3C_el("hRef_TOFACT3C_el", "; t_{1}-t_{0} [ns]; ACT3 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);

    // Act1cut
    TH2D hTOFACT0A_act1cut("hRef_TOFACT0A_act1cut", "; t_{1}-t_{0} [ns]; ACT0 Amplitude non-ele", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
    TH2D hTOFACT1A_act1cut("hRef_TOFACT1A_act1cut", "; t_{1}-t_{0} [ns]; ACT1 Amplitude non-ele", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
    TH2D hTOFACT2A_act1cut("hRef_TOFACT2A_act1cut", "; t_{1}-t_{0} [ns]; ACT2 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
    TH2D hTOFACT3A_act1cut("hRef_TOFACT3A_act1cut", "; t_{1}-t_{0} [ns]; ACT3 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);

    TH2D hTOFACT0C_act1cut("hRef_TOFACT0C_act1cut", "; t_{1}-t_{0} [ns]; ACT0 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);
    TH2D hTOFACT1C_act1cut("hRef_TOFACT1C_act1cut", "; t_{1}-t_{0} [ns]; ACT1 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);
    TH2D hTOFACT2C_act1cut("hRef_TOFACT2C_act1cut", "; t_{1}-t_{0} [ns]; ACT2 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);
    TH2D hTOFACT3C_act1cut("hRef_TOFACT3C_act1cut", "; t_{1}-t_{0} [ns]; ACT3 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);

    TH1D hTOF_act1cut("hTOF_act1cut", "; t_{1}-t_{0} [ns];", ntofbins, tofmin, tofmax);

    // Act2cut
    TH2D hTOFACT0A_act2cut("hRef_TOFACT0A_act2cut", "; t_{1}-t_{0} [ns]; ACT0 Amplitude non-ele", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
    TH2D hTOFACT1A_act2cut("hRef_TOFACT1A_act2cut", "; t_{1}-t_{0} [ns]; ACT1 Amplitude non-ele", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
    TH2D hTOFACT2A_act2cut("hRef_TOFACT2A_act2cut", "; t_{1}-t_{0} [ns]; ACT2 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
    TH2D hTOFACT3A_act2cut("hRef_TOFACT3A_act2cut", "; t_{1}-t_{0} [ns]; ACT3 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);

    TH2D hTOFACT0C_act2cut("hRef_TOFACT0C_act2cut", "; t_{1}-t_{0} [ns]; ACT0 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);
    TH2D hTOFACT1C_act2cut("hRef_TOFACT1C_act2cut", "; t_{1}-t_{0} [ns]; ACT1 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);
    TH2D hTOFACT2C_act2cut("hRef_TOFACT2C_act2cut", "; t_{1}-t_{0} [ns]; ACT2 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);
    TH2D hTOFACT3C_act2cut("hRef_TOFACT3C_act2cut", "; t_{1}-t_{0} [ns]; ACT3 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);

    TH1D hTOF_act2cut("hTOF_act2cut", "; t_{1}-t_{0} [ns];", 120, tofmin, tofmax);

    // 2D ACT charges
    TH2D hACT2CACT1C("hRef_ACT2CACT1C", "; ACT2 Charge; ACT1 Charge", 200, 0., actChargeMax, 200, 0., actChargeMax);
    TH2D hACT3CACT2C("hRef_ACT3CACT2C", "; ACT3 Charge; ACT2 Charge", 200, 0., actChargeMax, 200, 0., actChargeMax);
    TH2D hACT1CACT3C("hRef_ACT1CACT3C", "; ACT1 Charge; ACT3 Charge", 200, 0., actChargeMax, 200, 0., actChargeMax);

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
        TH1D temp5(name5.c_str(), title5.c_str(), 270, 0., 540.);
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


        double t00 = signalTime->at(8).at(indices.at(8));
	double t01 = signalTime->at(9).at(indices.at(9));
	double t02 = signalTime->at(10).at(indices.at(10));
	double t03 = signalTime->at(11).at(indices.at(11));
        hTimeDiffTOF01.Fill(t00 - t01); //carefull, this particular histogram is a difference of hit times,
	hTimeDiffTOF02.Fill(t00 - t02); //not a TOF per say but instead the cable length difference 
	hTimeDiffTOF03.Fill(t00 - t03); //plus the photon travel time though the panel

        double t10 = signalTime->at(12).at(indices.at(12));
        double t11 = signalTime->at(13).at(indices.at(13));
	double t12 = signalTime->at(14).at(indices.at(14));
	double t13 = signalTime->at(15).at(indices.at(15));
        hTimeDiffTOF11.Fill(t11 - t10);
	hTimeDiffTOF12.Fill(t11 - t12);
	hTimeDiffTOF13.Fill(t11 - t13);

	//acraplet
	//compare the hit times for the same event recorded by trigger PMTs on the same side (up/down, left/right)
	//assumption: the light travel time trough the trigger counter to the PMT should be about the same (if the beam is well aligned)
	//then we can check if the TOF is constant for a run
	//idea: using tof dofference we can triangulate the position of a given pulse on the trigger! could be a fun thing to check
	//this is after the calibration
	hTimeTOF0.Fill(t11 - t00); //positive tof
	hTimeTOF1.Fill(t13 - t02);
	hTimeTOF2.Fill(t10 - t01);
	hTimeTOF3.Fill(t12 - t03);
    	
	double t0 = (signalTime->at(8).at(indices.at(8)) + signalTime->at(9).at(indices.at(9)) + signalTime->at(10).at(indices.at(10)) + signalTime->at(11).at(indices.at(11)))/4.;
	double t1 = (signalTime->at(12).at(indices.at(12)) + signalTime->at(13).at(indices.at(13)) + signalTime->at(14).at(indices.at(14)) + signalTime->at(15).at(indices.at(15)))/4.;
	double tof = t1-t0;

	double act0a = peakVoltage->at(0).at(indices.at(0)) + peakVoltage->at(1).at(indices.at(1));
	double act1a = peakVoltage->at(2).at(indices.at(2)) + peakVoltage->at(3).at(indices.at(3));
	double act2a = peakVoltage->at(4).at(indices.at(4)) + peakVoltage->at(5).at(indices.at(5));
	double act3a = peakVoltage->at(6).at(indices.at(6)) + peakVoltage->at(7).at(indices.at(7));

	hTOFACT0A.Fill(tof, act0a);
        hTOFACT1A.Fill(tof, act1a);
        hTOFACT2A.Fill(tof, act2a);
        hTOFACT3A.Fill(tof, act3a);

	double act0c = intCharge->at(0).at(indices.at(0)) + intCharge->at(1).at(indices.at(1));
	double act1c = intCharge->at(2).at(indices.at(2)) + intCharge->at(3).at(indices.at(3));
	double act2c = intCharge->at(4).at(indices.at(4)) + intCharge->at(5).at(indices.at(5));
	double act3c = intCharge->at(6).at(indices.at(6)) + intCharge->at(7).at(indices.at(7));

	hTOFACT0C.Fill(tof, act0c);
  	hTOFACT1C.Fill(tof, act1c);
 	hTOFACT2C.Fill(tof, act2c);
  	hTOFACT3C.Fill(tof, act3c);

  	hACT1CACT3C.Fill(act1c, act3c);
  	hACT3CACT2C.Fill(act3c, act2c);
  	hACT2CACT1C.Fill(act2c, act1c); 

	//acraplet
	
		
        hTOFAll.Fill(tof);
	hTOFAllLow.Fill(tof);
	
        hT0.Fill(t0);
        hT1.Fill(t1);
	
  	bool pass = true;
  	bool isEl = false;
	
  	bool passed_act1a_cuts = false;
  	bool passed_act2a_cuts = false;

        switch(momentum)	  {

	case 320: { //; jk
	  // placeholder for momentum dependent cuts
	  break;
	} // 320

	  /*
	  // placeholder
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
	  if (act0a > 1.) {
	    isEl = true;  
	  }
	  if (act1a > 2.) {
	    passed_act1a_cuts = true;
	  }
	  if (act2a > 2.) {
	    passed_act1a_cuts = true; 
	  }
	  if (act3a > 1.) {
	    isEl = true;  
	  }
	} // default
	  
	} // case

        if (!pass) continue;
	
      
        if (isEl) {
	  // electrons
	  hTOFEl.Fill(tof);
	  hTOFElLow.Fill(tof);

	  hTOFACT1A_el.Fill(tof, act1a);
	  hTOFACT2A_el.Fill(tof, act2a);
	  hTOFACT3A_el.Fill(tof, act3a);
	      
	  hTOFACT1C_el.Fill(tof, act1c);
	  hTOFACT2C_el.Fill(tof, act2c);
	  hTOFACT3C_el.Fill(tof, act3c);
        }
        else {

	  // non-electrons
	  hTOFOther.Fill(tof);
	  hTOFOtherLow.Fill(tof);
	
	 
        } // non-electrons

	if (act1a < 1.) { // custom electron removal cut
	
	  hTOFAll.Fill(tof);
	  hTOFAllWide.Fill(tof);
	  
	  // TO VALIDATE the cut values!!!
	  if (passed_act1a_cuts) {
	    hTOFACT0A_act1cut.Fill(tof, act0a);
	    hTOFACT1A_act1cut.Fill(tof, act1a);
	    hTOFACT2A_act1cut.Fill(tof, act2a);
	    hTOFACT3A_act1cut.Fill(tof, act3a);

	    hTOFACT0C_act1cut.Fill(tof, act0c);
	    hTOFACT1C_act1cut.Fill(tof, act1c);
	    hTOFACT2C_act1cut.Fill(tof, act2c);
	    hTOFACT3C_act1cut.Fill(tof, act3c);

	    hTOF_act1cut.Fill(tof);
	  } // act1 cuts

	  // TO VALIDATE the cut values!!!
	  if (passed_act2a_cuts > 2.) { 
	    hTOFACT0A_act2cut.Fill(tof, act0a);
	    hTOFACT1A_act2cut.Fill(tof, act1a);
	    hTOFACT2A_act2cut.Fill(tof, act2a);
	    hTOFACT3A_act2cut.Fill(tof, act3a);

	    hTOFACT0C_act2cut.Fill(tof, act0c);
	    hTOFACT1C_act2cut.Fill(tof, act1c);
	    hTOFACT2C_act2cut.Fill(tof, act2c);
	    hTOFACT3C_act2cut.Fill(tof, act3c);

	    hTOF_act2cut.Fill(tof);
	  } // act2 cuts

	} // custom electron removal cut


	
    } // entries
    cout << "End of event loop!" << endl;
    cout << "Writing histograms to the output file..." << endl;
    outFile.Write();
    cout << "DONE!" << endl;
}


