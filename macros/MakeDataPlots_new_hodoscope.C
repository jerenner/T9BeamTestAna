#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSystem.h"

#include <string>
#include <vector>
#include <iostream>

#include "EventInfo.h"
#include "channelReadClass.h"


// new
using namespace std;

// Matej Pavin 2022
// modified by Jiri Kvita 2022
// modified for 2023 32 channels July 2023
// modified for different structures of the ntuples multiple times

// ______________________________________________________________


// https://docs.google.com/spreadsheets/d/1QBHKEbpC_roTHyY5QJSExFnMmDGn6aLyWtHHhrKORZA/edit?usp=sharing


const int nMaxChannels = 32;
const int nChannels = 31; 
TString treeNames[nMaxChannels] = {
  "ACT0L",    "ACT0R",
  "ACT1L",    "ACT1R",
  "ACT3L", 	"ACT3R",
  "TOF00", 	"TOF01", 	"TOF02", 	"TOF03",
  "TOF10", 	"TOF11", 	"TOF12", 	"TOF13",
  "PbGlass",
  "HDCH8",   "HDCH9",   "HDCH10",   "HDCH11",   "HDCH12", "HDCH13",   "HDCH14",
  "HDCH0",   "HDCH1",   "HDCH2",  "HDCH3",  "HDCH4",  "HDCH5",  "HDCH6",  "HDCH7", 


  //  "", 
};


// ______________________________________________________________
// ______________________________________________________________
// ______________________________________________________________


int getHighestPeakIndex(channelReadClass *reader)
 {
   int imax = -1;
   double maxA = -999;
   double a;
   for (int ipeak = 0; ipeak < reader -> nPeaks; ++ipeak) {
     a = reader -> PeakVoltage[ipeak];
     if (a > maxA) {
       maxA = a;
       imax = ipeak;
     }
   }
   return imax;
   //   return 0; 
 }

// ______________________________________________________________

double GetBeta(double mass, double momentum) {
  double bg = momentum/mass;
  double beta = sqrt(bg*bg/(1+bg*bg));
  return beta;
}


// ______________________________________________________________

void MakeDataPlots_new_hodoscope(string fileName, int momentum) {

  gSystem->Exec("mkdir -p histos/");

  TFile *infile = new TFile(fileName.c_str(), "READ");
  EventInfo *eventInfo = new EventInfo(infile, "EventInfo");
  
  int ent[nChannels];
  channelReadClass *reader[nChannels];
  TTree *trees[nChannels];
  map<TString, channelReadClass*> readerMap;

  for (int ich = 0; ich < nChannels; ++ich) {
    cout << "Initializing " <<  treeNames[ich] << endl;
    reader[ich] = new channelReadClass(infile, treeNames[ich]);
    trees[ich] = reader[ich] -> fChain;
    ent[ich] = trees[ich] -> GetEntries();
    readerMap[treeNames[ich]] = reader[ich];
  }

  double tofmin = 10.;
  double tofmax = 40.;
  int ntofbins = 200;

  double tofminlow = 10.;
  double tofmaxlow = 20.;
  int ntofbinslow = 100;

  int ntofbins2d = 400;

  double actChargeMax = 1.0;
  double actAmplitudeMax =  20.;

  TString peakModeTag = "";
  //if (peakMode != "")
  //  peakModeTag = "_" + peakMode;
  TString outFileName = TString(fileName.substr(0, fileName.size()-5).c_str()) + "_plots" + peakModeTag + ".root";
  outFileName = outFileName.ReplaceAll("output/", "histos/").ReplaceAll("ntuple_files/","histos/");

  TFile outFile(outFileName.Data(), "RECREATE");
  outFile.cd();

  // TOF 1D
  TH1D hTOFAll("hTOFAll", ";t_{TOF}^{All} [ns]", 120, tofmin, tofmax);
  TH1D hTOFAllWide("hTOFAllWide", ";t_{TOF}^{All} [ns]", 2*ntofbins, tofmin, 2*tofmax);
  TH1D hTOFEl("hTOFEl", ";t_{TOF}^{e} [ns]", ntofbins, tofmin, tofmax);
  TH1D hTOFOther("hTOFOther", ";t_{TOF}^{non-e} [ns]", ntofbins, tofmin, tofmax);

  TH1D hTOFAllLow("hTOFAllLow", ";t_{TOF}^{All} [ns]", ntofbinslow, tofminlow, tofmaxlow);
  TH1D hTOFElLow("hTOFElLow", ";t_{TOF}^{e} [ns]", ntofbinslow, tofminlow, tofmaxlow);
  TH1D hTOFOtherLow("hTOFOtherLow", ";t_{TOF}^{non-e} [ns]", ntofbinslow, tofminlow, tofmaxlow);

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
  TH1D hTimeTOF0("hTimeTOF0", "; hTimeTOF0", 100, 0.,50.);
  TH1D hTimeTOF1("hTimeTOF1", "; hTimeTOF1", 100, 0.,50.);
  TH1D hTimeTOF2("hTimeTOF2", "; hTimeTOF2", 100, 0.,50.);
  TH1D hTimeTOF3("hTimeTOF3", "; hTimeTOF3", 100, 0.,50.);

  TH2D hPbACT0A("hRef_pbA_act0A", "; Pb-glass Amplitude ; ACT0 Amplitude", 200, 0., actAmplitudeMax, 400, 0., actAmplitudeMax);
  TH2D hPbACT0C("hRef_pbC_act0C", "; Pb-glass Charge ; ACT1 Charge)", 200, 0., actChargeMax, 400, 0., actAmplitudeMax);
  TH2D hPbACT1A("hRef_pbA_act1A", "; Pb-glass Amplitude ; ACT1 Amplitude", 200, 0., actAmplitudeMax, 400, 0., actAmplitudeMax);
  TH2D hPbACT1C("hRef_pbC_act1C", "; Pb-glass Charge ; ACT1 Charge)", 200, 0., actChargeMax, 400, 0., actAmplitudeMax);



  // also ACT 0 and 1, separately:
  /* seems they were already defined below...
  TH2D hTOFACT0A("hRef_TOFACT0A", "; t_{1}-t_{0} [ns]; ACT0 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
  TH2D hTOFACT1A("hRef_TOFACT1A", "; t_{1}-t_{0} [ns]; ACT1 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
  TH2D hTOFACT0C("hRef_TOFACT0C", "; t_{1}-t_{0} [ns]; ACT0 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);
  TH2D hTOFACT1C("hRef_TOFACT1C", "; t_{1}-t_{0} [ns]; ACT1 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);
  */

  //TOF vs Pb-glass plots
  TH2D hPbATOF("hRef_PbATOF", "; Pb-glass Amplitude; t_{1}-t_{0} [ns]", 200, 0., actAmplitudeMax/2, ntofbins2d, tofmin, tofmax);
  TH2D hPbCTOF("hRef_PbCTOF", "; Pb-glass Charge; t_{1}-t_{0} [ns]", 200, 0., actChargeMax, ntofbins2d, tofmin, tofmax);
  TH2D hTOFPbA("hRef_TOFPbA", "; t_{1}-t_{0} [ns]; Pb-glass Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax/2);
  
  //acraplet - investigate "weird electrons"
  TH2D hHC0AHC1A("hweirdE_HC0AHC1A", "; Hole Counter 0 Amplitude; Hole Counter 1 Amplitude", 200, 0., 1000, 200, 0., 1000.);
  TH2D hHC0CHC1C("hweirdE_HC0CHC1C", "; Hole Counter 0 Charge; Hole Counter 1 Charge", 200, 0., 1., 200, 0., 1.);

  // standard
  vector<TH1D> hCharge;
  vector<TH1D> hVoltage;
  vector<TH1D> hPedestalSigma;
  vector<TH1D> hTime;
  vector<TH1D> hnPeaks;

  // no cuts
  TH2D hTOFACT0A("hRef_TOFACT0A", "; t_{1}-t_{0} [ns]; ACT0 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
  TH2D hTOFACT1A("hRef_TOFACT1A", "; t_{1}-t_{0} [ns]; ACT1 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
  TH2D hTOFACT3A("hRef_TOFACT3A", "; t_{1}-t_{0} [ns]; ACT3 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);

  TH2D hTOFACT0C("hRef_TOFACT0C", "; t_{1}-t_{0} [ns]; ACT0 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);
  TH2D hTOFACT1C("hRef_TOFACT1C", "; t_{1}-t_{0} [ns]; ACT1 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);
  TH2D hTOFACT3C("hRef_TOFACT3C", "; t_{1}-t_{0} [ns]; ACT3 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);


  // 2D ACT charges
  TH2D hACT1CACT3C("hRef_ACT1CACT3C", "; ACT1 Charge; ACT3 Charge", 200, 0., actChargeMax, 200, 0., actChargeMax);

  // nPeak 2D plots;)
  int nbn = 16.;
  double n1 = 0.;
  double n2 = 4.;

  TH2D hnPeaksToF1vsnPeaksToF0("hnPeaksToF1vsnPeaksToF0", "hnPeaksToF1vsnPeaksToF0;<n_{Peaks}^{ToF0}>;<n_{Peaks}^{ToF1}>", nbn, n1, n2, nbn, n1, n2);
  
  TH2D hnPeaksToFvsToF("hnPeaksToFvsToF", "hnPeaksToFvsToF;t_{TOF};<n_{Peaks}^{ToF}>", ntofbins2d/4, tofmin, tofmax, nbn, n1, n2);
  TH2D hnPeaksToFvsToFlow("hnPeaksToFvsToFlow", "hnPeaksToFvsToF;t_{TOF};<n_{Peaks}^{ToF}>", ntofbins2d/4, tofminlow, tofmaxlow, nbn, n1, n2);
  TH2D hnPeaksToFvsLeadGlassA("hnPeaksToFvsLeadGlassA", "hnPeaksToFvsLeadGlassA;lead glass A;<n_{Peaks}^{ToF}>", 100,  0., actAmplitudeMax/2., nbn, n1, n2);
  n1 = 0.;
  n2 = 10.;
  TH2D hnPeaksLeadGlassvsLeadGlassA("hnPeaksLeadGlassvsLeadGlassA", "hnPeaksLeadGlassvsLeadGlassA;lead glass A;n_{Peaks}^{Pb}", 100,  0., actAmplitudeMax/2., int(n2-n1), n1, n2);
  
  for(int i = 0; i < nChannels; i++) {
    string name1 = "hRef_Charge" + to_string(i);
    string name2 = "hRef_Voltage" + to_string(i);
    string name3 = "hRef_Hits" + to_string(i);
    string name4 = "hRef_PedestalSigma" + to_string(i);
    string name5 = "hRef_Time" + to_string(i);
    string name6 = "hRef_nPeaks" + to_string(i);

    string title1 = "Channel " + to_string(i) + "; Charge [nC]; Triggers";
    string title2 = "Channel " + to_string(i) + "; Total Amplitude [V]; Triggers";
    string title3 = "Channel " + to_string(i) + "; Hits per trigger; Triggers";
    string title4 = "Channel " + to_string(i) + "; #sigma_{ped} [V]; Triggers";
    string title5 = "Channel " + to_string(i) + "; Time [ns]; Triggers";
    string title6 = "Channel " + to_string(i) + "; Number of peaks; Triggers";

    TH1D temp1(name1.c_str(), title1.c_str(), 240, -0.16, 30*0.08);
    TH1D temp2(name2.c_str(), title2.c_str(), 320, 0., 13*0.8);
    TH1D temp3(name3.c_str(), title3.c_str(), 5, -0.5, 4.5);
    TH1D temp4(name4.c_str(), title4.c_str(), 200, 0., 0.01);
    TH1D temp5(name5.c_str(), title5.c_str(), 270, 0., 540.);
    TH1D temp6(name6.c_str(), title6.c_str(), 20, 0., 20.);

    hCharge.push_back(temp1);
    hVoltage.push_back(temp2);
    hPedestalSigma.push_back(temp4);
    hTime.push_back(temp5);
    hnPeaks.push_back(temp6);
  }


  
  // hodoscope runs
  // lead glass A vs hodoscope occupancy with some amplitude cuts
  // subject to mV callibration!!
  // jiri on shift 28.7.2023
  TString name = "LeadGlassPhotonAVsPositronHodoOcc";
  TString title = name + ";HD Channel ID;A^{#gamma}_{Pb}";
  TH2D LeadGlassPhotonAVsPositronHodoOcc(name, title, 15, 0, 15, 125, 0, 5);

  name = "LeadGlassPhotonAVsPositronMaxHodoOcc";
  title = name + ";HD Max. Channel ID;A^{#gamma}_{Pb}";
  TH2D LeadGlassPhotonAVsPositronMaxHodoOcc(name, title, 15, 0, 15, 125, 0, 5);

  name = "HodoOccScatter";
  title = name + ";HD channel;HD channel;entries";
  TH2D HodoOccScatter(name, title, 15, 0, 15, 15, 0, 15);
   
  name = "HodoOccScatterFrac";
  title = name + ";HD channel;HD channel;fractions";
  TH2D HodoOccScatterFrac(name, title, 15, 0, 15, 15, 0, 15);
   
  // +-------------------------------+
  // |         event loop            |
  // +-------------------------------+

  cout << "Event loop!" << endl;
  int verbose = 10000;
  int debug = 0;
  // TODO:
  // check also the number of entries in the trees?
  
  for(int ientry = 0; ientry < ent[0]; ientry++) {

    if (ientry % verbose == 0) {
      cout << "processing " << ientry << " / " << ent[0] << endl;
    }

    eventInfo -> GetEntry(ientry);
    Long64_t  RunNumber = eventInfo->RunNumber;
    Int_t EventNumber = eventInfo -> EventNumber;
    Int_t SpillNumber = eventInfo -> SpillNumber;
    //    cout << " RunNumber=" << RunNumber << " EventNumber=" << EventNumber << " SpillNumber=" << SpillNumber << endl;
    
    for (int ich = 0; ich < nChannels; ++ich) {
      if (debug)	cout << "getting entry for " <<  treeNames[ich] << endl;
      reader[ich] -> GetEntry(ientry);
    }
    if (debug)      cout << "done" << endl;


  
    //    vector<int> indices(nChannels, 0);

    bool onePeakInAllACTs = true;
    bool onePeakInAllToFs = true;
    bool onePeakInAll = true;

    bool moreThanOnePeakInAllACTs = true;
    bool moreThanOnePeakInAllToFs = true;
    bool moreThanOnePeakInAll = true;

    bool PbGlassAboveElectronLevel = true;

    double PbGlassElectronThreshA = 5;
    double PbGlassElectronUpperThreshA = 6.5;

    double ACTC23ElectronThreshA = 1.5;
    double ACTC23ElectronUpperThreshA = 3.5;

    if (debug)      cout << "point a" << endl;
    bool onePeakInPbGlass = (reader[18] -> nPeaks == 1);
    if (debug)      cout << "point b" << endl;



    
    // TOF trigger scintilators

    map<TString,int> PeakID;
    map<TString,double> Amplitudes; // amplitude
    map<TString,double> Charges; // charge
    map<TString,double> PeakTimes; // time

    // read all channels information for all waveforms!

    for (int ich = 0; ich < nChannels; ++ich) {
      TString chname = treeNames[ich];
      if (debug)      cout << "point c, " << chname.Data() << endl;

      PeakID[chname] = getHighestPeakIndex(readerMap[chname]);
      int ipeak = PeakID[chname];
      if ( ipeak >= 0 && ipeak < readerMap[chname] -> nPeaks) {
	Amplitudes[chname] = readerMap[chname] -> PeakVoltage[ipeak];
	Charges[chname] = readerMap[chname] -> IntCharge[ipeak];
	PeakTimes[chname] = readerMap[chname] -> SignalTime[ipeak];
	
	// histograms over all channels
	// can be simplified using the above maps
    	hCharge.at(ich).Fill(reader[ich] -> IntCharge[ipeak]);
	hVoltage.at(ich).Fill(reader[ich] -> PeakVoltage[ipeak]);
	hTime.at(ich).Fill(reader[ich] -> SignalTime[ipeak]);
	//hNbPeaks.at(ich).Fill(reader[ich] -> nPeaks);
	hPedestalSigma.at(ich).Fill(reader[ich] -> PedestalSigma);
	hnPeaks.at(ich).Fill(reader[ich] -> nPeaks);
	
      } else {
	Amplitudes[chname] = 0.;
	Charges[chname] = 0.;
	PeakTimes[chname] = 0.;
      }
      if (debug)      cout << "point d" << endl;

    }


    // jiri, hodoscope
    double pbanew = reader[16] -> PeakVoltage[0];
    int maxhdchid = -1;
    double maxA = -1;
    bool hits[15];
    for(int j = 0; j < 15; ++j)
      hits[j] = false;
    for(int j = 0; j < nChannels; j++) {
      if (j < 17)
	continue;
      int hdchid = -1;
      double Aj = reader[j] -> PeakVoltage[0];
      if (j >= 17 && j < 24) {
	hdchid = j - 9;
      } else {
	hdchid = j - 24;
      }
      if (Aj > 1.5) {
	hits[hdchid] = true;
	LeadGlassPhotonAVsPositronHodoOcc.Fill(hdchid, pbanew);
	if (Aj > maxA) {
	  maxA = Aj;
	  maxhdchid = hdchid;
	}
      }  // A cut
    } // channels
    
    if (maxhdchid >= 0)
      LeadGlassPhotonAVsPositronMaxHodoOcc.Fill(maxhdchid, pbanew);
    for(int j = 0; j < 15; ++j) {
      for(int k = 0; k <= j; ++k) {
	if (hits[j] && hits[k])
	  HodoOccScatter.Fill(j,k); 
      }
    }
    
    
    
    
    // TOF trigger scintilators

    double t00 = PeakTimes["TOF00"];
    double t01 = PeakTimes["TOF01"];
    double t02 = PeakTimes["TOF02"];
    double t03 = PeakTimes["TOF03"];

    double t10 = PeakTimes["TOF10"];
    double t11 = PeakTimes["TOF11"];
    double t12 = PeakTimes["TOF12"];
    double t13 = PeakTimes["TOF13"];
    
    // JK's time resolution of 2022
    // diagonal combinations
    double t0a = (t00 + t03) / 2.;
    double t0b = (t01 + t02) / 2.;
    double t1a = (t10 + t13) / 2.;
    double t1b = (t11 + t12) / 2.;

    // time diffs for time resolution histogramme 2022:
    double t0diff = t0a - t0b;
    double t1diff = t1a - t1b;

    // jiri
    // Fill resolution histograms of 2022:
    hTimeReso0.Fill(t0diff);
    hTimeReso1.Fill(t1diff);
    hTimeReso0_zoom.Fill(t0diff);
    hTimeReso1_zoom.Fill(t1diff);

    hTimeDiffTOF01.Fill(t00 - t01); // carefull, this particular histogram is a difference of hit times,
    hTimeDiffTOF02.Fill(t00 - t02); // not a TOF per se but instead the cable length difference
    hTimeDiffTOF03.Fill(t00 - t03); // plus the photon travel time though the panel
 
    hTimeDiffTOF11.Fill(t11 - t10);
    hTimeDiffTOF12.Fill(t11 - t12);
    hTimeDiffTOF13.Fill(t11 - t13);

    // acraplet
    // compare the hit times for the same event recorded by trigger PMTs on the same side (up/down, left/right)
    // assumption: the light travel time trough the trigger counter to the PMT should be about the same (if the beam is well aligned)
    // then we can check if the TOF is constant for a run
    // idea: using tof dofference we can triangulate the position of a given pulse on the trigger! could be a fun thing to check
    // this is after the calibration
    hTimeTOF0.Fill(t11 - t00); // positive tof
    hTimeTOF1.Fill(t13 - t02);
    hTimeTOF2.Fill(t10 - t01);
    hTimeTOF3.Fill(t12 - t03);

    // Time of flight!
    double t0 = ( t00 + t01 + t02 + t03 ) / 4.;
    double t1 = ( t10 + t11 + t12 + t13 ) / 4.;
    double tof = t1 - t0;

    // ACTs
    double act0c = Charges["ACT0L"] + Charges["ACT0R"];
    double act1c = Charges["ACT1L"] + Charges["ACT1R"];
    double act3c = Charges["ACT3L"] + Charges["ACT3R"];

    double act0a = Amplitudes["ACT0L"] + Amplitudes["ACT0R"];
    double act1a = Amplitudes["ACT1L"] + Amplitudes["ACT1R"];
    double act3a = Amplitudes["ACT3L"] + Amplitudes["ACT3R"];

    // hole counters and lead glass
    
    double hc0c = Charges["Hole0"];
    double hc0a = Amplitudes["Hole0"];

    double hc1c = Charges["Hole1"];
    double hc1a = Amplitudes["Hole1"];

    double pbc = Charges["PbGlass"];
    double pba = Amplitudes["PbGlass"];
    // HACK!
    //    pba = readerMap["PbGlass"] -> PeakVoltage[0];


    
    // acraplet, hodoscope setup
    std::vector<int> channelToHodoscope = {8, 9, 10, 11, 12, 13, 14, 0, 1, 2, 3, 4, 5, 6, 7};
    double threshHodoscopeHit = 1.5; //Threshold estimated by hand using the online hist as reference 
    //each channel has a different 1pe peak (in particular ch11 has a low one) 
    // the online analysis threshold of 400mV corresponds to 1.9 in these units but we are cutting a lot of hits especially in channel 11, using 1.5 is better there  
    //careful ! position of the detectors on the digitiser have moved!!!
    for (int i=0; i<15; i++){
      if (peakVoltage[17+i][0] >= threshHodoscopeHit){
	hnHitsHodoscope.Fill(channelToHodoscope[i]);
      }
    }
    
 
    // amplitudes vs tof
    hTOFACT0A.Fill(tof, act0a);
    hTOFACT1A.Fill(tof, act1a);
    hTOFACT3A.Fill(tof, act3a);

    // lead glass vs acts and tof
    
    hPbACT0A.Fill(pba, act0a);
    hPbACT1A.Fill(pba, act1a);
    
    hPbATOF.Fill(pba, tof);
    hTOFPbA.Fill(tof, pba);

    // charges vs tof
   
    hTOFACT0C.Fill(tof, act0c);
    hTOFACT1C.Fill(tof, act1c);
    hTOFACT3C.Fill(tof, act3c);

    // lead glass vs acts and tof
    hPbACT0C.Fill(pbc, act0c);
    hTOFACT0C.Fill(tof, act0c);
    hPbACT1C.Fill(pbc, act1c);
    hTOFACT1C.Fill(tof, act1c);
    
    hPbCTOF.Fill(pbc, tof);

    // act 2d plots
    hACT1CACT3C.Fill(act1c, act3c);

    // acraplet - weird electrons which do not see anything in the ACT
    if (act3 >= 1.5 && tof >= 13.5 && tof <= 16.5) {
      //      if (act23aAver != 1.5 && tof >= 13.5 && tof <= 16.5) {
      hHC0AHC1A.Fill(hc0a, hc1a);
      hHC0CHC1C.Fill(hc0c, hc1c);
    }

    // acraplet
    hTOFAll.Fill(tof);
    hTOFAllLow.Fill(tof);

    hT0.Fill(t0);
    hT1.Fill(t1);


    // 2D <nPeaks> studies
    int nPeaksToFAver = (readerMap["TOF00"]->nPeaks + readerMap["TOF01"]->nPeaks + readerMap["TOF02"]->nPeaks + readerMap["TOF03"]->nPeaks + readerMap["TOF10"]->nPeaks + readerMap["TOF11"]->nPeaks + readerMap["TOF12"]->nPeaks + readerMap["TOF13"]->nPeaks) / 8;

    int nPeaksToF0Aver = (readerMap["TOF00"]->nPeaks + readerMap["TOF01"]->nPeaks + readerMap["TOF02"]->nPeaks + readerMap["TOF03"]->nPeaks ) / 4;
    int nPeaksToF1Aver = (readerMap["TOF10"]->nPeaks + readerMap["TOF11"]->nPeaks + readerMap["TOF12"]->nPeaks + readerMap["TOF13"]->nPeaks) / 4;
    int nPeaksACT3Aver = (readerMap["ACT3L"]->nPeaks + readerMap["ACT3R"]->nPeaks) / 2;

    int nPeaksLeadGlass = readerMap["PbGlass"]->nPeaks;
   

    hnPeaksToF1vsnPeaksToF0.Fill(nPeaksToF1Aver, nPeaksToF0Aver);

    hnPeaksToFvsToF.Fill(tof, nPeaksToFAver);
    hnPeaksToFvsToFlow.Fill(tof, nPeaksToFAver);
    hnPeaksLeadGlassvsLeadGlassA.Fill(pba, nPeaksLeadGlass);

    hnPeaksToFvsLeadGlassA.Fill(pba, nPeaksToFAver);

    // selections:

    bool pass = true;
    bool isEl = false;

    bool passed_act1a_cuts = false;

    switch(momentum)	  {

      case 320: { //; jk
        // placeholder for momentum dependent cuts
        break;
      } // 320

      default: {
        if (ientry < 10)
          cout << "WARNING: Using default settings for the " << momentum << " MeV/c beam" << endl;
	// add ACT1 cuts?
	//        if ( (act2a + act3a)/2. > 3.) { // custom electron removal cut
	// isEl = true;
	// }
      } // default
    } // case

    if (!pass) continue;

    if (isEl) {
      // electrons
      hTOFEl.Fill(tof);
      hTOFElLow.Fill(tof);
    }
    else {
      // non-electrons
      hTOFOther.Fill(tof);
      hTOFOtherLow.Fill(tof);
    } // non-electrons

    hTOFAll.Fill(tof);
    hTOFAllWide.Fill(tof);

  } // entries
  cout << "End of event loop!" << endl;
  outFile.Write();


}





