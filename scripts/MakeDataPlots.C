#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSystem.h"

#include <string>
#include <vector>
#include <iostream>
// new
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
// peakMode: "", a, b, c, d, e, f

void MakeDataPlots(string fileName, int momentum, double x0_cut = 4.9, double y0_cut = 6, double act23_pi_minA =0.4, double act23_pi_maxA = 2.4, double pb_min = 0.1, double actThresh = .5, TString peakMode = "") {

  const int nMaxChannels = 32;
  const int nChannels = 19; // UPDATE THIS FOR HODOSCOPE PMTs! to e.g. 32!

  Double_t peakVoltage[nMaxChannels][1];
  Double_t peakTime[nMaxChannels][1];
  Double_t signalTime[nMaxChannels][1];
  Double_t intCharge[nMaxChannels][1];
  Double_t pedestal[nMaxChannels];
  Double_t pedestalSigma[nMaxChannels];
  Double_t nPeaks[nMaxChannels];

  gSystem->Exec("mkdir -p histos/");

  TFile inFile(fileName.c_str(), "READ");
  TTree *tree = (TTree*) inFile.Get("anaTree");

  tree->SetBranchAddress("PeakVoltage",&peakVoltage);
  tree->SetBranchAddress("PeakTime",&peakTime);
  tree->SetBranchAddress("SignalTime",&signalTime);
  tree->SetBranchAddress("IntCharge",&intCharge);
  tree->SetBranchAddress("Pedestal",&pedestal);
  tree->SetBranchAddress("PedestalSigma",&pedestalSigma);
  //tree->SetBranchAddress("NbPeaks",&nbPeaks);
  //tree->SetBranchAddress("PassThreshold",&passThreshold);
  tree->SetBranchAddress("nPeaks",&nPeaks);

  int ent = tree->GetEntries();

  double tofmin = 10.;
  double tofmax = 40.;
  int ntofbins = 400;

  double tofminlow = 10.;
  double tofmaxlow = 20.;
  int ntofbinslow = 100;

  int ntofbins2d = 400;

  double actChargeMin = -0.04;
  double actChargeMax = 2.0;
  double actAmplitudeMax =  20;

  TString peakModeTag = "";
  if (peakMode != "")
    peakModeTag = "_" + peakMode;
  TString outFileName = TString(fileName.substr(0, fileName.size()-5).c_str()) + "_plots" + peakModeTag + ".root";
  outFileName = outFileName.ReplaceAll("output/", "histos/");
  outFileName = outFileName.ReplaceAll("data/", "histos/");

  TFile outFile(outFileName.Data(), "RECREATE");
  outFile.cd();

  TH1D hTOFAll("hTOFAll", ";t_{TOF}^{All} [ns]", 120, tofmin, tofmax);
  TH1D hTOFAllWide("hTOFAllWide", ";t_{TOF}^{All} [ns]", 2*ntofbins, tofmin, 2*tofmax);
  TH1D hTOFEl("hTOFEl", ";t_{TOF}^{e} [ns]", ntofbins, tofmin, tofmax);
  TH1D hTOFOther("hTOFOther", ";t_{TOF}^{non-e} [ns]", ntofbins, tofmin, tofmax);

  double tofmax_particle = 40;

  TH1D hTOFMuACT23pb("hTOFMuACT23pb", ";t_{TOF} [ns]", ntofbins, tofmin, tofmax_particle);
  TH1D hTOFElACT23pb("hTOFElACT23pb", ";t_{TOF} [ns]", ntofbins, tofmin, tofmax_particle);
  TH1D hTOFPiACT23pb("hTOFPiACT23pb", ";t_{TOF} [ns]", ntofbins, tofmin, tofmax_particle);
  TH1D hTOFpACT23pb("hTOFpACT23pb", ";t_{TOF} [ns]", ntofbins, tofmin, tofmax_particle);
  TH1D hTOFdACT23pb("hTOFdACT23pb", ";t_{TOF} [ns]", ntofbins, tofmin, tofmax_particle);



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

  //lead glass vs act 2 and 3 - identify particles
  TH2D hPbACT23A("hRef_pbA_act23A", "; Pb-glass Amplitude ; (ACT2+ACT3)/2 Amplitude", 400, 0., actAmplitudeMax, 400, 0., actAmplitudeMax);
  TH2D hPbACT1A("hRef_pbA_act1A", "; Pb-glass Amplitude ; ACT1 mean Amplitude", 400, 0., actAmplitudeMax, 400, 0., actAmplitudeMax);
  TH2D hPbACT23C("hRef_pbC_act23C", "; Pb-glass Charge ; (ACT2+ACT3)/2 Charge)", 400, 0., actChargeMax, 400, 0., actAmplitudeMax);
  TH2D hPbCCT1C("hRef_pbA_act1C", "; Pb-glass Charge (nC) ; ACT1 mean Charge (nC)", 400,  actChargeMin, actChargeMax, 400, actChargeMin, actChargeMax);


  //ACT2+ACT3/2 vs TOF plots
  TH2D hTOFACT23A("hRef_TOFACT23A", "; t_{1}-t_{0} [ns]; (ACT2+ACT3)/2 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
  TH2D hTOFACT23C("hRef_TOFACT23C", "; t_{1}-t_{0} [ns]; (ACT2+ACT3)/2 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);


  //TOF vs Pb-glass plots
  TH2D hPbATOF("hRef_PbATOF", "; Pb-glass Amplitude; t_{1}-t_{0} [ns]", 400, 0., actAmplitudeMax/2, ntofbins2d, tofmin, tofmax);
  TH2D hPbCTOF("hRef_PbCTOF", "; Pb-glass Charge; t_{1}-t_{0} [ns]", 400, 0., actChargeMax, ntofbins2d, tofmin, tofmax);

  //acraplet - investigate "weird electrons"
  TH2D hHC0AHC1A("hweirdE_HC0AHC1A", "; Hole Counter 0 Amplitude; Hole Counter 1 Amplitude", 200, 0., 1000, 200, 0., 1000.);

  TH2D hHC0CHC1C("hweirdE_HC0CHC1C", "; Hole Counter 0 Charge; Hole Counter 1 Charge", 200, 0., 1., 200, 0., 1.);

  // standard
  vector<TH1D> hCharge;
  vector<TH1D> hVoltage;
  vector<TH1D> hHit;
  vector<TH1D> hPedestalSigma;
  vector<TH1D> hTime;
  vector<TH1D> hnPeaks;
  vector<TH1D> hPedestal;
  vector<TH2D> hPedestalNbPeaks;

  // no cuts
  TH2D hTOFACT0A("hRef_TOFACT0A", "; t_{1}-t_{0} [ns]; ACT0 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
  TH2D hTOFACT1A("hRef_TOFACT1A", "; t_{1}-t_{0} [ns]; ACT1 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
  TH2D hTOFACT2A("hRef_TOFACT2A", "; t_{1}-t_{0} [ns]; ACT2 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
  TH2D hTOFACT3A("hRef_TOFACT3A", "; t_{1}-t_{0} [ns]; ACT3 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);

  TH2D hTOFACT0C("hRef_TOFACT0C", "; t_{1}-t_{0} [ns]; ACT0 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);
  TH2D hTOFACT1C("hRef_TOFACT1C", "; t_{1}-t_{0} [ns]; ACT1 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);
  TH2D hTOFACT2C("hRef_TOFACT2C", "; t_{1}-t_{0} [ns]; ACT2 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);
  TH2D hTOFACT3C("hRef_TOFACT3C", "; t_{1}-t_{0} [ns]; ACT3 Charge", ntofbins2d, tofmin, tofmax, 200, 0., actChargeMax);


  // ACT2+ACT3 cut
  TH1D hTOF_act2act3cut("hTOF_act2act3cut", "; t_{1}-t_{0} [ns];", 120, tofmin, tofmax);

  // 2D ACT charges
  TH2D hACT2CACT1C("hRef_ACT2CACT1C", "; ACT2 Charge; ACT1 Charge", 200, 0., actChargeMax, 200, 0., actChargeMax);
  TH2D hACT3CACT2C("hRef_ACT3CACT2C", "; ACT3 Charge; ACT2 Charge", 200, 0., actChargeMax, 200, 0., actChargeMax);
  TH2D hACT1CACT3C("hRef_ACT1CACT3C", "; ACT1 Charge; ACT3 Charge", 200, 0., actChargeMax, 200, 0., actChargeMax);

  // nPeak 2D plots;)
  int nbn = 16.;
  double n1 = 0.;
  double n2 = 4.;

  TH2D hnPeaksACT23vsnPeaksToF("hnPeaksACT23vsnPeaksToF", "hnPeaksACT23vsnPeaksToF;<n_{Peaks}^{ToF}>;<n_{Peaks}^{ACT23}>", nbn, n1, n2, nbn, n1, n2);
  TH2D hnPeaksToF1vsnPeaksToF0("hnPeaksToF1vsnPeaksToF0", "hnPeaksToF1vsnPeaksToF0;<n_{Peaks}^{ToF0}>;<n_{Peaks}^{ToF1}>", nbn, n1, n2, nbn, n1, n2);
  TH2D hnPeaksACT3vsnPeaksACT2("hnPeaksACT3vsnPeaksACT2", "hnPeaksACT3vsnPeaksACT2;<n_{Peaks}^{ACT2}>;<n_{Peaks}^{ACT3}>", nbn/2, n1, n2, nbn/2, n1, n2);
  
  TH2D hnPeaksACT23vsToF("hnPeaksACT23vsToF", "hnPeaksACT23vsToF;t_{TOF};<n_{Peaks}^{ACT23}>", ntofbins2d/4, tofmin, tofmax, nbn, n1, n2);
  TH2D hnPeaksACT23vsToFlow("hnPeaksACT23vsToFlow", "hnPeaksACT23vsToF;t_{TOF};<n_{Peaks}^{ACT23}>", ntofbins2d/4, tofminlow, tofmaxlow, nbn, n1, n2);
  TH2D hnPeaksToFvsToF("hnPeaksToFvsToF", "hnPeaksToFvsToF;t_{TOF};<n_{Peaks}^{ToF}>", ntofbins2d/4, tofmin, tofmax, nbn, n1, n2);
  TH2D hnPeaksToFvsToFlow("hnPeaksToFvsToFlow", "hnPeaksToFvsToF;t_{TOF};<n_{Peaks}^{ToF}>", ntofbins2d/4, tofminlow, tofmaxlow, nbn, n1, n2);
  TH2D hnPeaksACT23vsLeadGlassA("hnPeaksACT23vsLeadGlassA", "hnPeaksACT23vsLeadGlassA;lead glass A;<n_{Peaks}^{ACT23}>", 100,  0., actAmplitudeMax/2., nbn, n1, n2);
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
    string name7 = "hRef_Pedestal" + to_string(i);
    string name8 = "hRef_PedestalNbPeaks" + to_string(i);

    string title1 = "Channel " + to_string(i) + "; Charge [nC]; Triggers";
    string title2 = "Channel " + to_string(i) + "; Total Amplitude [V]; Triggers";
    string title3 = "Channel " + to_string(i) + "; Hits per trigger; Triggers";
    string title4 = "Channel " + to_string(i) + "; #sigma_{ped} [V]; Triggers";
    string title5 = "Channel " + to_string(i) + "; Time [ns]; Triggers";
    string title6 = "Channel " + to_string(i) + "; Number of peaks; Triggers";
    string title7 = "Channel " + to_string(i) + "; Pedestal Amplitude; Triggers/1mV";
    string title8 = "Channel " + to_string(i) + "; Pedestal Amplitude; Number of peaks";

    TH1D temp1(name1.c_str(), title1.c_str(), 400, actChargeMin, actChargeMax);
    TH1D temp2(name2.c_str(), title2.c_str(), 400, 0., 2);
    TH1D temp3(name3.c_str(), title3.c_str(), 5, -0.5, 4.5);
    TH1D temp4(name4.c_str(), title4.c_str(), 200, 0., 0.01);
    TH1D temp5(name5.c_str(), title5.c_str(), 270, 0., 540.);
    TH1D temp6(name6.c_str(), title6.c_str(), 20, 0., 20.);
    TH1D temp7(name7.c_str(), title7.c_str(), 1000., 1.65, 1.65+1000*0.0012207);
    TH2D temp8(name8.c_str(), title8.c_str(), 200, 0., 3*0.8,  5, 0., 5.);

    hCharge.push_back(temp1);
    hVoltage.push_back(temp2);
    //    hHit.push_back(temp3);
    hPedestalSigma.push_back(temp4);
    hTime.push_back(temp5);
    hnPeaks.push_back(temp6);
    hPedestal.push_back(temp7);
    hPedestalNbPeaks.push_back(temp8);
  }

  // +-------------------------------+
  // |         event loop            |
  // +-------------------------------+

  cout << "Event loop!" << endl;

  for(int i = 0; i < ent; i++) {

    tree->GetEntry(i);
    vector<int> indices(nChannels, 0);;

    bool onePeakInAllACTs = true;
    bool onePeakInAllToFs = true;
    bool onePeakInAll = true;

    bool moreThanOnePeakInAllACTs = true;
    bool moreThanOnePeakInAllToFs = true;
    bool moreThanOnePeakInAll = true;

    bool PbGlassAboveElectronLevel = true;
    bool ACT23AboveElectronLevel = true;
    bool  noPeaksInACT1 = true;

    bool ACT1AboveElectronLevel = true;

    double PbGlassElectronThreshA = 5; //these values will have to be changed for the new
    //amplitude scale
    double PbGlassElectronUpperThreshA = 6.5;


    double ACTC23ElectronThreshA = 1.5;
    double ACTC23ElectronUpperThreshA = 3.5;

    double ACT1ElectronThreshA = 0.015;
    double ACT1ElectronThreshC = 0.01;

    for(int j = 0; j < nChannels; j++) {
      if (j < 16) {
        onePeakInAll = onePeakInAll && (nPeaks[j] == 1);
        moreThanOnePeakInAll = moreThanOnePeakInAll && (nPeaks[j] > 1);
      }
      if (j < 8) {
        onePeakInAllACTs = onePeakInAllACTs && (nPeaks[j] == 1);
        moreThanOnePeakInAllACTs = moreThanOnePeakInAllACTs && (nPeaks[j] > 1);
      }
      if (j >= 8 && j < 16) {
        onePeakInAllToFs = onePeakInAllToFs && (nPeaks[j] == 1);
        moreThanOnePeakInAllToFs = moreThanOnePeakInAllToFs && (nPeaks[j] > 1);
      }
      //dirty add-on to select electrons
      if (j==18){
        PbGlassAboveElectronLevel = PbGlassAboveElectronLevel && (peakVoltage[j][0] > PbGlassElectronThreshA);
        PbGlassAboveElectronLevel = PbGlassAboveElectronLevel && (peakVoltage[j][0] < PbGlassElectronUpperThreshA);
      }
      if (j==4){
        ACT23AboveElectronLevel = ACT23AboveElectronLevel && ((peakVoltage[j][0] + peakVoltage[j+1][0] + peakVoltage[j+2][0] + peakVoltage[j+3][0])/2. > ACTC23ElectronThreshA);

        ACT23AboveElectronLevel = ACT23AboveElectronLevel && ((peakVoltage[j][0] + peakVoltage[j+1][0] + peakVoltage[j+2][0] + peakVoltage[j+3][0])/2. < ACTC23ElectronUpperThreshA);
      }
      //dirty check to use ACT1 to 'veto' electrons
      if (j==2){
        ACT1AboveElectronLevel = ACT1AboveElectronLevel && ( (intCharge[j][0] + intCharge[j+1][0]) > ACT1ElectronThreshC);
      }

      if (j==2){
        noPeaksInACT1 = noPeaksInACT1 && (nPeaks[j] == 0 || nPeaks[j+1] == 0);

      }

      // if (j==signalTime[13][0]){
      //   ACT1AboveElectronLevel = ACT1AboveElectronLevel && ( (intCharge[j][0] + intCharge[j+1][0]) > ACT1ElectronThreshC);
      // }


    } // channels


    if (peakMode == "a" && ! (onePeakInAll) )
      continue;
    if (peakMode == "b" && (! moreThanOnePeakInAllACTs) )
      continue;
    if (peakMode == "c" && ! (onePeakInAllACTs && moreThanOnePeakInAllToFs) )
      continue;
    if (peakMode == "d" && ! (moreThanOnePeakInAllACTs && onePeakInAllToFs) )
      continue;
    if (peakMode == "e" && ! (onePeakInAllACTs) )
      continue;
    if (peakMode == "f" && ! (onePeakInAllToFs) )
      continue;
    if (peakMode == "g" && ( ! (onePeakInAllToFs) || ! (ACT23AboveElectronLevel) || ! (PbGlassAboveElectronLevel)))
      continue;
    if (peakMode == "h"  && ( ! (onePeakInAllToFs) || ACT1AboveElectronLevel))
      continue; //we continue ie do not take the data if we are above the level
    if (peakMode == "i"  && ( ! (onePeakInAllToFs) || !ACT1AboveElectronLevel))
      continue; //we continue ie do not take the data if we are below the level
    if (peakMode == "j"  && ( ! (onePeakInAllToFs) || !noPeaksInACT1))
      continue; //we continue ie do not take the data if we are below the level
    if (peakMode == "k"  && ( ! (onePeakInAllToFs) || noPeaksInACT1))
      continue; //we continue ie do not take the data if we are below the level


    for(int j = 0; j < nChannels; j++) {


      /*
       *	int ind = std::max_element(peakVoltage->at(j).begin(),peakVoltage->at(j).end()) - peakVoltage->at(j).begin();
       *	indices.at(j) = ind;
       *	hCharge.at(j).Fill(intCharge->at(j).at(indices.at(j)));
       *	hVoltage.at(j).Fill(peakVoltage->at(j).at(indices.at(j)));
       *	hTime.at(j).Fill(signalTime->at(j).at(indices.at(j)));
       *	hHit.at(j).Fill(peakVoltage->at(j).size());
       *	hPedestalSigma.at(j).Fill(pedestalSigma[j]);
       */
      hCharge.at(j).Fill(intCharge[j][0]);
      hVoltage.at(j).Fill(peakVoltage[j][0]);
      hTime.at(j).Fill(signalTime[j][0]);
      //hNbPeaks.at(j).Fill(nPeaks[j]);
      hPedestalSigma.at(j).Fill(pedestalSigma[j]);
      hnPeaks.at(j).Fill(nPeaks[j]);
      hPedestal.at(j).Fill(pedestal[j]);
      hPedestalNbPeaks.at(j).Fill(pedestal[j], nPeaks[j]);
    }



    // JK's time resolution of 2022
    double t0a = (signalTime[8][0]  + signalTime[11][0] ) / 2.;
    double t0b = (signalTime[9][0]  + signalTime[10][0] ) / 2.;
    double t1a = (signalTime[12][0] + signalTime[15][0] ) / 2.;
    double t1b = (signalTime[13][0] + signalTime[14][0] ) / 2.;

    // time diffs for time resolution histogramme 2022:
    double t0diff = t0a - t0b;
    double t1diff = t1a - t1b;

    // jiri
    // Fill resolution histograms of 2022:
    hTimeReso0.Fill(t0diff);
    hTimeReso1.Fill(t1diff);
    hTimeReso0_zoom.Fill(t0diff);
    hTimeReso1_zoom.Fill(t1diff);


    double t00 = signalTime[8][0];
    double t01 = signalTime[9][0];
    double t02 = signalTime[10][0];
    double t03 = signalTime[11][0];
    hTimeDiffTOF01.Fill(t00 - t01); //carefull, this particular histogram is a difference of hit times,
    hTimeDiffTOF02.Fill(t00 - t02); //not a TOF per say but instead the cable length difference
    hTimeDiffTOF03.Fill(t00 - t03); //plus the photon travel time though the panel

    double t10 = signalTime[12][0];
    double t11 = signalTime[13][0];
    double t12 = signalTime[14][0];
    double t13 = signalTime[15][0];
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

    double t0 = ( signalTime[8][0] +  signalTime[9][0] + signalTime[10][0] + signalTime[11][0]) / 4.;
    double t1 = (signalTime[12][0] + signalTime[13][0] + signalTime[14][0] + signalTime[15][0]) / 4.;
    double tof = t1-t0;

    double act0a = peakVoltage[0][0] + peakVoltage[1][0];
    double act1a = peakVoltage[2][0] + peakVoltage[3][0];
    double act2a = peakVoltage[4][0] + peakVoltage[5][0];
    double act3a = peakVoltage[6][0] + peakVoltage[7][0];

    double hc0a = peakVoltage[16][0];
    double hc1a = peakVoltage[17][0];
    double pba = peakVoltage[18][0];

    hTOFACT0A.Fill(tof, act0a);
    hTOFACT1A.Fill(tof, act1a);
    hTOFACT2A.Fill(tof, act2a);
    hTOFACT3A.Fill(tof, act3a);

    hPbACT23A.Fill(pba, (act2a+act3a) / 2.);
    hPbACT1A.Fill(pba, (act1a));
    hTOFACT23A.Fill(tof, (act2a+act3a) / 2.);
    hPbATOF.Fill(pba, tof);

    double act0c = intCharge[0][0] + intCharge[1][0];
    double act1c = intCharge[2][0] + intCharge[3][0];
    double act2c = intCharge[4][0] + intCharge[5][0];
    double act3c = intCharge[6][0] + intCharge[7][0];


    double hc0c = intCharge[16][0];
    double hc1c = intCharge[17][0];
    double pbc = intCharge[18][0];

    hTOFACT0C.Fill(tof, act0c);
    hTOFACT1C.Fill(tof, act1c);
    hTOFACT2C.Fill(tof, act2c);
    hTOFACT3C.Fill(tof, act3c);

    hPbCCT1C.Fill(pbc, (act1c));

    hPbACT23C.Fill(pbc, (act2c+act3c) / 2.);
    hTOFACT23C.Fill(tof, (act2c+act3c) / 2.);
    hPbCTOF.Fill(pbc, tof);

    hACT1CACT3C.Fill(act1c, act3c);
    hACT3CACT2C.Fill(act3c, act2c);
    hACT2CACT1C.Fill(act2c, act1c);

    //acraplet - weird electrons which do not see anthing in the ACT
    if ((act2a+act3a) / 2. != 1.5 && tof >= 13.5 && tof <= 16.5) {
      hHC0AHC1A.Fill(hc0a, hc1a);
      hHC0CHC1C.Fill(hc0c, hc1c);
    }

    //acraplet
    hTOFAll.Fill(tof);
    hTOFAllLow.Fill(tof);

    hT0.Fill(t0);
    hT1.Fill(t1);


    // 2D nPeaks
    double nPeaksToFAver = (nPeaks[8] + nPeaks[9] + nPeaks[10] + nPeaks[11] + nPeaks[12] + nPeaks[13] + nPeaks[14] + nPeaks[15]) / 8;
    double nPeaksACT23Aver = (nPeaks[4] + nPeaks[5] + nPeaks[6] + nPeaks[7]) / 4;

    double nPeaksToF0Aver = (nPeaks[8] + nPeaks[9] + nPeaks[10] + nPeaks[11] ) / 4;
    double nPeaksToF1Aver = (nPeaks[12] + nPeaks[13] + nPeaks[14] + nPeaks[15]) / 4;
    double nPeaksACT2Aver = (nPeaks[4] + nPeaks[5] ) / 2;
    double nPeaksACT3Aver = (nPeaks[6] + nPeaks[7]) / 2;


    hnPeaksACT23vsnPeaksToF.Fill(nPeaksToFAver, nPeaksACT23Aver);

    hnPeaksToF1vsnPeaksToF0.Fill(nPeaksToF1Aver, nPeaksToF0Aver);
    hnPeaksACT3vsnPeaksACT2.Fill(nPeaksACT2Aver, nPeaksACT3Aver);

    hnPeaksACT23vsToF.Fill(tof,  nPeaksACT23Aver);
    hnPeaksACT23vsToFlow.Fill(tof,  nPeaksACT23Aver);
    hnPeaksToFvsToF.Fill(tof, nPeaksToFAver);
    hnPeaksToFvsToFlow.Fill(tof, nPeaksToFAver);
    hnPeaksLeadGlassvsLeadGlassA.Fill(pba, nPeaks[18]);

    hnPeaksACT23vsLeadGlassA.Fill(pba, nPeaksACT23Aver);
    hnPeaksToFvsLeadGlassA.Fill(pba, nPeaksToFAver);

    // selections:

    bool pass = true;
    bool isEl = false;

    bool isElACT23pb = false;
    bool isMuACT23pb = false;
    bool isPiACT23pb = false;
    bool ispACT23pb = false;
    bool isdACT23pb = false;



    bool passed_act1a_cuts = false;
    bool passed_act2a_cuts = false;

    switch(momentum)	  {

    //new populations for number of particle analysis

      case 320: { //; jk
        // placeholder for momentum dependent cuts
        break;
      } // 320

      case 900: {

        if (tof >= y0_cut && onePeakInAllToFs) {
          isdACT23pb = true;
        }

        else if (tof >= x0_cut && onePeakInAllToFs){
          ispACT23pb = true;
        }

        else if ((act2a+act3a)/2 < act23_pi_minA && tof <= x0_cut && onePeakInAllToFs) {
          if (tof >= 13){
            std::cout << tof << std::endl;
          }
          isMuACT23pb = true;
        }

        else if ((act2a+act3a)/2 > act23_pi_minA && tof <= x0_cut && onePeakInAllToFs) {
          isElACT23pb = true;
        }


      }



      default: {
        // if (i < 10)
        //   cout << "WARNING: Using default settings for the " << momentum << " MeV/c beam" << endl;
        // if ( pbc > 0.9 || pbc < 0.1 || act1c > 0.25) { // custom electron removal cut (act2a + act3a)/2.
        //   isEl = true;
        // }
        //
        // if ((act2a+act3a)/2 > y0_cut - y0_cut/x0_cut * pba && pba > pb_min && onePeakInAllToFs) {
        //   std::cout << onePeakInAllToFs << std::endl;
        //   isElACT23pb = true;
        // }
        //
        // else if ( (act2a+act3a)/2 > act23_pi_maxA && pba > pb_min && onePeakInAllToFs) {
        //   isMuACT23pb = true;
        //
        // }
        // else if ( (act2a+act3a)/2 > act23_pi_minA && pba > pb_min && onePeakInAllToFs) {
        //   isPiACT23pb = true;
        //
        // }
        //
        // else if ( pba < pb_min && !onePeakInAllToFs) {
        //   ispACT23pb = true;
        // }
        //

      } // default
    } // case

    if (!pass) continue;

    if (isEl) {
      // electrons
      hTOFEl.Fill(tof);
      hTOFElLow.Fill(tof);
    }

    if (isElACT23pb) {
      // electrons
      hTOFElACT23pb.Fill(tof);
      //hTOFElLow.Fill(tof);
    }

    if (isMuACT23pb) {
      // muons
      hTOFMuACT23pb.Fill(tof);
      //hTOFElLow.Fill(tof);
    }

    if (isPiACT23pb) {
      // pions
      hTOFPiACT23pb.Fill(tof);
      //hTOFElLow.Fill(tof);
    }

    if (ispACT23pb) {
      hTOFpACT23pb.Fill(tof);
    }

    if (isdACT23pb) {
      hTOFdACT23pb.Fill(tof);
    }

    if (!isEl || !isElACT23pb) {
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





