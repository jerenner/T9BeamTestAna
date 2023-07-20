#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "../include/Enumerators.h"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <numeric>

using namespace std;

// Junjie Xia 2023
// macro to find the waveforms of events
///// 1. in different regions in the act vs Pb amplitude plane
///// or 2. with a given event number

// prescribed from MakeDataPlots.C and waveform_analysis.cc

const int nMaxChannels = 32;
const int nChannelsPerDigi = 8;
const double nanosecsPerSample = 2;
//const int nChannels = 19; // UPDATE THIS FOR HODOSCOPE PMTs! to e.g. 32!

const string rootfile_dir = "root_files";
const string root_treename = "midas_data_D30";
const string ntuplefile_dir = "ntuple_files";
const string ntuple_treename = "anaTree";


bool LinkDigiDet(){
  //fill digitizer channel map
  for (int ich = 0; ich < nMaxChannels; ich++){
    pair<int, int> dpair(ich/nChannelsPerDigi, ich%nChannelsPerDigi);    
    vec_Digitizers.push_back(dpair);
  }

  for (int idet = 0; idet < (int)DetChannels::nDetChannels; idet++){
    DCtoDGmap[DetChannels(idet)] = vec_Digitizers[idet];
  }
  return true;  
}


int EventSelection(double act, double pbg){
    // in the act2+3 vs. pb amplitude space
    // pb < 1.3 && act < 10 -> mu/pi
    // 1.3 < pb < 2.5 && 10 < act -> e
    // pb > 2.5 -> mu/pi in e pile-up?

  if (pbg < 1.3 && act < 10){
    return (int)EventRegion::mupi_1pulse;
  }
  else if (pbg > 1.3 && pbg < 2.5 && act > 10){
    return (int)EventRegion::e_npulses;
  }
  else if (pbg > 2.5){
    return (int)EventRegion::mupi_npulses;
  }
  else{
    return (int)EventRegion::not_defined;
  }
}

void CheckWaveForms(string data_dir, int runnum, int det_to_plot, const int nplots = 10, const int spillID = -1, const int eventID = -1) {
//void CheckWaveForms() {  
  if (!LinkDigiDet()) {
    cerr << "Cannot link digitizer channels to detectors" << endl;
    exit(-1);
  }
  
  TFile ntupleFile(Form("%s/%s/ntuple_%06d.root", data_dir.c_str(), ntuplefile_dir.c_str(), runnum), "READ"); // digitized events  
  TFile rootFile(Form("%s/%s/root_run_%06d.root", data_dir.c_str(), rootfile_dir.c_str(), runnum), "READ");

  if (ntupleFile.IsZombie() || rootFile.IsZombie() || det_to_plot < (int)DetChannels::act_e_veto0 || det_to_plot >= (int)DetChannels::nDetChannels || spillID*eventID < 0) {
    cout << "<<<<<USAGE>>>>>: " << endl;
    cout << "root -l -b -q 'CheckWaveForms(\"/your/data/dir/that/has/ntuple_files/and/root_files\", run number, detector id to show the waveform (0~18), number of plots in each category (default = 10), (optional) spill number to plot, (optional) event number to plot)'" << endl;
    exit(0);
  }

  bool foundspill = false;
  bool foundevent = false;
  
  gSystem->Exec(Form("mkdir -p %s/plots", data_dir.c_str()));//create plot dir in the data dir
  gSystem->Exec(Form("mkdir -p %s/root_hists", data_dir.c_str()));//create dir for root files of hists in the data dir
  
  Double_t peakVoltage[nMaxChannels][1];
  vector<double> *waveform = NULL;
  int eventNumber, spillNumber;

  pair<int, int> digipair = DCtoDGmap.at(DetChannels(det_to_plot));//digitizer ID and channels

  TTree *tree_ntuple = (TTree*)ntupleFile.Get(ntuple_treename.c_str());
  std::vector<TTree*> vec_root_trees;
  for (int ich = 0; ich < nMaxChannels/nChannelsPerDigi; ich++){
    vec_root_trees.push_back((TTree*)rootFile.Get(Form("%s%d", root_treename.c_str(), ich)));
  }

  int nevt = tree_ntuple->GetEntries();
  assert(nevt == vec_root_trees[digipair.first]->GetEntries());
  tree_ntuple->AddFriend(vec_root_trees[digipair.first]); // only add one digitizer branch
  
  tree_ntuple->SetBranchAddress("PeakVoltage",&peakVoltage);
  tree_ntuple->SetBranchAddress(Form("Channel%d", digipair.second), &waveform);
  tree_ntuple->SetBranchAddress("eventNumber", &eventNumber);
  tree_ntuple->SetBranchAddress("spillNumber", &spillNumber);

  vector<int> nplots_region(int(EventRegion::EventRegionCount), 0);
  TH1D* hist_waveform_cuts[int(EventRegion::EventRegionCount)][nplots]; //waveform for events selected by cuts 
  TH1D* hist_waveform;// waveform for events specified in the input

  for (int ievt = 0; ievt < nevt; ievt++){
      tree_ntuple->GetEntry(ievt);
      if (spillID >= 0 && eventID >= 0){
	if (spillNumber != spillID) { // selected spill
	  continue;
	}
	else {foundspill = true;}
	if (eventNumber != eventID){ // selected event
	  continue;
	}
	else {foundevent = true;}
      }
      double actamp = 0.5*(peakVoltage[(int)DetChannels::act_2_R][0] +
			   peakVoltage[(int)DetChannels::act_2_R][0] +
			   peakVoltage[(int)DetChannels::act_3_R][0] +
			   peakVoltage[(int)DetChannels::act_3_R][0]);
      
      double pba = peakVoltage[(int)DetChannels::pbg][0];
      int iregion = EventSelection(actamp, pba); // find the category given act and pb amplitudes
      
      if ((spillID < 0 || eventID < 0) && (iregion < 0 || nplots_region[iregion] >= nplots)){
	//not interested in these events out of the cut regions
	continue;
      }

      int nsamples = waveform->size();
      //get vector of x-axis values
      vector<double> fill_bins(nsamples);
      int start_bin = 1;
      std::generate(fill_bins.begin(), fill_bins.end(), [&start_bin]{ return start_bin+=nanosecsPerSample; });

      if (spillID < 0 || eventID < 0){
	char* histname_cuts = Form("Spill%d_Evt%d_%s_%s", spillNumber, eventNumber, EventRegionNames.at(EventRegion(iregion)).c_str(), DetChNames.at(DetChannels(det_to_plot)).c_str());//decimals are not accepted in hist name
	char* histtitle_cuts = Form("Spill%d_Evt%d_%s_%s_ACT%.2f_Pb%.2f", spillNumber, eventNumber, EventRegionNames.at(EventRegion(iregion)).c_str(), DetChNames.at(DetChannels(det_to_plot)).c_str(), actamp, pba);
		
	//initialize hists
	hist_waveform_cuts[iregion][nplots_region[iregion]] =
	  new TH1D(histname_cuts, histtitle_cuts, nsamples, 0, nsamples*nanosecsPerSample);
	hist_waveform_cuts[iregion][nplots_region[iregion]]->SetDirectory(0);
	hist_waveform_cuts[iregion][nplots_region[iregion]]->FillN(nsamples, &fill_bins[0], &(*waveform)[0]);
	nplots_region[iregion]++;
      }
      else{
	char* histname = Form("Spill%d_Evt%d_%s", spillNumber, eventNumber, DetChNames.at(DetChannels(det_to_plot)).c_str());
	char* histtitle = Form("Spill%d_Evt%d_%s_ACT%.2f_Pb%.2f", spillNumber, eventNumber, DetChNames.at(DetChannels(det_to_plot)).c_str(), actamp, pba);
	hist_waveform = new TH1D(histname, histtitle, nsamples, 0, nsamples*nanosecsPerSample);
	hist_waveform->SetDirectory(0);
	hist_waveform->FillN(nsamples, &fill_bins[0], &(*waveform)[0]);
      }
      
  }

  
  if (spillID < 0 || eventID < 0){
    for (int iregion = 0; iregion < (int)EventRegion::EventRegionCount; iregion++){
      printf("Number of plots plotted in region %s is %d.\n", EventRegionNames.at(EventRegion(iregion)).c_str(), nplots_region[iregion]);
    }
  }
  else if (!(foundspill&&foundevent)){
    printf("Cannot find spill %d event %d in run %d.\n", spillID, eventID, runnum);
    exit(1);
  }
 
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  if (eventID < 0 || spillID < 0){//events selected by cuts
    c1->Print(Form("%s/plots/%s_waveform_cuts_%d_evts.pdf[", data_dir.c_str(), DetChNames.at(DetChannels(det_to_plot)).c_str(), nplots));
    TFile *outfile = new TFile(Form("%s/root_hists/%s_waveform_cuts_%d_evts.root", data_dir.c_str(), DetChNames.at(DetChannels(det_to_plot)).c_str(), nplots),"recreate");
    for (int iregion = 0; iregion < (int)EventRegion::EventRegionCount; iregion++){
      for (int iplot = 0; iplot < nplots; iplot++){
	outfile->cd();
	hist_waveform_cuts[iregion][iplot]->Write();
	c1->cd();	
	hist_waveform_cuts[iregion][iplot]->Draw("hist c");
	c1->Update();
	auto stat = dynamic_cast<TPaveStats*>(hist_waveform_cuts[iregion][iplot]->FindObject("stats"));
	if (stat) {
	  stat->SetX1NDC(0.7);
	  stat->SetX2NDC(0.9);
	  stat->SetY1NDC(0.1);
	  stat->SetY2NDC(0.3);
	  stat->Draw();
	}
	c1->Modified();
	c1->Update();
	c1->Print(Form("%s/plots/%s_waveform_cuts_%d_evts.pdf", data_dir.c_str(), DetChNames.at(DetChannels(det_to_plot)).c_str(), nplots));
      }    
    }    
    c1->Print(Form("%s/plots/%s_waveform_cuts_%d_evts.pdf]", data_dir.c_str(), DetChNames.at(DetChannels(det_to_plot)).c_str(), nplots));
    outfile->Close();
  }
  else{
    TFile *outfile = new TFile(Form("%s/root_hists/%s_spill_%d_evt_%d_waveform.root", data_dir.c_str(), DetChNames.at(DetChannels(det_to_plot)).c_str(), spillID, eventID), "recreate");
    outfile->cd();
    hist_waveform->Write();
    c1->cd();
    hist_waveform->Draw("hist c");
    c1->Update();
    auto stat = dynamic_cast<TPaveStats*>(hist_waveform->FindObject("stats"));
    if (stat) {
      stat->SetX1NDC(0.7);
      stat->SetX2NDC(0.9);
      stat->SetY1NDC(0.1);
      stat->SetY2NDC(0.3);
      stat->Draw();
    }    
    c1->Modified();
    c1->Update();
    c1->Print(Form("%s/plots/%s_spill_%d_evt_%d_waveform.png", data_dir.c_str(), DetChNames.at(DetChannels(det_to_plot)).c_str(), spillID, eventID));
    c1->Print(Form("%s/plots/%s_spill_%d_evt_%d_waveform.pdf", data_dir.c_str(), DetChNames.at(DetChannels(det_to_plot)).c_str(), spillID, eventID));
    outfile->Close();
  }
  
}

