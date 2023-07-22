
// WCTE TB July 2023
// 20th July 2023
// author mh
// updated by jk


#include "runs.h"



void compare_LG(){

  //int run_numbers[] = {435,436,437,449,438,439,440,441}; ///< n=1.02, positive
  //int run_numbers[] = {399,440};
  //int run_numbers[] = {447,446,445,444,443,442};  ///<n=1.02, negative
  //int run_numbers[] = {459, 458, 457, 456, 455, 454, 453, 451, 452}; ///< n=1.015, positive

  //int run_numbers[] = {460, 461, 462, 463, 469, 470, 471, 459, 458, 457}; ///< n=1.015, negative
  int run_numbers[] = {490, 491, 492};
  
  
  int lsts[] = {1,1,1,2,2,2,
		1,1,1,1,1,1,1,1,1,1,1,1,1};

  int cols[] = {1, kGreen+2, 4, 1, kGreen+2, 4, 3, 5, 6, 7, 8};
  
  
  int nfiles = sizeof(run_numbers) / sizeof(int);

    
  TFile *fi[nfiles];
  TH1D *hists[nfiles];

  TLegend *leg = new TLegend(0.60,0.6,0.88,0.88);
  leg->SetBorderSize(0.);
  leg->SetFillColor(0);

	
  TString hnames[] = {"DT730_02_Amplitude_2",
		      "DT730_01_Amplitude_4",
		      "DT730_00_Amplitude_0",
		      "DT730_00_Amplitude_1"};
  int nhistos = sizeof(hnames) / sizeof(TString);
  // HACK
  //nhistos = 1;
	
  TString path = "/media/wcte/T7/WCTE/data/2023/hist_files/";

  TCanvas *cans[nhistos];

  int ndraw = 0;
  for (int ih = 0; ih < nhistos; ++ih) {

    TString hname = hnames[ih];
    TString canname = "can_" + hname;
    cans[ih] = new TCanvas(canname, canname, ih*50, ih*50, 1000, 800);

	       
    for (int i=0; i < nfiles; i++) {

      TString fname = Form("%soutput00000%d.root", path.Data(), run_numbers[i]);
      fi[i] = new TFile(fname);
      if (!fi[i] || fi[i]->IsZombie()) {
	cerr << "ERROR opening file " << fname.Data() << endl;
	continue;
      }

      hists[i] = (TH1D*)fi[i]->Get(hname);
      if (!hists[i]) {
	cerr << "ERROR getting histo " << hname.Data() << endl;
	continue;
      }
      
      //hists[i]->SetLineColor(i+1);

      hists[i]->SetLineColor(cols[i]);
      hists[i]->SetLineWidth(3);
      hists[i]->SetLineStyle(lsts[i]);
      
      hists[i]->Rebin(4);
      int b1 = 0;
      int b2 = hists[i] -> GetNbinsX() + 1;
      if (hname.Contains("DT730_02_Amplitude_2"))
      	b1 = hists[i] -> GetXaxis()->FindBin(60.);
      hists[i]->Scale(1.0/hists[i]->Integral(b1, b2));

      hists[i]->SetStats(false);

      
      if (i==0) hists[i]->Draw("hist");
      else hists[i]->Draw("same hist");

      if (ih == 0) {
	leg->AddEntry(hists[i],Form("Run %d %+i: MeV/c",run_numbers[i], runsDict[run_numbers[i]]),"l");
      }
      ndraw++;
    } // files
  } // histos


  gSystem -> Exec("mkdir -p png/ pdf/");
  for (auto can : cans) {
    can -> cd();
    leg -> Draw();
    if (false) {
      gPad -> SetLogy(1);	
      can -> Print(TString("png/") + TString(can -> GetName()) + "_logy.png");
      can -> Print(TString("png/") + TString(can -> GetName()) + "_logy.png");
      can -> Print(TString("pdf/") + TString(can -> GetName()) + "_logy.png");
      can -> Print(TString("pdf/") + TString(can -> GetName()) + "_logy.png");
      gPad -> SetLogy(0);
      can -> Print(TString("png/") + TString(can -> GetName()) + "_liny.png");
      can -> Print(TString("png/") + TString(can -> GetName()) + "_liny.png");
      can -> Print(TString("pdf/") + TString(can -> GetName()) + "_liny.png");
      can -> Print(TString("pdf/") + TString(can -> GetName()) + "_liny.png");
    }
  }

}
