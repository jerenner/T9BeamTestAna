#define EventInfo_cxx
#define PMT_cxx

#include "EventInfo.h"
#include "PMT.h"

void EventInfo::Loop() {}
void PMT::Loop() {}

void triggerTimeDrift(string input = "/neut/datasrv2a/jrenner/ntuple_files/ntuple_000435.root",
                      string slopes_file = "triggerTimeDrift.txt")
{

  gStyle->SetOptStat(1);
  gStyle->SetFuncWidth(1);
  gStyle->SetMarkerStyle(6);

  bool saveplots = false;

  // input file
  TFile * f = new TFile(input.c_str());
  cout << f->GetName() << endl;

  //extract run number
  string run = input.substr(input.size()-8,3);
  string plotdir(".");
  string plotout(Form("%s/run%s",plotdir.c_str(),run.c_str()));
  string title(Form("Run %s",run.c_str()));
  cout << "run number " << run << endl;

  // input trees
  vector<string> tree_name;
  tree_name.push_back("ACT0L");//0
  tree_name.push_back("ACT0R");
  tree_name.push_back("ACT1L");//2
  tree_name.push_back("ACT1R");
  tree_name.push_back("ACT2L");//4
  tree_name.push_back("ACT2R");
  tree_name.push_back("ACT3L");//6
  tree_name.push_back("ACT3R");
  tree_name.push_back("TOF00");//8
  tree_name.push_back("TOF01");
  tree_name.push_back("TOF02");//10
  tree_name.push_back("TOF03");
  tree_name.push_back("TOF10");//12
  tree_name.push_back("TOF11");
  tree_name.push_back("TOF12");//14
  tree_name.push_back("TOF13");
  tree_name.push_back("Hole0");//16
  tree_name.push_back("Hole1");
  tree_name.push_back("PbGlass");//18
  tree_name.push_back("EventInfo");
  const int npmts = tree_name.size()-1;

  const int ndigitizers = 3;
  unsigned long offset = 2147483648;//2^31 steps = ~17s
  double to_s  = 8e-9;

  // get trees
  EventInfo * info;
  vector<PMT*> pmt;
  vector<TTree*> tree;
  for (int itree=0; itree<tree_name.size(); itree++) {
    TTree * t = (TTree*)f->Get(tree_name[itree].c_str());
    tree.push_back(t);
    if (itree<npmts) {
      pmt.push_back(new PMT(t));
    }
    else {
      info = new EventInfo(t);
    }
  }

  // all digitizers should have
  // the same number of entries
  int nentries = tree[0]->GetEntries();
  cout << "entries " << nentries << endl;
  for (int itree=0; itree<tree_name.size(); itree++) {
    if (tree[itree]->GetEntries()!=nentries ) {
      cout << "different entries in " << tree[itree]->GetName();
      cout << " " << tree[itree]->GetEntries() << endl;
    }
  }

  // plots
  TGraph * gdrift[ndigitizers];
  TGraph * gdriftcor[ndigitizers];
  TGraph * gdriftspill[ndigitizers];
  TGraph * gdriftspillsingle[ndigitizers];
  for (int dg=0; dg<ndigitizers; dg++) {
    gdrift[dg] = new TGraph();
    gdriftcor[dg] = new TGraph();
    gdriftspill[dg] = new TGraph();
    gdriftspillsingle[dg] = new TGraph();
  }
  TH1D * hslope[ndigitizers];
  for (int dg=0; dg<ndigitizers; dg++) {
    if (dg==0) {
      hslope[dg] = new TH1D(Form("hslope%i",dg),
                            Form(";Digitizer %i slope (us/s)",dg),
                            100,-2.08,-2.0);
    }
    else {
      hslope[dg] = new TH1D(Form("hslope%i",dg),
                            Form(";Digitizer %i slope (us/s)",dg),
                            100,0.5,0.68);
    }
  }

  // save for each spill
  vector< map<int,double> > slopespill(ndigitizers);

  // store previous trigger information
  unsigned int spillNumber_pre = 0;
  vector<unsigned long> timeStamp(ndigitizers,0);
  vector<unsigned long> timeStamp_pre(ndigitizers,0);
  vector<unsigned long> timeStamp_init(ndigitizers,0);
  vector<unsigned long> triggerTime(ndigitizers,0);
  vector<unsigned long> triggerTime_pre(ndigitizers,0);
  vector<unsigned long> triggerTime_full_pre(ndigitizers,0);
  vector<unsigned long> triggerTime_off(ndigitizers,0);
  vector<unsigned long> triggerTime_init(ndigitizers,0);
  vector<unsigned long> triggerTime_spill(ndigitizers,0);

  // loop over entries
  for (int ientry=0; ientry<nentries; ientry++) {

    // load all trees
    info->GetEntry(ientry);
    for (int ipmt=0; ipmt<npmts; ipmt++) {
      pmt[ipmt]->GetEntry(ientry);
    }

    // calculate slope for previous spill
    if (info->EventNumber==0 && info->SpillNumber>0) {

     // for (int d=0; d<ndigitizers; d++) {

      if (gdriftspill[0]->GetN()>100) {
        gdriftspill[0]->Fit("pol1","Q");
        TF1 * fitspill0 = gdriftspill[0]->GetFunction("pol1");
        slopespill[0][spillNumber_pre] = fitspill0->GetParameter(1);
        hslope[0]->Fill(fitspill0->GetParameter(1));
        delete fitspill0;

      }

      //slopespill[1][spillNumber_pre] = 0;

      if (gdriftspill[2]->GetN()>100) {
        gdriftspill[2]->Fit("pol1","Q");
        TF1 * fitspill2 = gdriftspill[2]->GetFunction("pol1");
        slopespill[2][spillNumber_pre] = fitspill2->GetParameter(1);
        hslope[2]->Fill(fitspill2->GetParameter(1));
        delete fitspill2;
      }

      delete gdriftspill[0];
      delete gdriftspill[1];
      delete gdriftspill[2];

      gdriftspill[0] = new TGraph();
      gdriftspill[1] = new TGraph();
      gdriftspill[2] = new TGraph();

    }

    // trigger time and time stamp for each digitizer
    triggerTime[0] = pmt[0] ->triggerTime;
    triggerTime[1] = pmt[8] ->triggerTime;
    triggerTime[2] = pmt[16]->triggerTime;

    timeStamp[0] = pmt[0] ->timeStamp;
    timeStamp[1] = pmt[8] ->timeStamp;
    timeStamp[2] = pmt[16]->timeStamp;

    for (int dg=0; dg<ndigitizers; dg++) {

      // initilize
      if (ientry==0) {
        cout << "digitizer " << dg;
        cout << " first TT "     << triggerTime[dg];
        cout << " TT (s) " << triggerTime[dg]*to_s;
        cout << endl;
        spillNumber_pre          = info->SpillNumber;
        timeStamp_pre[dg]        = timeStamp[dg];
        timeStamp_init[dg]       = timeStamp[dg];
        triggerTime_pre[dg]      = triggerTime[dg];
        triggerTime_full_pre[dg] = triggerTime[dg];
        triggerTime_init[dg]     = triggerTime[dg];
        triggerTime_spill[dg]    = triggerTime[dg];
      }

      // trigger time resets after 2^31
      // so add an offset for every time it resets
      //if (triggerTime[dg]<triggerTime_pre[dg]) {
      //  triggerTime_off[dg] += offset;
      //}

      // increase trigger time offset until the total trigger time
      // difference matches the time stamp difference
      unsigned long diffTT = triggerTime[dg]+triggerTime_off[dg]-triggerTime_full_pre[dg];
      unsigned long diffTS = timeStamp[dg]-timeStamp_pre[dg];
      while (fabs(diffTS-diffTT*to_s)>2) {
        triggerTime_off[dg] += offset;
        diffTT = triggerTime[dg]+triggerTime_off[dg]-triggerTime_full_pre[dg];
      }

      // update spill time
      // if this is a new spill
      if (info->SpillNumber!=spillNumber_pre) {
        triggerTime_spill[dg] = triggerTime[dg] + triggerTime_off[dg];
      }

    }//ndigitizer

    // add total offset to trigger time
    unsigned long triggerTimeFull[3];
    triggerTimeFull[0] = triggerTime[0] + triggerTime_off[0];
    triggerTimeFull[1] = triggerTime[1] + triggerTime_off[1];
    triggerTimeFull[2] = triggerTime[2] + triggerTime_off[2];

    // trigger time from run start
    unsigned long triggerTimeFullInit[3];
    triggerTimeFullInit[0] = triggerTimeFull[0] - triggerTime_init[0];
    triggerTimeFullInit[1] = triggerTimeFull[1] - triggerTime_init[1];
    triggerTimeFullInit[2] = triggerTimeFull[2] - triggerTime_init[2];

    // trigger time from run start in seconds
    double timeFullInit[3];
    timeFullInit[0] = triggerTimeFullInit[0] * to_s;
    timeFullInit[1] = triggerTimeFullInit[1] * to_s;
    timeFullInit[2] = triggerTimeFullInit[2] * to_s;

    // trigger time from spill start
    unsigned long triggerTimeFullSpill[3];
    triggerTimeFullSpill[0] = triggerTimeFull[0] - triggerTime_spill[0];
    triggerTimeFullSpill[1] = triggerTimeFull[1] - triggerTime_spill[1];
    triggerTimeFullSpill[2] = triggerTimeFull[2] - triggerTime_spill[2];

    // trigger time from spill start in seconds
    double timeFullSpill[3];
    timeFullSpill[0] = triggerTimeFullSpill[0] * to_s;
    timeFullSpill[1] = triggerTimeFullSpill[1] * to_s;
    timeFullSpill[2] = triggerTimeFullSpill[2] * to_s;

    for (int dg=0; dg<ndigitizers; dg++) {

      // drift since run start
      long ttFullInitDiff  = triggerTimeFullInit[dg]  - triggerTimeFullInit[1];
      double tFullInitDiff  = ttFullInitDiff  * to_s * 1e6;//us
      gdrift[dg]->SetPoint(gdrift[dg]->GetN(),timeFullInit[1],tFullInitDiff);

      // drift since spill start
      long ttFullSpillDiff = triggerTimeFullSpill[dg] - triggerTimeFullSpill[1];
      double tFullSpillDiff = ttFullSpillDiff * to_s * 1e6;//us
      gdriftspill[dg]->SetPoint(gdriftspill[dg]->GetN(),timeFullSpill[1],tFullSpillDiff);

      // 1 spill
      if (info->SpillNumber==1) {
        gdriftspillsingle[dg]->SetPoint(gdriftspillsingle[dg]->GetN(),timeFullSpill[1],tFullSpillDiff);
      }

    }//digitizers

    // save times for next loop
    for (int dg=0; dg<ndigitizers; dg++) {
      timeStamp_pre[dg] = timeStamp[dg];
      triggerTime_pre[dg] = triggerTime[dg];
      triggerTime_full_pre[dg] = triggerTime[dg] + triggerTime_off[dg];
    }//digitizers

    // save spill number
    if (info->SpillNumber!=spillNumber_pre) {
      spillNumber_pre = info->SpillNumber;
    }

  }// end input tree

  cout << "last spill number " << info->SpillNumber << endl;
  for (int dg=0; dg<ndigitizers; dg++) {
    cout << "digitizer " << dg;
    cout << " TT since run start (s) ";
    cout << (triggerTime[dg]+triggerTime_off[dg]-triggerTime_init[dg])*to_s;
    cout << endl;
  }//digitizers

  double slope[ndigitizers] = {0,0,0};
  double period[ndigitizers] = {8,8,8};

  for (int dg=0; dg<ndigitizers; dg++) {

    if (dg==1) continue;

    // draw plots
    TCanvas * c;
    double min, max;

    // all spills
    c = new TCanvas();      
    gdrift[dg]->Draw("ap");
    gdrift[dg]->GetYaxis()->SetTitle(Form("TT%i-TT1 (us)",dg));
    gdrift[dg]->GetXaxis()->SetTitle("TT1 (s)");
    min = gdrift[dg]->GetHistogram()->GetMinimum();
    max = gdrift[dg]->GetHistogram()->GetMaximum();
    gdrift[dg]->GetYaxis()->SetRangeUser(min,max*(1.4));
    gdrift[dg]->SetTitle(Form("%s",title.c_str()));
    gdrift[dg]->Fit("pol1","Q");
    gdrift[dg]->Draw("p");
    if (saveplots) c->Print(Form("%s_drift_%i.png",plotout.c_str(),dg));

    // last spill
    c = new TCanvas();      
    gdriftspill[dg]->Draw("ap");
    gdriftspill[dg]->GetYaxis()->SetTitle(Form("TT%i-TT1 (us)",dg));
    gdriftspill[dg]->GetXaxis()->SetTitle("TT1 (s)");
    gdriftspill[dg]->SetTitle(Form("%s spill %i",title.c_str(),info->SpillNumber));
    gdriftspill[dg]->Fit("pol1","Q");
    min = gdriftspill[dg]->GetHistogram()->GetMinimum();
    max = gdriftspill[dg]->GetHistogram()->GetMaximum();
    gdriftspill[dg]->GetYaxis()->SetRangeUser(min,max*(1.4));
    gdriftspill[dg]->SetTitle(Form("%s",title.c_str()));
    gdriftspill[dg]->Fit("pol1","Q");
    gdriftspill[dg]->Draw("p");
    if (saveplots) c->Print(Form("%s_driftspill_%i.png",plotout.c_str(),dg));

    c = new TCanvas();      
    gdriftspillsingle[dg]->Draw("ap");
    gdriftspillsingle[dg]->GetYaxis()->SetTitle(Form("TT%i-TT1 (us)",dg));
    gdriftspillsingle[dg]->GetXaxis()->SetTitle("TT1 (s)");
    gdriftspillsingle[dg]->SetTitle(Form("%s spill 1",title.c_str()));
    gdriftspillsingle[dg]->Fit("pol1","Q");
    min = gdriftspillsingle[dg]->GetHistogram()->GetMinimum();
    max = gdriftspillsingle[dg]->GetHistogram()->GetMaximum();
    gdriftspillsingle[dg]->GetYaxis()->SetRangeUser(min,max*(1.4));
    gdriftspillsingle[dg]->Draw("p");
    if (saveplots) c->Print(Form("%s_driftspillsingle_%i.png",plotout.c_str(),dg));

    // remove drift slope
    for (int p=0; p<gdrift[dg]->GetN(); p++) {
      double diff1,tt1;
      gdrift[dg]->GetPoint(p,tt1,diff1);
      double diff1new = diff1 - hslope[dg]->GetMean()*tt1;
      gdriftcor[dg]->SetPoint(p,tt1,diff1new);
    }

    c = new TCanvas();      
    gdriftcor[dg]->Draw("ap");
    gdriftcor[dg]->GetYaxis()->SetTitle(Form("(TT%i-TT1)-fit (us)",dg));
    gdriftcor[dg]->GetXaxis()->SetTitle("TT1 (min)");
    gdriftcor[dg]->SetTitle(Form("%s corrected",title.c_str()));
    gdriftcor[dg]->Draw("p");
    if (saveplots) c->Print(Form("%s_driftcor_%i.png",plotout.c_str(),dg));

    // get fit values
    TF1 * fit = gdrift[dg]->GetFunction("pol1");
    TF1 * fitspill = gdriftspill[dg]->GetFunction("pol1");
    TF1 * fitspillsingle = gdriftspillsingle[dg]->GetFunction("pol1");

    // add last spill result
    //cout << "filling last spill " << info->SpillNumber << endl;
    slopespill[dg][info->SpillNumber] = fitspill->GetParameter(1);
    hslope[dg]->Fill(fitspill->GetParameter(1));

    //cout << "saved slopes and periods size " << slopespill[dg].size() << endl;
    //for (map<int,double>::iterator sp=slopespill[dg].begin(); sp!=slopespill[dg].end(); sp++) {
    //  cout << "  spill " << sp->first << " slope " << sp->second << endl;
    //}

    // draw drift slope
    c = new TCanvas();               
    hslope[dg]->Draw();             
    if (saveplots) c->Print(Form("%s_slope_%i.png",plotout.c_str(),dg));

    //slope[dg] = slopespill[dg][0];//slope spill 0
    slope[dg] = hslope[dg]->GetMean();//mean slope
    period[dg] = 8/(1+1e-6*slope[dg]);

    cout << "digitizer " << dg << endl;
    cout << " slope run " << fit->GetParameter(1) << endl;
    cout << " slope spill " << fitspill->GetParameter(1) << endl;
    cout << " slope spill single " << fitspillsingle->GetParameter(1) << endl;
    cout << " slope mean " << hslope[dg]->GetMean() << endl;
    cout << " slope dev  " << hslope[dg]->GetStdDev() << endl;
    cout << " period-8ns " << period[dg]-8 << endl;
                            
  }// ndigitizers

  // copy results to a text file
  ofstream fout(slopes_file,ofstream::app);
  fout << run << " ";
  fout << setprecision(20);
  fout << slope[0] << " ";
  fout << slope[2] << " ";
  fout << slopespill[0].size() << " ";
  fout << slopespill[2].size() << endl;
  fout.close();

}
