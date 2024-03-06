#define EventInfo_cxx
#define PMT_cxx

#include "EventInfo.h"
#include "PMT.h"

void EventInfo::Loop() {}
void PMT::Loop() {}

void triggerTimeDrift(string input = "singlePE_-16ns_45ns_run462.root",
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

  // is low momentum or tagged gamma configuration?
  int run_number = stoi(run);
  bool isLM = (run_number<=579) ? true : false;

  // input trees
  vector<string> tree_name;
  if (isLM) {
    tree_name.push_back("ACT0L");//0 - digitizer 0
    tree_name.push_back("ACT0R");
    tree_name.push_back("ACT1L");//2
    tree_name.push_back("ACT1R");
    tree_name.push_back("ACT2L");//4
    tree_name.push_back("ACT2R");
    tree_name.push_back("ACT3L");//6
    tree_name.push_back("ACT3R");
    tree_name.push_back("TOF00");//8 - digitizer 1
    tree_name.push_back("TOF01");
    tree_name.push_back("TOF02");//10
    tree_name.push_back("TOF03");
    tree_name.push_back("TOF10");//12
    tree_name.push_back("TOF11");
    tree_name.push_back("TOF12");//14
    tree_name.push_back("TOF13");
    tree_name.push_back("Hole0");//16 - digitizer 2
    tree_name.push_back("Hole1");
    tree_name.push_back("PbGlass");//18
  }
  else {
    tree_name.push_back("ACT0L");//0 - digitizer 0
    tree_name.push_back("ACT0R");
    tree_name.push_back("ACT1L");//2
    tree_name.push_back("ACT1R");
    tree_name.push_back("ACT3L");//4
    tree_name.push_back("ACT3R");
    tree_name.push_back("TriggerScint");//6
    tree_name.push_back("TOF00");//7 - digitizer 1
    tree_name.push_back("TOF01");//8
    tree_name.push_back("TOF02");
    tree_name.push_back("TOF03");//10
    tree_name.push_back("TOF10");
    tree_name.push_back("TOF11");//12
    tree_name.push_back("TOF12");
    tree_name.push_back("TOF13");//14
    tree_name.push_back("PbGlass");//15 - digitizer 2
    tree_name.push_back("HD8");//16
    tree_name.push_back("HD9");
    tree_name.push_back("HD10");//18
    tree_name.push_back("HD11");
    tree_name.push_back("HD12");//20
    tree_name.push_back("HD13");
    tree_name.push_back("HD14");//22
    tree_name.push_back("HD0");//23 - digitizer 3
    tree_name.push_back("HD1");//24
    tree_name.push_back("HD2");
    tree_name.push_back("HD3");//26
    tree_name.push_back("HD4");
    tree_name.push_back("HD5");//28
    tree_name.push_back("HD6");
    tree_name.push_back("HD7");//30
  }
  tree_name.push_back("EventInfo");
  const int npmts = tree_name.size()-1;

  const int ndigitizers = isLM ? 3 : 4;
  // all pmt channels in a digitizer have
  // the same trigger time and time stamp
  // so pick one to represent each digitizer
  vector<int> dg_pmt;
  if (isLM) {
    dg_pmt.push_back(0);//ACT0L for digitizer 0
    dg_pmt.push_back(8);//TOF00 for digitizer 1
    dg_pmt.push_back(16);//PbGlass for digitizer 2
  }
  else {
    dg_pmt.push_back(0);//ACT0L for digitizer 0
    dg_pmt.push_back(7);//TOF00 for digitizer 1
    dg_pmt.push_back(15);//PbGlass for digitizer 2
    dg_pmt.push_back(23);//HD0 digitizer 3
  }

  unsigned long offset = 2147483648;//2^31 steps = ~17s
  double to_s  = 8e-9;
  int mintrigs = 100;

  // get trees
  EventInfo * info;
  vector<PMT*> pmt;
  vector<TTree*> tree;
  for (int itree=0; itree<tree_name.size(); itree++) {
    TTree * t = (TTree*)f->Get(tree_name[itree].c_str());
    if (t==0) {
      cout << tree_name[itree] << " tree doesn't exist" << endl;
      return;
    }
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
                            100,-2.1,-1.9);
    }
    else if (dg==2) {
      hslope[dg] = new TH1D(Form("hslope%i",dg),
                            Form(";Digitizer %i slope (us/s)",dg),
                            100,0.5,0.7);
    }
    else if (dg==3) {
      hslope[dg] = new TH1D(Form("hslope%i",dg),
                            Form(";Digitizer %i slope (us/s)",dg),
                            100,1.2,1.5);
    }
  }

  // save for each spill
  vector< map<int,double> > slopespill(ndigitizers);

  // trigger times containers for all digitizers
  vector<unsigned long> timeStamp(ndigitizers,0);
  vector<unsigned long> triggerTime(ndigitizers,0);
  vector<unsigned long> triggerTimeFull(ndigitizers,0);
  vector<unsigned long> triggerTimeFullInit(ndigitizers,0);
  vector<double>        timeFullInit(ndigitizers,0);
  vector<unsigned long> triggerTimeFullSpill(ndigitizers,0);
  vector<double>        timeFullSpill(ndigitizers,0);

  unsigned int spillNumber_pre = 0;
  vector<unsigned long> timeStamp_pre(ndigitizers,0);
  vector<unsigned long> timeStamp_init(ndigitizers,0);
  vector<unsigned long> triggerTime_pre(ndigitizers,0);
  vector<unsigned long> triggerTime_full_pre(ndigitizers,0);
  vector<unsigned long> triggerTime_off(ndigitizers,0);
  vector<unsigned long> triggerTime_init(ndigitizers,0);
  vector<unsigned long> triggerTime_spill(ndigitizers,0);

  // loop over entries
  for (int ientry=0; ientry<nentries+1; ientry++) {

    // load all trees
    info->GetEntry(ientry);
    for (int ipmt=0; ipmt<npmts; ipmt++) {
      pmt[ipmt]->GetEntry(ientry);
    }

    // initialize spill number
    // because it doesn't always start at 0
    if (ientry==0) spillNumber_pre = info->SpillNumber;

    // calculate slope for previous spill
    if ((info->SpillNumber!=spillNumber_pre) ||
        ientry==nentries ) {

      for (int d=0; d<ndigitizers; d++) {
        // do fit if there are at least
        // 100 triggers in this spill
        // skip digitizer 1 (TOF)
        if (d==1) {
          slopespill[d][spillNumber_pre] = 0;
        }
        else {
          if (gdriftspill[d]->GetN()>mintrigs) {
            gdriftspill[d]->Fit("pol1","Q");
            TF1 * fitspill = gdriftspill[d]->GetFunction("pol1");
            slopespill[d][spillNumber_pre] = fitspill->GetParameter(1);
            hslope[d]->Fill(fitspill->GetParameter(1));
            delete fitspill;
          }
          else {
            cout << "spill " << spillNumber_pre;
            cout << " digitizer " << d;
            cout << " triggers " << gdriftspill[d]->GetN();
            double min = TMath::MinElement(gdriftspill[d]->GetN(),gdriftspill[d]->GetX());
            double max = TMath::MaxElement(gdriftspill[d]->GetN(),gdriftspill[d]->GetX());
            cout << " range " << max-min;
            cout << endl;
          }
        }

        // reset graph
        delete gdriftspill[d];
        gdriftspill[d] = new TGraph();
      }

      if (ientry==nentries) break;

    }

    // calculate an offset to trigger time to
    // remove effect of trigger time overflow
    // trigger time resets after 2^31
    for (int dg=0; dg<ndigitizers; dg++) {

      // trigger time and time stamp for each digitizer
      triggerTime[dg] = pmt[dg_pmt[dg]] ->triggerTime;
      timeStamp[dg] = pmt[dg_pmt[dg]] ->timeStamp;

      // initilize
      if (ientry==0) {
        cout << "digitizer " << dg;
        cout << " first TT "     << triggerTime[dg];
        cout << " TT (s) " << triggerTime[dg]*to_s;
        cout << endl;
        timeStamp_pre[dg]        = timeStamp[dg];
        timeStamp_init[dg]       = timeStamp[dg];
        triggerTime_pre[dg]      = triggerTime[dg];
        triggerTime_full_pre[dg] = triggerTime[dg];
        triggerTime_init[dg]     = triggerTime[dg];
        triggerTime_spill[dg]    = triggerTime[dg];
      }

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

      // trigger time plus total offset
      triggerTimeFull[dg] = triggerTime[dg] + triggerTime_off[dg];

      // trigger time from run start
      triggerTimeFullInit[dg] = triggerTimeFull[dg] - triggerTime_init[dg];
      timeFullInit[dg] = triggerTimeFullInit[dg] * to_s;

      // trigger time from spill start
      triggerTimeFullSpill[dg] = triggerTimeFull[dg] - triggerTime_spill[dg];
      timeFullSpill[dg] = triggerTimeFullSpill[dg] * to_s;

      // save times for next loop
      timeStamp_pre[dg] = timeStamp[dg];
      triggerTime_pre[dg] = triggerTime[dg];
      triggerTime_full_pre[dg] = triggerTime[dg] + triggerTime_off[dg];

    }//ndigitizer

    // save spill number
    if (info->SpillNumber!=spillNumber_pre) {
      spillNumber_pre = info->SpillNumber;
    }

    // compare with digitizer 1
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
      if (info->SpillNumber==0) {
        gdriftspillsingle[dg]->SetPoint(gdriftspillsingle[dg]->GetN(),timeFullSpill[1],tFullSpillDiff);
      }

    }//digitizers

  }// end input tree

  cout << "last spill number " << info->SpillNumber << endl;
  cout << "saved spills " << slopespill[1].size() << endl;
  cout << "lost spills " << int(info->SpillNumber+1) - slopespill[1].size() << endl;

  for (int dg=0; dg<ndigitizers; dg++) {
    cout << "digitizer " << dg;
    cout << " TT since run start (s) ";
    cout << (triggerTime[dg]+triggerTime_off[dg]-triggerTime_init[dg])*to_s;
    cout << endl;
  }//digitizers

  for (int dg=0; dg<ndigitizers; dg++) {

    if (dg==1) continue;

    // draw plots
    TCanvas * c;
    double min, max;

    // all spills
    TF1 * fit = nullptr;
    if (gdrift[dg]->GetN()>mintrigs) {
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
      fit = gdrift[dg]->GetFunction("pol1");
    }

    // spill 0 fit
    TF1 * fitspillsingle = nullptr;
    if (gdriftspillsingle[dg]->GetN()>mintrigs) {
      c = new TCanvas();
      gdriftspillsingle[dg]->Draw("ap");
      gdriftspillsingle[dg]->GetYaxis()->SetTitle(Form("TT%i-TT1 (us)",dg));
      gdriftspillsingle[dg]->GetXaxis()->SetTitle("TT1 (s)");
      gdriftspillsingle[dg]->SetTitle(Form("%s spill 0",title.c_str()));
      gdriftspillsingle[dg]->Fit("pol1","Q");
      min = gdriftspillsingle[dg]->GetHistogram()->GetMinimum();
      max = gdriftspillsingle[dg]->GetHistogram()->GetMaximum();
      gdriftspillsingle[dg]->GetYaxis()->SetRangeUser(min,max*(1.4));
      gdriftspillsingle[dg]->Draw("p");
      if (saveplots) c->Print(Form("%s_driftspillsingle_%i.png",plotout.c_str(),dg));
      fitspillsingle = gdriftspillsingle[dg]->GetFunction("pol1");
    }

    // remove drift
    for (int p=0; p<gdrift[dg]->GetN(); p++) {
      double diff1,tt1;
      gdrift[dg]->GetPoint(p,tt1,diff1);
      double diff1new = diff1 - hslope[dg]->GetMean()*tt1;
      gdriftcor[dg]->SetPoint(p,tt1,diff1new);
    }

    if (gdriftcor[dg]->GetN()>1) {
      c = new TCanvas();
      gdriftcor[dg]->Draw("ap");
      gdriftcor[dg]->GetYaxis()->SetTitle(Form("(TT%i-TT1)-fit (us)",dg));
      gdriftcor[dg]->GetXaxis()->SetTitle("TT1 (min)");
      gdriftcor[dg]->SetTitle(Form("%s corrected",title.c_str()));
      gdriftcor[dg]->Draw("p");
      if (saveplots) c->Print(Form("%s_driftcor_%i.png",plotout.c_str(),dg));
    }

    //cout << "saved slopes and periods size " << slopespill[dg].size() << endl;
    //for (map<int,double>::iterator sp=slopespill[dg].begin(); sp!=slopespill[dg].end(); sp++) {
    //  cout << "  spill " << sp->first << " slope " << sp->second << endl;
    //}

    // draw drift slope
    c = new TCanvas();               
    hslope[dg]->Draw();             
    if (saveplots) c->Print(Form("%s_slope_%i.png",plotout.c_str(),dg));

    cout << "digitizer " << dg << endl;
    cout << " entries full run " << gdrift[dg]->GetN() << endl;
    cout << " entries spill 0 " << gdriftspillsingle[dg]->GetN() << endl;
    if (fit) cout << " slope full run " << fit->GetParameter(1) << endl;
    if (fitspillsingle) cout << " slope 1 spill " << fitspillsingle->GetParameter(1) << endl;
    cout << " slope mean " << hslope[dg]->GetMean() << endl;
    cout << " slope dev  " << hslope[dg]->GetStdDev() << endl;
    cout << " slope entries " << hslope[dg]->GetEntries() << endl;
  }// ndigitizers


  // copy results to a text file
  ofstream fout(slopes_file,ofstream::app);
  fout << run;
  fout << setprecision(10);
  for (int dg=0; dg<ndigitizers; dg++) {
    if (dg==1) continue;//skip tof
    fout << " " << hslope[dg]->GetMean();
    fout << " " << hslope[dg]->GetStdDev();
    fout << " " << hslope[dg]->GetEntries();
  }
  fout << endl;
  fout.close();

}
