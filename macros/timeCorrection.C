#define EventInfo_cxx
#define PMT_cxx
#include <TTreeReader.h>
#include <TTreeReaderValue.h>


#include "EventInfo.h"
#include "PMT.h"
#include "runs.h"

void EventInfo::Loop() {}
void PMT::Loop() {}

void timeCorrection(string input = "singlePE_-16ns_45ns_run462.root",
                    string slopes_file = "triggerTimeDrift.txt",
                    string output = "test_run462_corrected.root")
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

  //find momentum
  int run_number = stoi(run);
  int mom = runsDict[run_number];
  cout << "momentum " << mom << endl;
  if (mom==0) return;

  // is low momentum or tagged gamma configuration?
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

  // compare these channels times with TOF10
  vector<int> chan;
  if (isLM) {
    for (int i=0;i<8;i++) chan.push_back(i);//ACTs
    chan.push_back(18);//PbGlass
  }
  else {
    for (int i=0;i<6;i++) chan.push_back(i);//ACTs
    chan.push_back(6);//TriggerScint
    chan.push_back(15);//PbGlass
    for (int i=16;i<31;i++) chan.push_back(i);//HD0-14
  }
  const int nchans = chan.size();

  // leadglas channel
  int leadchan;
  if (isLM) leadchan = 18;
  else      leadchan = 15;

  unsigned long offset = 2147483648;//2^31 steps = ~17s
  double to_s  = 8e-9;

  // read slopes from text file
  vector<double> slope(ndigitizers,0);
  ifstream fslopes(slopes_file);
  while (fslopes.good()) {
    string line;
    getline(fslopes,line);
    if (line.size()==0) continue;
    istringstream ss(line);
    string runread;
    ss >> runread;
    if (runread==run) {
      for (int dg=0; dg<ndigitizers; dg++) {
        if (dg==1) continue;//skip tof
        ss >> slope[dg];
        double dev;
        int num;
        ss >> dev;//skip slope stddev
        ss >> num;//skip # triggers
      }
    }
  }

  bool allok = true;
  for (int dg=0; dg<ndigitizers; dg++) {
    if (dg==1) continue;//skip tof
    if (slope[dg]==0) allok = false;
  }

  if (allok) {
    cout << "slopes found" << endl;
    for (int dg=0; dg<ndigitizers; dg++) {
      cout << "digitizer " << dg << " slope " << slope[dg] << endl;
    }
  }
  else {
    cout << "slopes not found for run " << run << endl;
    return;
  }

  vector<double> period(ndigitizers,0);
  for (int dg=0; dg<ndigitizers; dg++) {
    period[dg] = 8/(1+1e-6*slope[dg]);//ns
    cout << "digitizer " << dg << " period-8ns " << period[dg]-8 << endl;
  }

  // ACT thresholds
  vector<double> thresholdv;
  thresholdv.push_back(0.20);//ACT0
  thresholdv.push_back(0.20);
  thresholdv.push_back(0.05);//ACT1
  thresholdv.push_back(0.05);
  thresholdv.push_back(0.10);//ACT2
  thresholdv.push_back(0.05);
  thresholdv.push_back(0.10);//ACT3
  thresholdv.push_back(0.05);

  vector<double> thresholdq;
  thresholdq.push_back(0.20);//ACT0
  thresholdq.push_back(0.20);
  thresholdq.push_back(0.05);//ACT1
  thresholdq.push_back(0.05);
  thresholdq.push_back(0.10);//ACT2
  thresholdq.push_back(0.05);
  thresholdq.push_back(0.10);//ACT3
  thresholdq.push_back(0.05);

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

  // voltage and charge histos
  TH1D * hv[nchans];
  TH1D * hq[nchans];
  for (int chi=0; chi<nchans; chi++) {
    hv[chi] = new TH1D(Form("%s_V",tree_name[chan[chi]].c_str()),"",100,0,1);
    hq[chi] = new TH1D(Form("%s_Q",tree_name[chan[chi]].c_str()),"",100,0,0.25);
    hv[chi]->SetTitle(Form("Run %d %d MeV/c;V",run_number,mom));
    hq[chi]->SetTitle(Form("Run %d %d MeV/c;Q",run_number,mom));
  }

  // tof vs lg
  TH2D * htoflg = new TH2D("toflg","",100,10,15,100,0,0.5);
  htoflg->SetTitle(Form("Run %d %d MeV/c",run_number,mom));
  htoflg->SetXTitle("TOF (ns)");
  htoflg->SetYTitle("Lead Glass Charge");

  // tof
  TH1D * htof0 = new TH1D("tof0","",100,0,200);
  TH1D * htof1 = new TH1D("tof1","",100,0,200);
  TH1D * htof00 = new TH1D("tof00","",100,0,200);
  TH1D * htof10 = new TH1D("tof10","",100,0,200);

  htof0->SetTitle(Form("Run %d %d MeV/c;Signal time (ns)",run_number,mom));
  htof1->SetTitle(Form("Run %d %d MeV/c;Signal time (ns)",run_number,mom));
  htof00->SetTitle(Form("Run %d %d MeV/c;Signal time (ns)",run_number,mom));
  htof10->SetTitle(Form("Run %d %d MeV/c;Signal time (ns)",run_number,mom));

  // time difference
  TH1D * htdiff[nchans];
  TH1D * htdifftmp[nchans];
  TH1D * htdiffcor[nchans];
  TH1D * hoffmean [nchans];
  TH1D * hoffstd[nchans];
  TH2D * hoffcorre[nchans];

  for (int chi=0;chi<nchans;chi++) {

    htdiff[chi] = new TH1D(Form("%s_TOF10",tree_name[chan[chi]].c_str()),"",100,-100,100);
    htdiff[chi]->SetXTitle(Form("%s-TOF10 (ns)",tree_name[chan[chi]].c_str()));
    htdiff[chi]->SetTitle("e-like");

    htdifftmp[chi] = new TH1D(Form("%s_TOF10_tmp",tree_name[chan[chi]].c_str()),"",100,-100,100);
    htdifftmp[chi]->SetXTitle(Form("%s-TOF10 (ns)",tree_name[chan[chi]].c_str()));
    htdifftmp[chi]->SetTitle("e-like");

    htdiffcor[chi] = new TH1D(Form("%s_TOF10_cor",tree_name[chan[chi]].c_str()),"",100,-100,200);
    htdiffcor[chi]->SetXTitle(Form("%s-TOF10 (ns)",tree_name[chan[chi]].c_str()));
    htdiffcor[chi]->SetTitle("e-like");

    hoffmean[chi] = new TH1D(Form("%soffmean",tree_name[chi].c_str()),
                           Form(";%s-TOF10 mean per spill (ns)",tree_name[chi].c_str()),
                           100,-100,100);

    hoffstd[chi] = new TH1D(Form("%soffstd",tree_name[chan[chi]].c_str()),
                                 Form(";%s-TOF10 stddev per spill (ns)",tree_name[chan[chi]].c_str()),
                                 100,0,20);

    hoffcorre[chi] = new TH2D(Form("%soffcorre",tree_name[chan[chi]].c_str()),
                                 Form(";%s-TOF01 (ns);ACT3R-TOF10",tree_name[chan[chi]].c_str()),
                                 100,-100,100,100,-100,100);
  }

  // act and lg trigger time correction

  vector<TH1D*> httcor(ndigitizers);
  for (int dg=0; dg<ndigitizers; dg++) {
    httcor[dg] = new TH1D(Form("ttcor%i",dg),Form(";TT%i-TT1 (ns)",dg),100,-100,100);
  }

  // save for each spill
  vector< map<int,double> > off_mean(nchans);
  vector< map<int,double> > off_std(nchans);

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

    // calculate time offset for previous spill
    if ((info->SpillNumber!=spillNumber_pre) ||
        ientry==nentries ) {

      for (int chi=0;chi<nchans;chi++) {
        hoffmean[chi]->Fill(htdifftmp[chi]->GetMean());
        hoffstd[chi]->Fill(htdifftmp[chi]->GetStdDev());
        off_mean[chi][spillNumber_pre] = htdifftmp[chi]->GetMean();
        off_std[chi][spillNumber_pre] = htdifftmp[chi]->GetStdDev();
        if (htdifftmp[chi]->GetBinContent(0)>0 ||
            htdifftmp[chi]->GetBinContent(htdifftmp[chi]->GetNbinsX())>0) {
          cout << "spill number " << spillNumber_pre << endl;
          cout << tree_name[chi] << " ";
          cout << htdifftmp[chi]->GetMean() << " ";
          cout << htdifftmp[chi]->GetStdDev() << " ";
          cout << htdifftmp[chi]->GetBinContent(0) << " ";
          cout << htdifftmp[chi]->GetBinContent(htdifftmp[chi]->GetNbinsX());
          cout << endl;
        }
        htdifftmp[chi]->Reset();
      }

      for (int chi=0;chi<8;chi++) {
        hoffcorre[chi]->Fill(off_mean[chi][spillNumber_pre],off_mean[7][spillNumber_pre]);
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

    // trigger times
    vector<double> tt(ndigitizers,0);
    for (int dg=0; dg<ndigitizers; dg++) {
      tt[dg] = triggerTimeFullSpill[dg]*period[dg];
    }

    // trigger time correction for each digitizer
    vector<double> ttcor(ndigitizers,0);
    for (int dg=0; dg<ndigitizers; dg++) {
      ttcor[dg] = tt[dg]-tt[1];
      httcor[dg]->Fill(ttcor[dg]);
    }

    // one peak in TOF and LeadGlass
    bool cut = true;
    for (int ch=8; ch<16; ch++) cut = cut && (pmt[ch]->nPeaks==1);
    cut = cut && pmt[leadchan]->nPeaks==1;

    // TOF and LeadGlass signal time in first bunch
    for (int ch=8; ch<16; ch++) cut = cut && (pmt[ch]->SignalTime[0]<200);
    cut = cut && pmt[leadchan]->SignalTime[0]<200;

    if (cut==true) {

      double tof0 = (pmt[8] ->SignalTime[0]+
                     pmt[9] ->SignalTime[0]+
                     pmt[10]->SignalTime[0]+
                     pmt[11]->SignalTime[0])/4.;
      double tof1 = (pmt[12]->SignalTime[0]+
                     pmt[13]->SignalTime[0]+
                     pmt[14]->SignalTime[0]+
                     pmt[15]->SignalTime[0])/4.;
      double tof = tof1-tof0;
      double lg  = pmt[leadchan]->IntCharge[0];
      double tref = pmt[12]->SignalTime[0];

      htoflg->Fill(tof,lg);
      htof0->Fill(tof0);
      htof00->Fill(pmt[8]->SignalTime[0]);
      htof1->Fill(tof1);
      htof10->Fill(pmt[12]->SignalTime[0]);

      // e-like selection
      bool is_electron = false;
      if (abs(mom)<400) {
        is_electron = tof<12;
      }
      else if(abs(mom)<700) {
        is_electron = tof<12 && lg>0.2;
      }
      else {
        is_electron = tof<12 && lg>0.3;
      }

      if (is_electron) {

        // pulse charge and voltage histos
        for (int chi=0; chi<nchans; chi++) {
          hv[chi]->Fill(pmt[chan[chi]]->PeakVoltage[0]);
          hq[chi]->Fill(pmt[chan[chi]]->IntCharge[0]);
        }

        // act digitizer
        for (int acti=0; acti<8; acti++) {
          if (pmt[acti]->nPeaks==1 &&
              pmt[acti]->SignalTime[0]<200 &&
              pmt[acti]->PeakVoltage[0]>thresholdv[acti]) {
            double actdiff = pmt[acti]->SignalTime[0]-tref;
            htdiff[acti]->Fill(actdiff);
            htdifftmp[acti]->Fill(actdiff+ttcor[0]);
          }
        }

        // lead glass digitizer
        double lgdiff = pmt[leadchan]->SignalTime[0]-tref;
        htdiff[8]->Fill(lgdiff);
        htdifftmp[8]->Fill(lgdiff+ttcor[2]);

      }
    }

  }// end input tree

  cout << "last spill number " << info->SpillNumber << endl;
  for (int dg=0; dg<ndigitizers; dg++) {
    cout << "digitizer " << dg;
    cout << " TT since run start (s) ";
    cout << (triggerTime[dg]+triggerTime_off[dg]-triggerTime_init[dg])*to_s;
    cout << endl;
  }//digitizers

  if (saveplots) {

    TCanvas * c;

    c = new TCanvas();
    htoflg->SetStats(0);
    htoflg->Draw("colz");
    c->Print(Form("%s_tof_lg.png",plotout.c_str()));

    c = new TCanvas();
    htof0->SetStats(0);
    htof0->SetLineColor(9);
    htof0->Draw();
    htof00->SetLineColor(46);
    htof00->Draw("same");
    TLegend * leg = new TLegend(0.55,0.7,0.85,0.85);
    leg->AddEntry(htof0,Form("TOF0 StdDev=%.2f ns",htof0->GetStdDev()),"l");
    leg->AddEntry(htof00,Form("TOF00 StdDev=%.2f ns",htof00->GetStdDev()),"l");
    leg->Draw();
    c->Print(Form("%s_tof0.png",plotout.c_str()));

    c = new TCanvas();
    htof1->SetStats(0);
    htof1->SetLineColor(9);
    htof1->Draw();
    htof10->SetLineColor(46);
    htof10->Draw("same");
    leg = new TLegend(0.55,0.7,0.85,0.85);
    leg->AddEntry(htof1,Form("TOF1 StdDev=%.2f ns",htof1->GetStdDev()),"l");
    leg->AddEntry(htof10,Form("TOF10 StdDev=%.2f ns",htof10->GetStdDev()),"l");
    leg->Draw();
    c->Print(Form("%s_tof1.png",plotout.c_str()));

    for (int chi=0; chi<nchans; chi++) {
      c = new TCanvas();
      hv[chi]->SetStats(1);
      hv[chi]->Draw();
      c->Print(Form("%s_%s_V.png",plotout.c_str(),tree_name[chan[chi]].c_str()));
    
      c = new TCanvas();
      hq[chi]->SetStats(1);
      hq[chi]->Draw();
      c->Print(Form("%s_%s_Q.png",plotout.c_str(),tree_name[chan[chi]].c_str()));
    }

    for (int dg=0; dg<ndigitizers; dg++) {
      c = new TCanvas();
      httcor[dg]->Draw();
      c->Print(Form("%s_ttcor%i.png",plotout.c_str(),dg));
    }

    for (int chi=0;chi<nchans;chi++) {
      c = new TCanvas();
      hoffmean[chi]->Draw();
      c->Print(Form("%s_%s_T10_offmean.png",plotout.c_str(),tree_name[chi].c_str()));

      c = new TCanvas();
      hoffstd[chi]->Draw();
      c->Print(Form("%s_%s_T10_offstd.png",plotout.c_str(),tree_name[chi].c_str()));

      c = new TCanvas();
      hoffcorre[chi]->Draw("colz");
      c->Print(Form("%s_%s_offcorre.png",plotout.c_str(),tree_name[chi].c_str()));
    }

  }

  // mean and std corrections
  //for (int i=0; i<9; i++) {
  //  cout << "channel " << i << endl;
  //  cout << "saved spill corrections size " << off_mean[i].size() << endl;
  //  for (map<int,double>::iterator sp=off_mean[i].begin(); sp!=off_mean[i].end(); sp++) {
  //    cout << "  spill " << sp->first << " " << sp->second << endl;
  //  }
  //}

  // clone input tree and replace Signal time with
  // new values
  // save cloned tree in a different file
  TFile * newfile = new TFile(output.c_str(),"recreate");
  int nPeaksMax = 100;
  double signalTimeCor[npmts][nPeaksMax];
  for (int ipmt=0; ipmt<npmts; ipmt++) {
    for (int i=0;i<nPeaksMax;i++){
      //setting nans instead of 0s to keep the correct structure and only fill the relevant peaks
      signalTimeCor[ipmt][i] = std::numeric_limits<double>::quiet_NaN();
    }
  }


  // TTreeReader reader(inputTree);


  double signalTimeCorVal = 0;

  vector<TTree*> newtree;
  for (int itree=0; itree<tree.size(); itree++) {
    newfile->cd();
    TTree * t;
    if (itree<npmts) {
      t = tree[itree]->CloneTree(0);
      //Modified the new signalTimeCorrected branch to actually hold the corrected timings
      //previously it contains the base signal Times
      t->SetBranchAddress("SignalTimeCorrected",&signalTimeCor[itree]);

    }
    else {
      t = tree[itree]->CloneTree(0);
    }
    newtree.push_back(t);
  }

  // reset trigger times containers
  timeStamp.assign(ndigitizers,0);
  triggerTime.assign(ndigitizers,0);
  triggerTimeFull.assign(ndigitizers,0);
  triggerTimeFullInit.assign(ndigitizers,0);
  timeFullInit.assign(ndigitizers,0);
  triggerTimeFullSpill.assign(ndigitizers,0);
  timeFullSpill.assign(ndigitizers,0);

  spillNumber_pre = 0;
  timeStamp_pre.assign(ndigitizers,0);
  timeStamp_init.assign(ndigitizers,0);
  triggerTime_pre.assign(ndigitizers,0);
  triggerTime_full_pre.assign(ndigitizers,0);
  triggerTime_off.assign(ndigitizers,0);
  triggerTime_init.assign(ndigitizers,0);
  triggerTime_spill.assign(ndigitizers,0);

  // loop over entries again
  for (int ientry=0; ientry<nentries; ientry++) {

    // load all trees
    info->GetEntry(ientry);
    for (int ipmt=0; ipmt<npmts; ipmt++) {
      pmt[ipmt]->GetEntry(ientry);

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
        spillNumber_pre          = info->SpillNumber;
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

    }//digitizers

    // save current spill number
    if (info->SpillNumber!=spillNumber_pre) {
      spillNumber_pre = info->SpillNumber;
    }

    // trigger times
    vector<double> tt(ndigitizers,0);
    for (int dg=0; dg<ndigitizers; dg++) {
      tt[dg] = triggerTimeFullSpill[dg]*period[dg];
    }

    // trigger time correction for each digitizer
    vector<double> ttcor(ndigitizers,0);
    for (int dg=0; dg<ndigitizers; dg++) {
      ttcor[dg] = tt[dg]-tt[1];
      httcor[dg]->Fill(ttcor[dg]);
    }

    // full correction of signal time for all channels
    for (int ipmt=0; ipmt<npmts; ipmt++) {
      for (int ipeak=0; ipeak<pmt[ipmt]->nPeaks; ipeak++) {
        // ACT channels
        if (ipmt<8) {
          signalTimeCor[ipmt][ipeak] = pmt[ipmt]->SignalTime[ipeak] + ttcor[0] - off_mean[6][info->SpillNumber];
        }
        // TOF channels
        else if (ipmt<16) {
          signalTimeCor[ipmt][ipeak] = pmt[ipmt]->SignalTime[ipeak];
        }
        // Hole and LG channels
        else {
          signalTimeCor[ipmt][ipeak] = pmt[ipmt]->SignalTime[ipeak] + ttcor[2] - off_mean[8][info->SpillNumber];
        }
        // signalTimeCorVal = signalTimeCor[ipmt][ipeak];

      }
      // std::cout << signalTimeCor[ipmt][0] << " " <<   signalTimeCor[ipmt][1] << " " <<  signalTimeCor[ipmt][2] << std::endl;

      newtree[ipmt]->Fill();
    }

    // e-like selection
    // with full corrections
    // one peak cut
    bool cut = true;
    for (int ch=8; ch<16; ch++) cut = cut && (pmt[ch]->nPeaks==1);
    cut = cut && pmt[leadchan]->nPeaks==1;

    // signal time in first bunch  cut
    for (int ch=8; ch<16; ch++) cut = cut && (pmt[ch]->SignalTime[0]<200);
    cut = cut && pmt[leadchan]->SignalTime[0]<200;

    if (cut==true) {

      // TOF and LG
      double tof0 = (pmt[8] ->SignalTime[0]+
                     pmt[9] ->SignalTime[0]+
                     pmt[10]->SignalTime[0]+
                     pmt[11]->SignalTime[0])/4.;
      double tof1 = (pmt[12]->SignalTime[0]+
                     pmt[13]->SignalTime[0]+
                     pmt[14]->SignalTime[0]+
                     pmt[15]->SignalTime[0])/4.;
      double tof = tof1-tof0;
      double lg  = pmt[leadchan]->IntCharge[0];

      // e-like selection
      if (tof<12 && lg>0.2) {

        for (int acti=0; acti<8; acti++) {
          if (pmt[acti]->nPeaks==1 &&
              pmt[acti]->SignalTime[0]<200 &&
              pmt[acti]->PeakVoltage[0]>thresholdv[acti]) {
            htdiffcor[acti]->Fill(signalTimeCor[acti][0]-signalTimeCor[12][0]);
          }
        }

        htdiffcor[8]->Fill(signalTimeCor[leadchan][0]-signalTimeCor[12][0]);

      }

    }

  }// end input tree

  cout << "last spill number " << info->SpillNumber << endl;
  for (int dg=0; dg<ndigitizers; dg++) {
    cout << "digitizer " << dg;
    cout << " TT since run start (s) ";
    cout << (triggerTime[dg]+triggerTime_off[dg]-triggerTime_init[dg])*to_s;
    cout << endl;
  }//digitizers

  gStyle->SetOptTitle(0);

  if (saveplots) {
    for (int chi=0; chi<nchans; chi++) {
      TCanvas * c = new TCanvas();
      THStack * s = new THStack();
      htdiff[chi]->SetStats(0);
      htdiffcor[chi]->SetStats(0);
      htdiff[chi]->SetLineColor(kBlue);
      htdiffcor[chi]->SetLineColor(kRed);
      s->Add(htdiff[chi]);
      s->Add(htdiffcor[chi]);
      htdiff[chi]->SetTitle("e-like");
      htdiffcor[chi]->SetTitle("e-like corr.");
      s->Draw("nostack");
      s->GetXaxis()->SetTitle(htdiff[chi]->GetXaxis()->GetTitle());
      s->GetYaxis()->SetTitle(htdiff[chi]->GetYaxis()->GetTitle());
      gPad->BuildLegend(0.6,0.7,0.88,0.89);
      //gPad->SetLogy();          
      c->Print(Form("%s_%s_T10_cor_stack.png",plotout.c_str(),tree_name[chan[chi]].c_str()));
      //c->Print(Form("%s_%s_T10_cor_stack_logy.png",plotout.c_str(),tnames[chi].c_str()));
    }
  }

  // copy results to a text file
  ofstream fout("timeCorrection.txt",ofstream::app);
  fout << run;
  for (int chi=0;chi<nchans;chi++) {
    fout << " " << htdiff[chi]->GetStdDev();
    fout << " " << htdiffcor[chi]->GetStdDev();
  }
  fout << endl;
  fout.close();

  // save new ntuple
  newfile->Write();

}
