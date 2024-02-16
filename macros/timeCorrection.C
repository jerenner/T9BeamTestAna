#define EventInfo_cxx
#define PMT_cxx

#include "EventInfo.h"
#include "PMT.h"

void EventInfo::Loop() {}
void PMT::Loop() {}

void timeCorrection(string input = "/neut/datasrv2a/jrenner/ntuple_files/ntuple_000435.root",
                    string slopes_file = "triggerTimeDrift.txt",
                    string output = "ntuple_000435.root")
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
  unsigned long offset = 2147483648;//2^31
  double to_s  = 8e-9;

  // read slopes from text file
  double slope[ndigitizers] = {0,0,0};
  ifstream fslopes(slopes_file);
  while (fslopes.good()) {
    string line;
    getline(fslopes,line);
    if (line.size()==0) continue;
    istringstream ss(line);
    string runread;
    ss >> runread;
    if (runread==run) {
      ss >> slope[0];
      ss >> slope[2];
      break;
    }
  }

  if (slope[0]!=0 || slope[2]!=0) {
    cout << "slopes found" << endl;
    cout << "digitizer 0 slope " << slope[0] << endl;
    cout << "digitizer 2 slope " << slope[2] << endl;
  }
  else {
    cout << "no slopes found" << endl;
    return;
  }

  double period[ndigitizers];
  period[1] = 8;//ns
  period[0] = 8/(1+1e-6*slope[0]);//ns
  period[2] = 8/(1+1e-6*slope[2]);//ns
  cout << "digitizer 0 period-8ns " << period[0]-8 << endl;
  cout << "digitizer 2 period-8ns " << period[2]-8 << endl;

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
  TH1D * hv[npmts];
  TH1D * hq[npmts];
  for (int ipmt=0; ipmt<npmts; ipmt++) {
    hv[ipmt] = new TH1D(Form("%s_V",tree_name[ipmt].c_str()),";V",100,0,0.5);
    hq[ipmt] = new TH1D(Form("%s_Q",tree_name[ipmt].c_str()),";Q",100,0,0.3);
  }

  // tof vs lg
  TH2D * htoflg = new TH2D("toflg","",100,10,15,100,0,0.5);
  htoflg->SetXTitle("TOF (ns)");
  htoflg->SetYTitle("Lead Glass Charge");

  // lg vs act23
  TH2D * hlgact23 = new TH2D("lgact23","",100,0,0.5,100,0,1);
  hlgact23->SetXTitle("Lead Glass Charge");
  hlgact23->SetYTitle("ACT23 Charge");

  // compare these channels times with TOF0
  vector<int> chan;
  chan.push_back(0);
  chan.push_back(1);
  chan.push_back(2);
  chan.push_back(3);
  chan.push_back(4);
  chan.push_back(5);
  chan.push_back(6);
  chan.push_back(7);
  chan.push_back(18);
  const int nchans = chan.size();

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
  TH1D * httcoract = new TH1D("ttcoract",
                              ";TT0-TT1 (ns)",
                              100,-100,100);
  TH1D * httcorlg  = new TH1D("ttcorlg",
                              ";TT2-TT1 (ns)",
                              100,-100,100);

  // save for each spill
  vector< map<int,double> > off_mean(nchans);
  vector< map<int,double> > off_std(nchans);

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

    // calculate time offset for previous spill
    if (info->EventNumber==0 && info->SpillNumber>0) {

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
    }

    // trigger time and time stamp 
    // for each digitizer
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

    // save current spill number
    if (info->SpillNumber!=spillNumber_pre) {
      spillNumber_pre = info->SpillNumber;
    }

    unsigned long triggerTimeFull[3];
    triggerTimeFull[0] = triggerTime[0] + triggerTime_off[0];
    triggerTimeFull[1] = triggerTime[1] + triggerTime_off[1];
    triggerTimeFull[2] = triggerTime[2] + triggerTime_off[2];

    unsigned long triggerTimeFullInit[3];
    triggerTimeFullInit[0] = triggerTimeFull[0] - triggerTime_init[0];
    triggerTimeFullInit[1] = triggerTimeFull[1] - triggerTime_init[1];
    triggerTimeFullInit[2] = triggerTimeFull[2] - triggerTime_init[2];

    double timeFullInit[3];
    timeFullInit[0] = triggerTimeFullInit[0] * to_s;
    timeFullInit[1] = triggerTimeFullInit[1] * to_s;
    timeFullInit[2] = triggerTimeFullInit[2] * to_s;

    unsigned long triggerTimeFullSpill[3];
    triggerTimeFullSpill[0] = triggerTimeFull[0] - triggerTime_spill[0];
    triggerTimeFullSpill[1] = triggerTimeFull[1] - triggerTime_spill[1];
    triggerTimeFullSpill[2] = triggerTimeFull[2] - triggerTime_spill[2];

    // time correction between tof and act digitizers
    double tt_act = triggerTimeFullSpill[0]*period[0];
    double tt_tof = triggerTimeFullSpill[1]*period[1];
    double tt_lg  = triggerTimeFullSpill[2]*period[2];

    // trigger time correction
    double ttcoract = tt_act-tt_tof;
    double ttcorlg  = tt_lg -tt_tof;

    httcoract->Fill(ttcoract);
    httcorlg ->Fill(ttcorlg);

    // e-like selection
    // one peak cut
    bool cut = true;
    for (int ch=8; ch<16; ch++) cut = cut && (pmt[ch]->nPeaks==1);
    cut = cut && pmt[18]->nPeaks==1;

    // signal time in first bunch  cut
    for (int ch=8; ch<16; ch++) cut = cut && (pmt[ch]->SignalTime[0]<200);
    cut = cut && pmt[18]->SignalTime[0]<200;

    // pulse charge and voltage histos
    if (cut==true) {
      for (int chi=0; chi<19; chi++) {
        hv[chi]->Fill(pmt[chi]->PeakVoltage[0]);
        hq[chi]->Fill(pmt[chi]->IntCharge[0]);
      }
    }

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
      double lg  = pmt[18]->IntCharge[0];
      htoflg->Fill(tof,lg);

      // ACT23 charge
      double act2  = pmt[4]->IntCharge[0]+pmt[5]->IntCharge[0];
      double act3  = pmt[6]->IntCharge[0]+pmt[7]->IntCharge[0];
      double act23 = act2+act3;
      hlgact23->Fill(lg,act23);

      if (tof<12 && lg>0.2) {

        // act digitizer
        for (int acti=0; acti<8; acti++) {
          if (pmt[acti]->nPeaks==1 &&
              pmt[acti]->SignalTime[0]<200 &&
              pmt[acti]->PeakVoltage[0]>thresholdv[acti]) {
              //pmt[acti]->IntCharge[0]>thresholdq[acti]) {
            double actdiff = pmt[acti]->SignalTime[0]-pmt[12]->SignalTime[0];
            htdiff[acti]->Fill(actdiff);
            htdifftmp[acti]->Fill(actdiff+ttcoract);
          }
        }

        // lead glass digitizer
        double lgdiff = pmt[18]->SignalTime[0]-pmt[12]->SignalTime[0];
        htdiff[8]->Fill(lgdiff);
        htdifftmp[8]->Fill(lgdiff+ttcorlg);

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

  for (int chi=0;chi<nchans;chi++) {
    hoffmean[chi]->Fill(htdifftmp[chi]->GetMean());
    hoffstd[chi]->Fill(htdifftmp[chi]->GetStdDev());
    off_mean[chi][info->SpillNumber] = htdifftmp[chi]->GetMean();
    off_std[chi][info->SpillNumber] = htdifftmp[chi]->GetStdDev();
    cout << tree_name[chi] << " ";
    cout << "spill " << spillNumber_pre;
    cout << " mean (ns) " << htdifftmp[chi]->GetMean();
    cout << " std (ns) " <<  htdifftmp[chi]->GetStdDev();
    cout << " entries " <<  htdifftmp[chi]->GetEntries() << endl;
    cout << tree_name[chi] << " ";
    cout << "spill " << spillNumber_pre;
    cout << " mean (ns) " << htdiff[chi]->GetMean();
    cout << " std (ns) " <<  htdiff[chi]->GetStdDev();
    cout << " entries " <<  htdiff[chi]->GetEntries() << endl;
    if (htdifftmp[chi]->GetBinContent(0)>0 ||
        htdifftmp[chi]->GetBinContent(htdifftmp[chi]->GetNbinsX())>0) {
      cout << "spill number " << info->SpillNumber << endl;
      cout << tree_name[chi] << " ";
      cout << htdifftmp[chi]->GetMean() << " ";
      cout << htdifftmp[chi]->GetStdDev() << " ";
      cout << htdifftmp[chi]->GetBinContent(0) << " ";
      cout << htdifftmp[chi]->GetBinContent(htdifftmp[chi]->GetNbinsX());
      cout << endl;
    }
    htdifftmp[chi]->Reset();
  }

  if (saveplots) {

    TCanvas * c;

    c = new TCanvas();
    htoflg->SetStats(0);
    htoflg->Draw("colz");
    c->Print(Form("%s_tof_lg.png",plotout.c_str()));

    c = new TCanvas();
    hlgact23->SetStats(0);
    hlgact23->Draw("colz");
    c->Print(Form("%s_lg_act23.png",plotout.c_str()));

    for (int chi=0; chi<19; chi++) {
      c = new TCanvas();
      hv[chi]->SetStats(1);
      hv[chi]->Draw();
      c->Print(Form("%s_%s_V.png",plotout.c_str(),tree_name[chi].c_str()));
    
      c = new TCanvas();
      hq[chi]->SetStats(1);
      hq[chi]->Draw();
      c->Print(Form("%s_%s_Q.png",plotout.c_str(),tree_name[chi].c_str()));
    }

    c = new TCanvas();
    httcoract->Draw();
    c->Print(Form("%s_ttcoract.png",plotout.c_str()));

    c = new TCanvas();
    httcorlg->Draw();
    c->Print(Form("%s_ttcorlg.png",plotout.c_str()));

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

  double signalTimeCor[npmts][100];
  for (int ipmt=0; ipmt<npmts; ipmt++) {
    for (int i=0;i<100;i++){
      signalTimeCor[ipmt][i] = 0;
    }
  }

  vector<TTree*> newtree;
  for (int itree=0; itree<tree.size(); itree++) {
    newfile->cd();
    TTree * t;
    if (itree<npmts) {
      t = tree[itree]->CloneTree(0);
      t->SetBranchAddress("SignalTime",&signalTimeCor[itree]);
    }
    else {
      t = tree[itree]->CloneTree();
    }
    newtree.push_back(t);
  }

  // loop over entries again
  spillNumber_pre = 0;
  timeStamp.assign(ndigitizers,0);
  timeStamp_pre.assign(ndigitizers,0);
  timeStamp_init.assign(ndigitizers,0);
  triggerTime.assign(ndigitizers,0);
  triggerTime_pre.assign(ndigitizers,0);
  triggerTime_full_pre.assign(ndigitizers,0);
  triggerTime_off.assign(ndigitizers,0);
  triggerTime_init.assign(ndigitizers,0);
  triggerTime_spill.assign(ndigitizers,0);

  for (int ientry=0; ientry<nentries; ientry++) {

    // load all trees
    info->GetEntry(ientry);
    for (int ipmt=0; ipmt<npmts; ipmt++) {
      pmt[ipmt]->GetEntry(ientry);
    }

    // trigger time and time stamp 
    // for each digitizer
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

    }//digitizers

    // save current spill number
    if (info->SpillNumber!=spillNumber_pre) {
      spillNumber_pre = info->SpillNumber;
    }

    unsigned long triggerTimeFull[3];
    triggerTimeFull[0] = triggerTime[0] + triggerTime_off[0];
    triggerTimeFull[1] = triggerTime[1] + triggerTime_off[1];
    triggerTimeFull[2] = triggerTime[2] + triggerTime_off[2];

    unsigned long triggerTimeFullInit[3];
    triggerTimeFullInit[0] = triggerTimeFull[0] - triggerTime_init[0];
    triggerTimeFullInit[1] = triggerTimeFull[1] - triggerTime_init[1];
    triggerTimeFullInit[2] = triggerTimeFull[2] - triggerTime_init[2];

    double timeFullInit[3];
    timeFullInit[0] = triggerTimeFullInit[0] * to_s;
    timeFullInit[1] = triggerTimeFullInit[1] * to_s;
    timeFullInit[2] = triggerTimeFullInit[2] * to_s;

    unsigned long triggerTimeFullSpill[3];
    triggerTimeFullSpill[0] = triggerTimeFull[0] - triggerTime_spill[0];
    triggerTimeFullSpill[1] = triggerTimeFull[1] - triggerTime_spill[1];
    triggerTimeFullSpill[2] = triggerTimeFull[2] - triggerTime_spill[2];

    double tt_act = triggerTimeFullSpill[0]*period[0];
    double tt_tof = triggerTimeFullSpill[1]*period[1];
    double tt_lg  = triggerTimeFullSpill[2]*period[2];

    // trigger time correction
    double ttcoract = tt_act-tt_tof;
    double ttcorlg  = tt_lg -tt_tof;

    // full correction of signal time for all channels
    for (int ipmt=0; ipmt<npmts; ipmt++) {
      for (int ipeak=0; ipeak<pmt[ipmt]->nPeaks; ipeak++) {
        // ACT channels
        if (ipmt<8) {
          signalTimeCor[ipmt][ipeak] = pmt[ipmt]->SignalTime[ipeak] + ttcoract - off_mean[6][info->SpillNumber];
        }
        // TOF channels
        else if (ipmt<16) {
          signalTimeCor[ipmt][ipeak] = pmt[ipmt]->SignalTime[ipeak];
        }
        // Hole and LG channels
        else {
          signalTimeCor[ipmt][ipeak] = pmt[ipmt]->SignalTime[ipeak] + ttcorlg - off_mean[8][info->SpillNumber];
        }
      }
      newtree[ipmt]->Fill();
    }

    // e-like selection
    // with full corrections
    // one peak cut
    bool cut = true;
    for (int ch=8; ch<16; ch++) cut = cut && (pmt[ch]->nPeaks==1);
    cut = cut && pmt[18]->nPeaks==1;

    // signal time in first bunch  cut
    for (int ch=8; ch<16; ch++) cut = cut && (pmt[ch]->SignalTime[0]<200);
    cut = cut && pmt[18]->SignalTime[0]<200;

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
      double lg  = pmt[18]->IntCharge[0];

      // e-like selection
      if (tof<12 && lg>0.2) {

        for (int acti=0; acti<8; acti++) {
          if (pmt[acti]->nPeaks==1 &&
              pmt[acti]->SignalTime[0]<200 &&
              pmt[acti]->PeakVoltage[0]>thresholdv[acti]) {
            htdiffcor[acti]->Fill(signalTimeCor[acti][0]-signalTimeCor[12][0]);
          }
        }

        htdiffcor[8]->Fill(signalTimeCor[18][0]-signalTimeCor[12][0]);

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
