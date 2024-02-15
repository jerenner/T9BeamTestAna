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

  TH1D * hactt10[nchans];
  TH1D * hactt10cor[nchans];
  TH1D * hactt10cortmp[nchans];
  TH1D * hspillcoract [nchans];
  TH1D * hspillstdact [nchans];
  TH2D * hactcor[nchans];

  for (int chi=0;chi<nchans;chi++) {

    hactt10[chi] = new TH1D(Form("%s_TOF10",tree_name[chan[chi]].c_str()),"",100,-100,100);
    hactt10[chi]->SetXTitle(Form("%s-TOF10 (ns)",tree_name[chan[chi]].c_str()));
    hactt10[chi]->SetTitle("e-like");

    hactt10cor[chi] = new TH1D(Form("%s_TOF10_cor",tree_name[chan[chi]].c_str()),"",100,-100,200);
    hactt10cor[chi]->SetXTitle(Form("%s-TOF10 (ns)",tree_name[chan[chi]].c_str()));
    hactt10cor[chi]->SetTitle("e-like");

    hactt10cortmp[chi] = new TH1D(Form("%s_TOF10_cor_tmp",tree_name[chan[chi]].c_str()),"",100,-100,100);
    hactt10cortmp[chi]->SetXTitle(Form("%s-TOF10 (ns)",tree_name[chan[chi]].c_str()));
    hactt10cortmp[chi]->SetTitle("e-like");

    hspillcoract[chi] = new TH1D(Form("%sspillcor",tree_name[chi].c_str()),
                           Form(";%s-TOF10 mean per spill (ns)",tree_name[chi].c_str()),
                           100,-100,100);

    hspillstdact[chi] = new TH1D(Form("%sspillstd",tree_name[chan[chi]].c_str()),
                                 Form(";%s-TOF10 stddev per spill (ns)",tree_name[chan[chi]].c_str()),
                                 100,0,20);

    hactcor[chi] = new TH2D(Form("%sactcor",tree_name[chan[chi]].c_str()),
                                 Form(";%s-TOF01 (ns);ACT3R-TOF10",tree_name[chan[chi]].c_str()),
                                 100,-100,100,100,-100,100);
  }

  // act and lg correction
  TH1D * httcoract = new TH1D("ttcoract",
                              ";TT0-TT1 (ns)",
                              100,-100,100);
  TH1D * httcorlg  = new TH1D("ttcorlg",
                              ";TT2-TT1 (ns)",
                              100,-100,100);

  // save for each spill
  vector< map<int,double> > cor_mean(nchans);
  vector< map<int,double> > cor_std(nchans);

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
        hspillcoract[chi]->Fill(hactt10cortmp[chi]->GetMean());
        hspillstdact[chi]->Fill(hactt10cortmp[chi]->GetStdDev());
        cor_mean[chi][spillNumber_pre] = hactt10cortmp[chi]->GetMean();
        cor_std[chi][spillNumber_pre] = hactt10cortmp[chi]->GetStdDev();
        if (hactt10cortmp[chi]->GetBinContent(0)>0 ||
            hactt10cortmp[chi]->GetBinContent(hactt10cortmp[chi]->GetNbinsX())>0) {
          cout << "spill number " << spillNumber_pre << endl;
          cout << tree_name[chi] << " ";
          cout << hactt10cortmp[chi]->GetMean() << " ";
          cout << hactt10cortmp[chi]->GetStdDev() << " ";
          cout << hactt10cortmp[chi]->GetBinContent(0) << " ";
          cout << hactt10cortmp[chi]->GetBinContent(hactt10cortmp[chi]->GetNbinsX());
          cout << endl;
        }
        hactt10cortmp[chi]->Reset();
      }

      for (int chi=0;chi<8;chi++) {
        hactcor[chi]->Fill(cor_mean[chi][spillNumber_pre],cor_mean[7][spillNumber_pre]);
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

    double ttcoract = tt_act-tt_tof;
    double ttcorlg  = tt_lg -tt_tof;

    httcoract->Fill(ttcoract);
    httcorlg ->Fill(ttcorlg);

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

      // e-like selection
      if (tof<12 && lg>0.2) {

        // act digitizer
        for (int acti=0; acti<8; acti++) {
          if (pmt[acti]->nPeaks==1 &&
              pmt[acti]->SignalTime[0]<200 &&
              pmt[acti]->PeakVoltage[0]>thresholdv[acti]) {
              //pmt[acti]->IntCharge[0]>thresholdq[acti]) {
            double actdiff = pmt[acti]->SignalTime[0]-pmt[12]->SignalTime[0];
            hactt10[acti]->Fill(actdiff);
            hactt10cortmp[acti]->Fill(actdiff+ttcoract);
          }
        }

        // lead glass digitizer
        double lgdiff = pmt[18]->SignalTime[0]-pmt[12]->SignalTime[0];
        hactt10[8]->Fill(lgdiff);
        hactt10cortmp[8]->Fill(lgdiff+ttcorlg);

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
    hspillcoract[chi]->Fill(hactt10cortmp[chi]->GetMean());
    hspillstdact[chi]->Fill(hactt10cortmp[chi]->GetStdDev());
    cor_mean[chi][info->SpillNumber] = hactt10cortmp[chi]->GetMean();
    cor_std[chi][info->SpillNumber] = hactt10cortmp[chi]->GetStdDev();
    cout << tree_name[chi] << " ";
    cout << "spill " << spillNumber_pre;
    cout << " mean (ns) " << hactt10cortmp[chi]->GetMean();
    cout << " std (ns) " <<  hactt10cortmp[chi]->GetStdDev();
    cout << " entries " <<  hactt10cortmp[chi]->GetEntries() << endl;
    cout << tree_name[chi] << " ";
    cout << "spill " << spillNumber_pre;
    cout << " mean (ns) " << hactt10[chi]->GetMean();
    cout << " std (ns) " <<  hactt10[chi]->GetStdDev();
    cout << " entries " <<  hactt10[chi]->GetEntries() << endl;
    if (hactt10cortmp[chi]->GetBinContent(0)>0 ||
        hactt10cortmp[chi]->GetBinContent(hactt10cortmp[chi]->GetNbinsX())>0) {
      cout << "spill number " << info->SpillNumber << endl;
      cout << tree_name[chi] << " ";
      cout << hactt10cortmp[chi]->GetMean() << " ";
      cout << hactt10cortmp[chi]->GetStdDev() << " ";
      cout << hactt10cortmp[chi]->GetBinContent(0) << " ";
      cout << hactt10cortmp[chi]->GetBinContent(hactt10cortmp[chi]->GetNbinsX());
      cout << endl;
    }
    hactt10cortmp[chi]->Reset();
  }

  TCanvas * c;

  c = new TCanvas();
  htoflg->SetStats(0);
  htoflg->Draw("colz");
  //c->Print(Form("%s_tof_lg.png",plotout.c_str()));

  c = new TCanvas();
  hlgact23->SetStats(0);
  hlgact23->Draw("colz");
  //c->Print(Form("%s_lg_act23.png",plotout.c_str()));

  for (int chi=0; chi<19; chi++) {
    c = new TCanvas();
    hv[chi]->SetStats(1);
    hv[chi]->Draw();
    //c->Print(Form("%s_%s_V.png",plotout.c_str(),tree_name[chi].c_str()));
  
    c = new TCanvas();
    hq[chi]->SetStats(1);
    hq[chi]->Draw();
    //c->Print(Form("%s_%s_Q.png",plotout.c_str(),tree_name[chi].c_str()));
  }

  c = new TCanvas();
  httcoract->Draw();
  //c->Print(Form("%s_ttcoract.png",plotout.c_str()));

  c = new TCanvas();
  httcorlg->Draw();
  //c->Print(Form("%s_ttcorlg.png",plotout.c_str()));

  for (int chi=0;chi<nchans;chi++) {
    c = new TCanvas();
    hspillcoract[chi]->Draw();
    //c->Print(Form("%s_%s_T10_spillcor.png",plotout.c_str(),tree_name[chi].c_str()));

    c = new TCanvas();
    hspillstdact[chi]->Draw();
    //c->Print(Form("%s_%s_T10_spillstd.png",plotout.c_str(),tree_name[chi].c_str()));

    c = new TCanvas();
    hactcor[chi]->Draw("colz");
    //c->Print(Form("%s_%s_T10_actcor.png",plotout.c_str(),tree_name[chi].c_str()));
  }

  // mean and std corrections
  //for (int i=0; i<9; i++) {
  //  cout << "channel " << i << endl;
  //  cout << "saved spill corrections size " << cor_mean[i].size() << endl;
  //  for (map<int,double>::iterator sp=cor_mean[i].begin(); sp!=cor_mean[i].end(); sp++) {
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

    double t_tof = pmt[12]->SignalTime[0];
    double t_act[8];
    for (int acti=0; acti<8; acti++) {
      t_act[acti] = pmt[acti]->SignalTime[0] + ttcoract - cor_mean[acti][info->SpillNumber];
    }
    double t_lg = pmt[18]->SignalTime[0] + ttcorlg - cor_mean[8][info->SpillNumber];

    for (int ipmt=0; ipmt<npmts; ipmt++) {
      for (int ipeak=0; ipeak<pmt[ipmt]->nPeaks; ipeak++) {
        if (ipmt<8) {
          signalTimeCor[ipmt][ipeak] = pmt[ipmt]->SignalTime[ipeak] + ttcoract - cor_mean[6][info->SpillNumber];
        }
        else if (ipmt<16) {
          signalTimeCor[ipmt][ipeak] = t_tof;
        }
        else {
          signalTimeCor[ipmt][ipeak] = pmt[18]->SignalTime[ipeak] + ttcorlg - cor_mean[8][info->SpillNumber];
        }
      }
      newtree[ipmt]->Fill();
    }

    // TOF cuts
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
            hactt10cor[acti]->Fill(t_act[acti]-t_tof);
          }
        }

        hactt10cor[8]->Fill(t_lg-t_tof);

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

  for (int chi=0; chi<nchans; chi++) {

    c = new TCanvas();
    THStack * s = new THStack();
    hactt10[chi]->SetStats(0);
    hactt10cor[chi]->SetStats(0);
    hactt10[chi]->SetLineColor(kBlue);
    hactt10cor[chi]->SetLineColor(kRed);
    s->Add(hactt10[chi]);
    s->Add(hactt10cor[chi]);
    hactt10[chi]->SetTitle("e-like");
    hactt10cor[chi]->SetTitle("e-like corr.");
    s->Draw("nostack");
    s->GetXaxis()->SetTitle(hactt10[chi]->GetXaxis()->GetTitle());
    s->GetYaxis()->SetTitle(hactt10[chi]->GetYaxis()->GetTitle());
    gPad->BuildLegend(0.6,0.7,0.88,0.89);
    //gPad->SetLogy();          
    c->Print(Form("%s_%s_T10_cor_stack.png",plotout.c_str(),tree_name[chan[chi]].c_str()));
    //c->Print(Form("%s_%s_T10_cor_stack_logy.png",plotout.c_str(),tnames[chi].c_str()));

  }

  //c = new TCanvas();
  //glgt10->Draw("ap");
  //glgt10->GetYaxis()->SetTitle("LG-T10 time (ns)");
  //glgt10->GetXaxis()->SetTitle("Sample number");
  //glgt10->Fit("pol1","Q");
  //glgt10->Draw("p");
  ////c->Print(Form("%s_LG_T10_drift.png",plotout.c_str()));

  //c = new TCanvas();
  //hlgt10->SetStats(1);
  //hlgt10->Draw();
  //c->Print(Form("%s_LG_T10.png",plotout.c_str()));
  ////gPad->SetLogy();          
  ////c->Print(Form("%s_LG_T10_logy.png",plotout.c_str()));
/*
  //c = new TCanvas();
  //hlgt10cor->SetStats(0);
  //hlgt10cor->SetTitle("");
  //hlgt10cor->SetLineColor(kRed);
  //hlgt10cor2->SetLineColor(kBlack);
  hlgt10cor->Draw();
  hlgt10cor2->Draw("same");
  gPad->BuildLegend(0.6,0.7,0.88,0.89);
  gPad->SetLogy();          
  c->Print(Form("%s_LG_T10_cor.png",plotout.c_str()));

  c = new TCanvas();
  THStack * s = new THStack();
  hlgt10->SetStats(0);
  hlgt10cor->SetStats(0);
  hlgt10->SetLineColor(kBlue);
  hlgt10cor->SetLineColor(kRed);
  s->Add(hlgt10);
  s->Add(hlgt10cor);
  hlgt10->SetTitle("e-like");
  hlgt10cor->SetTitle("e-like corr.");
  s->Draw("nostack");
  s->GetXaxis()->SetTitle(hlgt10->GetXaxis()->GetTitle());
  s->GetYaxis()->SetTitle(hlgt10->GetYaxis()->GetTitle());
  gPad->BuildLegend(0.6,0.7,0.88,0.89);
  //gPad->SetLogy();          
  c->Print(Form("%s_LG_T10_cor_stack.png",plotout.c_str()));
*/
  // copy results to a text file
  ofstream fout("timeCorrection.txt",ofstream::app);
  fout << run;
  for (int chi=0;chi<nchans;chi++) {
    fout << " " << hactt10[chi]->GetStdDev();
    fout << " " << hactt10cor[chi]->GetStdDev();
  }
  fout << endl;
  fout.close();

  // save new ntuple
  newfile->Write();

}
