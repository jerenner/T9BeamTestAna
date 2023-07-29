#ifndef MakeAllDataPlots_h
#define MakeAllDataPlots_h


#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSystem.h"
#include "TString.h"

#include <string>
#include <vector>
#include <iostream>


#include "EventInfo.h"
#include "channelReadClass.h"

#include <string>

using namespace std;

const int nMaxChannels = 32;


class MakeAllDataPlots
{


 private:

  bool _isHodoscopeRun;
  
  Double_t peakVoltage[nMaxChannels][1];
  Double_t peakTime[nMaxChannels][1];
  Double_t signalTime[nMaxChannels][1];
  Double_t intCharge[nMaxChannels][1];
  Double_t pedestal[nMaxChannels];
  Double_t pedestalSigma[nMaxChannels];
  Double_t nPeaks[nMaxChannels];

  map<TString,int> PeakID;
  map<TString,double> Amplitudes; // amplitude
  map<TString,double> Charges; // charge
  map<TString,double> SignalTimes; // time
    
  
  string _fileName;
  int _momentum;
  TString _peakMode;
  TFile *_infile;
  TFile *_outFile;
  EventInfo *_eventInfo;

  int _nChannels;
  vector<TString> _treeNames;
  
  // standard per channel
  vector<TH1D> hCharge;
  vector<TH1D> hVoltage;
  vector<TH1D> hPedestalSigma;
  vector<TH1D> hTime;
  vector<TH1D> hnPeaks;
  
  map<TString,TH1D*> _histos1d;
  map<TString,TH2D*> _histos2d;

  int _ent[nMaxChannels];
  channelReadClass *_reader[nMaxChannels];
  TTree *_trees[nMaxChannels];
  map<TString, channelReadClass*> _readerMap;

  int _debug;
  double _t0;
  double _t1;
  double _tof;
  
 public:
  
  MakeAllDataPlots(string fileName, int momentum, bool isHodoscopeRun, TString peakMode = "");
  ~MakeAllDataPlots();

  int getHighestPeakIndex(channelReadClass *reader);
  void Init();
  void InitReaders();
  void InitGeneralHistos();
  void InitChargedHistos();
  void InitHodoscopeHistos();

  void ReadChannels();
  bool PassedPeakCuts();
  void FillTofHistos();
  void FillChargedHistos();
  void FillHodoscopeHistos();

  
  void Loop( int verbose = 10000, int debug = 0);
  void Terminate();














};

#endif
