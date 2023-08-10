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
  bool _noAct1Cuts; // for ToF fits
  
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
  map<TString,double> NPeaks; // time

  // ranges
  double tofmin;
  double tofmax;
  int ntofbins;

  double tofminlow;
  double tofmaxlow;
  int ntofbinslow;

  int ntofbins2d;

  double actChargeMin;
  double actChargeMax;
  double actAmplitudeMax;

  
  // cuts
  map<int, map<TString,double > > _cutsMap;

  // IO
  
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

  // for peak multiplicity:
  bool _onePeakInAllACTs;
  bool _onePeakInAllToFs;
  bool _onePeakInAll;
  bool _onePeakInPbGlass;
  
  bool _moreThanOnePeakInAllACTs;
  bool _moreThanOnePeakInAllToFs;
  bool _moreThanOnePeakInAll;
  
  bool _PbGlassAboveElectronLevel;
  bool _ACT23AboveElectronLevel;

  // for tof fit:
  bool _isdACT23pb;
  bool _ispACT23pb;
  bool _isMuACT23pb;
  bool _isElACT23pb;
  
 public:
  
  MakeAllDataPlots(string fileName, int momentum, bool isHodoscopeRun, TString peakMode = "");
  ~MakeAllDataPlots();

  int getHighestPeakIndex(channelReadClass *reader);
  void Init(bool noAct1Cuts);
  void InitReaders();
  void InitTofHistos();
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
