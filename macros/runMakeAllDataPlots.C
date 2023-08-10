#include "MakeAllDataPlots.C"

// jk 29.7.2023

void runMakeAllDataPlots(string fileName, int momentum, bool isHodoscopeRun, bool noAct1Cuts, TString peakMode = "") {

  MakeAllDataPlots *analysis = new MakeAllDataPlots(fileName, momentum, isHodoscopeRun, peakMode);
  analysis -> Init(noAct1Cuts);
  analysis -> InitGeneralHistos();
  analysis -> InitTofHistos();
  if (!isHodoscopeRun) {
    cout << "Initializing charged particle analysis histos..." << endl;
    analysis ->  InitChargedHistos();
  } else {
    cout << "Initializing hodoscope analysis histos..." << endl;
    analysis ->  InitHodoscopeHistos();
  }
  analysis -> InitReaders();
  // verbose, debug:
  //  analysis -> Loop(10000, 1);
  analysis -> Loop();
  analysis -> Terminate();
}
