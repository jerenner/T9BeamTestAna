#include "MakeAllDataPlots.C"

void runMakeAllDataPlots(string fileName, int momentum, bool isHodoscopeRun, TString peakMode = "") {

  MakeAllDataPlots *analysis = new MakeAllDataPlots(fileName, momentum, isHodoscopeRun, peakMode);
  analysis -> Init();
  analysis -> InitGeneralHistos();
  if (true) {
    cout << "Initializing charged particle analysis histos..." << endl;
    analysis ->  InitChargedHistos();
  } else {
    cout << "Initializing hodoscope analysis histos..." << endl;
    analysis ->  InitHodoscopeHistos();
  }
  analysis -> InitReaders();
  //  analysis -> Loop(10000, 1);
  analysis -> Loop();
  analysis -> Terminate();
}
