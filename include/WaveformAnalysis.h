#ifndef WaveformAnalysis_h
#define WaveformAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "TH1D.h"
#include "TF1.h"

class WaveformAnalysis {

  public :

   WaveformAnalysis();
   ~WaveformAnalysis();
   void SetHistogram(TH1D *hist);
   void RunAnalysis();
   void SetAnalysisBinWindow(double t0, double t1);
   void SetPedestalBinWindow(double t0, double t1);
   void SetThreshold(double voltage){ fThreshold = voltage; };
   void SetPolarity(int pol){ fPolarity = pol; };
   void SetVoltageScale(double scale){ fVoltScale = scale; };
   void SetChargeMeasurement(bool val){fMeasureCharge = val;}
   void SetTimeMeasurement(bool val){fMeasureTime = val;}

   std::vector<double>& GetIntegratedCharge(){return fIntegratedCharge;};
   std::vector<double>& GetPeakVoltage(){return fPeakVoltage;};
   std::vector<double>& GetPeakTime(){return fPeakTime;};
   std::vector<double>& GetSignalTime(){return fSignalTime;};
   double GetPedestal(){return fPedestal;}
   double GetPedestalSigma(){return fPedestalSigma;}
   std::vector<bool>& GetThresholdStatus(){return fIsOverThreshold;}


  private :

   void FindPeaks();
   void CalculateSignalTime();
   void IntegrateCharge();
   std::vector<int> fPeakBins;
   
   TH1D *anaHist;

   double fAnaWindowT0, fAnaWindowT1;
   double fPedWindowT0, fPedWindowT1;
   int fAnaWindowT0Bin, fAnaWindowT1Bin;
   int fPedWindowT0Bin, fPedWindowT1Bin;
   double fThreshold;
   double fPolarity;
   
   std::vector<double> fIntegratedCharge;
   std::vector<double> fPeakVoltage;
   std::vector<double> fPeakTime;
   std::vector<double> fSignalTime;
   std::vector<bool> fIsOverThreshold;
   
   double fVoltScale;
   double fPedestal;
   double fPedestalSigma;
   

  TF1 *gaussFunc;

  bool fMeasureCharge;
  bool fMeasureTime;

};

#endif

