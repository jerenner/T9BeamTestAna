#define WaveformAnalysis_cxx
#include "WaveformAnalysis.h"
#include <fstream>
#include <string>
#include <TStyle.h>
#include <TMath.h>
#include <TGraph.h>
#include "math.h"
#include "stdlib.h"
#include "TH1D.h"

using namespace std;

// _______________________________________________________________________________

WaveformAnalysis::WaveformAnalysis() {

    anaHist = NULL;
    fAnaWindowT0Bin=-1;
    fAnaWindowT1Bin=-1;
    fPedWindowT0Bin=-1;
    fPedWindowT1Bin=-1;

    fAnaWindowT0=-1.;
    fAnaWindowT1=-1.;
    fPedWindowT0=-1.;
    fPedWindowT1=-1.;

    fThreshold = -1.;
    fPolarity = 0; //negative polarity by default

    fGlobalChannelID = -1;
 
    fPedestal = 0.;

    fVoltScale = 2.5/4096.;

    fMeasureCharge = true;
    fMeasureTime = false;
    gaussFunc = new TF1("gaussFunc","[0]+[1]*TMath::Exp(-0.5*pow((x-[2])/[3],2))",0,3000);

}

// _______________________________________________________________________________

WaveformAnalysis::~WaveformAnalysis(){

}


// _______________________________________________________________________________

void WaveformAnalysis::SetHistogram(TH1D *hist, int chid){
  fGlobalChannelID = chid;
  anaHist = hist;
  if(!anaHist){
    std::cout << "WaveformAnalysis::SetHistogram: No histogram" << std::endl;
    return;
  }

  //Calculate the analysis and pedestal windows if not already
  //
  int maxbin = anaHist->GetNbinsX();
   double lowt = lowt;
   double hight = anaHist->GetXaxis()->GetBinUpEdge(maxbin);
  
    if(fAnaWindowT0>=0. && fAnaWindowT1>0.){
        int maxbin = anaHist->GetNbinsX();
        if( fAnaWindowT0<lowt ||fAnaWindowT0>hight ) {
	  if (fGlobalChannelID >= 0) std::cout << "Channel " << fGlobalChannelID << " Analysis window T0 " << fAnaWindowT0 << " is out of range " << lowt << "," << hight << " , setting to first bin" << std::endl;
            fAnaWindowT0Bin = 1;
        }
        if( fAnaWindowT1<lowt || fAnaWindowT1>hight ) {
	  if (fGlobalChannelID >= 0) std::cout << "Channel " << fGlobalChannelID << " Analysis window T1 " << fAnaWindowT1 << " is out of range " << lowt << "," << hight << " , setting to last bin" << std::endl;
	  fAnaWindowT1Bin = maxbin;
        }
        fAnaWindowT0Bin = anaHist->GetXaxis()->FindBin(fAnaWindowT0);
        fAnaWindowT1Bin = anaHist->GetXaxis()->FindBin(fAnaWindowT1);
    }

    if(fPedWindowT0>=0. && fPedWindowT1>0.){
//        int maxbin = anaHist->GetNbinsX();
        if( fPedWindowT0<lowt ||
            fPedWindowT0>hight ) {
	  if (fGlobalChannelID >= 0) std::cout << "Channel " << fGlobalChannelID << " Pedestal window T0 " << fAnaWindowT0 << " is out of range " << lowt << "," << hight << " , setting to first bin" << std::endl;
            fPedWindowT0Bin = 1;
        }
        if( fPedWindowT1<lowt ||
            fPedWindowT1>hight ) {
           // std::cout<< fPedWindowT1 << " " << lowt << " " << hight<< std::endl;
	  if (fGlobalChannelID >= 0) std::cout << "Channel " << fGlobalChannelID << " Pedestal window T1 " << fAnaWindowT1 << " is out of range " << lowt << "," << hight << " , setting to 1/3 range" << std::endl;
            fPedWindowT1Bin = maxbin/3;
        }
        fPedWindowT0Bin = anaHist->GetXaxis()->FindBin(fPedWindowT0);
        fPedWindowT1Bin = anaHist->GetXaxis()->FindBin(fPedWindowT1);
    }

    fPedestal = 0.;
    fPedestalSigma = 0;
    //cout << "++++++++++++++++++++++++++++++" << endl;
    //cout << anaHist->GetNbinsX() << endl;
    //cout << fPedWindowT0 << " " << fPedWindowT1 << " " << fPedWindowT0Bin << " " << fPedWindowT1Bin << endl;
    //cout << fAnaWindowT0 << " " << fAnaWindowT1 << " " << fAnaWindowT0Bin << " " << fAnaWindowT1Bin << endl;
    for(int i = fPedWindowT0Bin; i < fPedWindowT1Bin; i++){
        fPedestal += fVoltScale*anaHist->GetBinContent(i)/(fPedWindowT1Bin-fPedWindowT0Bin);
    }
    for(int i = fPedWindowT0Bin; i < fPedWindowT1Bin; i++){
        fPedestalSigma += TMath::Power(fVoltScale*anaHist->GetBinContent(i)-fPedestal, 2)/(fPedWindowT1Bin-fPedWindowT0Bin);
    }
    fPedestalSigma = sqrt(fPedestalSigma);
    //cout << fPedestal << endl;
}


// _______________________________________________________________________________

void WaveformAnalysis::FindPeaks(){

    fPeakBins.clear();
    fPeakVoltage.clear();
    fPeakTime.clear();
    // fNbPeaks.clear();
    fNbPeaks = 0;
    int nbPeaks = 0; //default value to be overwritten
    double maxVal = 0;
    int max = -1;
    int len = 0;
    /*for(int i = fAnaWindowT0Bin; i < fAnaWindowT1Bin; i++){
        double voltage = anaHist->GetBinContent(i)*fVoltScale-fPedestal;
        
        if(fPolarity==0) voltage = -1.0*voltage;
        
        if(voltage > 10*fPedestalSigma){
            len++;
            if(voltage > maxVal){
                maxVal = voltage;
                max = i;
                
            }
        }
        else if(max > 0){
            if(len>6){
                fPeakBins.push_back(max);
                fPeakVoltage.push_back(maxVal);
                fPeakTime.push_back(anaHist->GetXaxis()->GetBinCenter(max));
            }
            len = 0;
            max = -1;
            maxVal = 0;
        }
    }*/


    //counting the number of peaks
    bool goingDown = true;
    double pedestalFluctuations = 5 * fPedestalSigma; //arbirary number - need to change that
    double previousVoltage = fPedestal;
 
    if(fPeakBins.size() == 0){

        double peakVoltage = 0;
        int maxbin = 0;
        nbPeaks = 0;

        for(int i = fAnaWindowT0Bin; i < fAnaWindowT1Bin; i++){
            double voltage = anaHist->GetBinContent(i)*fVoltScale-fPedestal;
            
            if(fPolarity==0) voltage = -1.0*voltage;


            if(voltage>peakVoltage){
                peakVoltage = voltage;
                maxbin = i;
            }
            //we want to count the number of peaks
            //go through the waveform, first going down

            //if (voltage - previousVoltage > pedestalFluctuations and goingDown = true) continue;
            if (voltage - previousVoltage >= pedestalFluctuations and goingDown == true){
                //here, the current voltage is smaller than the previous voltage by more than the pedestal fluctuations: we have found a peak
                nbPeaks += 1;
                goingDown = false;

            }
            if (voltage - previousVoltage <= -pedestalFluctuations and goingDown == false){
                //here, the current voltage is larger than the previous voltage by more than the pedestal fluctuations: we have finished a peak, will need to look for the next one
                goingDown = true;

            }


            //std::cout << "In WaveformAnalysis.cc: peakFinder info: " << nbPeaks << std::endl;
            previousVoltage = voltage;

        }
        

        fPeakBins.push_back(maxbin);
        fNbPeaks = nbPeaks;
        fPeakVoltage.push_back(peakVoltage);
        fPeakTime.push_back(anaHist->GetXaxis()->GetBinCenter(maxbin));
    }   


}

// _______________________________________________________________________________

void WaveformAnalysis::CalculateSignalTime(){
    fSignalTime.clear();
    double fraction = 0.4;
    //cout << fPeakBins.size() << endl;
    for(int i = 0; i < fPeakBins.size(); i++){
        int maxbin = fPeakBins.at(i);
        if(fMeasureTime){
            bool set = false;
            for(int j = maxbin; j > 0; j--){
                double amp = anaHist->GetBinContent(j)*fVoltScale-fPedestal;
                if(fPolarity==0) amp = -1.0*amp;
                if(amp <  fraction*fPeakVoltage.at(i)){
                    double x0 = anaHist->GetBinCenter(j);
                    double x1 = anaHist->GetBinCenter(j+1);
                    double y0 = fPedestal - fVoltScale*anaHist->GetBinContent(j);
                    double y1 = fPedestal - fVoltScale*anaHist->GetBinContent(j+1);
                    double a = (y1-y0)/(x1-x0);
                    double b = y0 - a*x0;

		    double uncalibt = (fraction*fPeakVoltage.at(i)-b)/a;
		    double calibt = uncalibt;

		    bool isToFChannel = (fGlobalChannelID >= 8 && fGlobalChannelID <= 15);
		    if (isToFChannel) {                                    // channels  8  9 10 11 12 13 14 15
		      int itof = (fGlobalChannelID - 8) / 4;               //     itof  0  0  0  0  1  1  1  1
		      int jtof = (fGlobalChannelID - 8) % 4;               //     jtof  0  1  2  3  0  1  2  3
		      if (itof == 0 && jtof != 0) {                        //              x  x  x
			int ktof = jtof - 1;                                   //     ktof     0  1  2
			calibt += ToF_channels_offset_calibration[itof][ktof]; //       ik    00 01 02
		      } else if (itof == 1 && jtof != 1) {                 //                       x     x  x
			int ktof = jtof;                                       //     ktof              0     2  3
			if (ktof > 1)                                          //                             x  x
			  ktof -= 1;                                           //     ktof              0     1  2
			calibt += ToF_channels_offset_calibration[itof][ktof]; //       ik             10    11 12
		      }
		    }
		    
                    fSignalTime.push_back(calibt);
                    set = true;
                    break;
                }
            }
            if(!set) fSignalTime.push_back(0);
        }
        else 
            fSignalTime.push_back(0);
        
    }
}

// _______________________________________________________________________________

void WaveformAnalysis::IntegrateCharge(){
    fIntegratedCharge.clear();
    double impedance = 50;
    for(int i = 0; i < fPeakBins.size(); i++){
        if(fMeasureCharge){
           
            
            int counter = fPeakBins.at(i);
            double voltage = anaHist->GetBinContent(fPeakBins.at(i))*fVoltScale-fPedestal;
            if(fPolarity==0) voltage = -1.0*voltage;
            
            double period = anaHist->GetXaxis()->GetBinWidth(1);
            double integratedCharge = 0;
            while(voltage > 3*fPedestalSigma && counter > 1){
                integratedCharge += voltage*period/impedance;
                counter--; 
                voltage = anaHist->GetBinContent(counter)*fVoltScale-fPedestal;  
                if(fPolarity==0) voltage = -1.0*voltage;
            }

            voltage = anaHist->GetBinContent(fPeakBins.at(i)+1)*fVoltScale-fPedestal;
            if(fPolarity==0) voltage = -1.0*voltage;
            counter = fPeakBins.at(i)+1;
            while(voltage > 3*fPedestalSigma && counter < anaHist->GetNbinsX()){
                integratedCharge += voltage*period/impedance;
                counter++; 
                voltage = anaHist->GetBinContent(counter)*fVoltScale-fPedestal;  
                if(fPolarity==0) voltage = -1.0*voltage;
            }  

            fIntegratedCharge.push_back(integratedCharge);
        }   
        else  
            fIntegratedCharge.push_back(0);
    }

}


// _______________________________________________________________________________

void WaveformAnalysis::RunAnalysis(){
    FindPeaks();
    CalculateSignalTime();
    IntegrateCharge();
    
    for(int i = 0; i < fPeakBins.size(); i++){
        if(fPeakVoltage.at(i) > fThreshold){
            fIsOverThreshold.push_back(true);
        }
        else fIsOverThreshold.push_back(false);
    }

}

// _______________________________________________________________________________

void WaveformAnalysis::SetAnalysisBinWindow(double t0, double t1){

  fAnaWindowT0 = t0;
  fAnaWindowT1 = t1;

  if(anaHist){
    int maxbin = anaHist->GetNbinsX();
    if( fAnaWindowT0<anaHist->GetXaxis()->GetBinLowEdge(1) ||
        fAnaWindowT0>anaHist->GetXaxis()->GetBinUpEdge(maxbin) ) {
      if (fGlobalChannelID >= 0) std::cout << "Channel " << fGlobalChannelID << " Analysis window T0 is out of range, setting to first bin" << std::endl;
        fAnaWindowT0Bin = 1;
    }
    if( fAnaWindowT1<anaHist->GetXaxis()->GetBinLowEdge(1) ||
        fAnaWindowT1>anaHist->GetXaxis()->GetBinUpEdge(maxbin) ) {
      if (fGlobalChannelID >= 0) std::cout << "Channel " << fGlobalChannelID << " Analysis window T1 is out of range, setting to last bin" << std::endl;
        fAnaWindowT1Bin = maxbin;
    }
    fAnaWindowT0Bin = anaHist->GetXaxis()->FindBin(fAnaWindowT0);
    fAnaWindowT1Bin = anaHist->GetXaxis()->FindBin(fAnaWindowT1);
  } else {
    if (fGlobalChannelID >= 0) std::cout << "WaveformAnalysis::SetAnalysisBinWindow: Channel " << fGlobalChannelID << " No histogram, so bin range will be set when histogram is loaded" << std::endl;
  }


}

// _______________________________________________________________________________

void WaveformAnalysis::SetPedestalBinWindow(double t0, double t1){

  fPedWindowT0 = t0;
  fPedWindowT1 = t1;

  if(anaHist){
    int maxbin = anaHist->GetNbinsX();
    if( fPedWindowT0<anaHist->GetXaxis()->GetBinLowEdge(1) ||
        fPedWindowT0>anaHist->GetXaxis()->GetBinUpEdge(maxbin) ) {
      if (fGlobalChannelID >= 0) std::cout << "Channel " << fGlobalChannelID << " Pedestal window T0 is out of range, setting to first bin" << std::endl;
        fPedWindowT0Bin = 1;
    }
    if( fPedWindowT1<anaHist->GetXaxis()->GetBinLowEdge(1) ||
        fPedWindowT1>anaHist->GetXaxis()->GetBinUpEdge(maxbin) ) {
      if (fGlobalChannelID >= 0)  std::cout << "Channel " << fGlobalChannelID << " Pedestal window T1 is out of range, setting to 1/3 of range" << std::endl;
        fPedWindowT1Bin = maxbin/3;
    }
    fPedWindowT0Bin = anaHist->GetXaxis()->FindBin(fPedWindowT0);
    fPedWindowT1Bin = anaHist->GetXaxis()->FindBin(fPedWindowT1);
  } else {
    if (fGlobalChannelID >= 0) std::cout << "WaveformAnalysis::SetPedestalBinWindow: Channel " << fGlobalChannelID << " No histogram, so bin range will be set when histogram is loaded" << std::endl;
  }


}

// _______________________________________________________________________________
// _______________________________________________________________________________
// _______________________________________________________________________________
