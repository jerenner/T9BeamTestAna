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
  fGlobalChannelID =  chid;
  anaHist = hist;
  if(!anaHist){
    std::cout << "No histogram" << std::endl;
    return;
  }

  //Calculate the analysis and pedestal windows if not already
  //
    if(fAnaWindowT0>=0. && fAnaWindowT1>0.){
        int maxbin = anaHist->GetNbinsX();
        if( fAnaWindowT0<anaHist->GetXaxis()->GetBinLowEdge(1) ||
            fAnaWindowT0>anaHist->GetXaxis()->GetBinUpEdge(maxbin) ) {
            std::cout << "Analysis window T0 is out of range, setting to first bin" << std::endl;
            fAnaWindowT0Bin = 1;
        }
        if( fAnaWindowT1<anaHist->GetXaxis()->GetBinLowEdge(1) ||
            fAnaWindowT1>anaHist->GetXaxis()->GetBinUpEdge(maxbin) ) {
            std::cout << "Analysis window T1 is out of range, setting to last bin" << std::endl;
            fAnaWindowT1Bin = maxbin;
        }
        fAnaWindowT0Bin = anaHist->GetXaxis()->FindBin(fAnaWindowT0);
        fAnaWindowT1Bin = anaHist->GetXaxis()->FindBin(fAnaWindowT1);
    }

    if(fPedWindowT0>=0. && fPedWindowT1>0.){
        int maxbin = anaHist->GetNbinsX();
        if( fPedWindowT0<anaHist->GetXaxis()->GetBinLowEdge(1) ||
            fPedWindowT0>anaHist->GetXaxis()->GetBinUpEdge(maxbin) ) {
            std::cout << "Pedestal window T0 is out of range, setting to first bin" << std::endl;
            fPedWindowT0Bin = 1;
        }
        if( fPedWindowT1<anaHist->GetXaxis()->GetBinLowEdge(1) ||
            fPedWindowT1>anaHist->GetXaxis()->GetBinUpEdge(maxbin) ) {
           // std::cout<< fPedWindowT1 << " " << anaHist->GetXaxis()->GetBinLowEdge(1) << " " << anaHist->GetXaxis()->GetBinUpEdge(maxbin)<< std::endl;
            std::cout << "Pedestal window T1 is out of range, setting to 1/3 of range" << std::endl;
            fPedWindowT1Bin = maxbin/3;
        }
        fPedWindowT0Bin = anaHist->GetXaxis()->FindBin(fPedWindowT0);
        fPedWindowT1Bin = anaHist->GetXaxis()->FindBin(fPedWindowT1);
    }

    fPedestal = 0.;
    fPedestalSigma=0;
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
 
    if(fPeakBins.size() == 0){
        double peakVoltage = 0;
        int maxbin = 0;
        for(int i = fAnaWindowT0Bin; i < fAnaWindowT1Bin; i++){
            double voltage = anaHist->GetBinContent(i)*fVoltScale-fPedestal;
            
            if(fPolarity==0) voltage = -1.0*voltage;


            if(voltage>peakVoltage){
                peakVoltage = voltage;
                maxbin = i;
            }
        }
        

        fPeakBins.push_back(maxbin);    
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
		    if (isToFChannel) {
		      int itof = (fGlobalChannelID - 8) / 4;
		      int jtof = (fGlobalChannelID - 8) % 4;
		      if (itof == 0 && jtof != 0) {
			int ktof = jtof - 1;
			calibt += ToF_channels_offset_calibration[itof][ktof];
		      } else if (itof == 1 && jtof != 1) {
			int ktof = jtof;
			if (ktof > 1)
			  ktof -= 1;
			calibt += ToF_channels_offset_calibration[itof][ktof];
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
            }
            
            voltage = anaHist->GetBinContent(fPeakBins.at(i)+1)*fVoltScale-fPedestal;
            counter = fPeakBins.at(i)+1;
            while(voltage > 3*fPedestalSigma && counter < anaHist->GetNbinsX()){
                integratedCharge += voltage*period/impedance;
                counter++; 
                voltage = anaHist->GetBinContent(counter)*fVoltScale-fPedestal;  
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
        std::cout << "Analysis window T0 is out of range, setting to first bin" << std::endl;
        fAnaWindowT0Bin = 1;
    }
    if( fAnaWindowT1<anaHist->GetXaxis()->GetBinLowEdge(1) ||
        fAnaWindowT1>anaHist->GetXaxis()->GetBinUpEdge(maxbin) ) {
        std::cout << "Analysis window T1 is out of range, setting to last bin" << std::endl;
        fAnaWindowT1Bin = maxbin;
    }
    fAnaWindowT0Bin = anaHist->GetXaxis()->FindBin(fAnaWindowT0);
    fAnaWindowT1Bin = anaHist->GetXaxis()->FindBin(fAnaWindowT1);
  } else {
    std::cout << "No histogram, so bin range will be set when histogram is loaded" << std::endl;
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
        std::cout << "Pedestal window T0 is out of range, setting to first bin" << std::endl;
        fPedWindowT0Bin = 1;
    }
    if( fPedWindowT1<anaHist->GetXaxis()->GetBinLowEdge(1) ||
        fPedWindowT1>anaHist->GetXaxis()->GetBinUpEdge(maxbin) ) {
        std::cout << "Pedestal window T1 is out of range, setting to 1/3 of range" << std::endl;
        fPedWindowT1Bin = maxbin/3;
    }
    fPedWindowT0Bin = anaHist->GetXaxis()->FindBin(fPedWindowT0);
    fPedWindowT1Bin = anaHist->GetXaxis()->FindBin(fPedWindowT1);
  } else {
    std::cout << "No histogram, so bin range will be set when histogram is loaded" << std::endl;
  }


}

// _______________________________________________________________________________
// _______________________________________________________________________________
// _______________________________________________________________________________
