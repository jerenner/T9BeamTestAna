#include <iostream>
#include <fstream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TString.h"

double GetBeta(double mass, double mom){
    double bg = mom/mass;
    double beta = sqrt(bg*bg/(1+bg*bg));
    return beta;
}

void FitTOF(string fileName, int p){
    double m[4]={0.511, 105.66, 139.57, 938.470};
    double c = 0.299792458;
    double l = 2.9;
    TFile *inFile = new TFile(fileName.c_str(), "READ");
    
    TH1D* h1 = (TH1D*) inFile->Get("hTOFAll");
    TH1D* h2 = (TH1D*) inFile->Get("hTOFEl");
    TH1D* h3 = (TH1D*) inFile->Get("hTOFOther");
    
    string title = "p = " + to_string(p) + " MeV/c";
    h3->SetTitle(title.c_str());
    h3->GetXaxis()->SetTitle("t_{1}-t_{0} [ns]");
    h3->GetYaxis()->SetTitle("Events");
    h3->SetStats(kFALSE);
    
    TF1* func = new TF1("func", "[0]*TMath::Gaus(x, [2], [3], 1) + [1]*TMath::Gaus(x, [5], [7], 1) + [4]*TMath::Gaus(x, [6], [7], 1)", 38, 43.5);
    
    func->SetParName(0, "norm1");
    func->SetParameter(0, 2000);
    func->SetParLimits(0, 0, 20000);
    
    func->SetParName(1, "norm2");
    func->SetParameter(1, 200);
    func->SetParLimits(1, 0, 10000);
    
    
    func->SetParName(2, "mean_e");
    func->SetParameter(2, h2->GetMean());
    func->SetParLimits(2, 0.99*h2->GetMean(), 1.01*h2->GetMean());
    
    func->SetParName(3, "sigma_e");
    //func->FixParameter(2, 0.215);
    func->SetParameter(3, h2->GetStdDev());
    func->SetParLimits(3, 0.4*h2->GetStdDev(), 1.2*h2->GetStdDev()); 
    
    func->SetParName(4, "norm3");
    func->SetParameter(4, 200);
    func->SetParLimits(4, 0, 1000);

    double val1 = h2->GetMean()+l/c/GetBeta(m[1], p)-l/c/GetBeta(m[0], p);
    func->SetParName(5, "mean_mu");
    //func->FixParameter(5, val1);
    func->SetParLimits(5, 0.99*val1, 1.01*val1);       
    
    /*func->SetParName(6, "sigma_mu");
    func->SetParameter(6, h2->GetStdDev());
    func->SetParLimits(6, 0.9*h2->GetStdDev(), 1.1*h2->GetStdDev()); */
    
    double val2 = h2->GetMean()+l/c/GetBeta(m[2], p)-l/c/GetBeta(m[0], p);
    func->SetParName(6, "mean_pi");
    //func->FixParameter(6, val2);
    func->SetParLimits(6, val2-0.2, val2+0.6);       
    
    func->SetParName(7, "sigma_mupi");
    //func->FixParameter(2, 0.215);
    func->SetParameter(7, h2->GetStdDev());
    func->SetParLimits(7, 0.6*h2->GetStdDev(), 1.2*h2->GetStdDev()); 
    
    h3->Fit(func);
    
    TCanvas *cv1 = new TCanvas("cv1", "", 1200, 800);
    cv1->cd()->SetTickx(kTRUE);
    cv1->cd()->SetTicky(kTRUE);
    //cv1->cd()->SetLogy(kTRUE);
    h3->SetMarkerStyle(20);
    h3->Draw("ep");
    
    func->SetLineWidth(2);
    func->SetLineColor(kBlack);
    func->Draw("csame");
    
    TF1* func1 = new TF1("func1", "[0]*TMath::Gaus(x, [1], [2], 1)", 38, 43.5);
    TF1* func2 = new TF1("func2", "[0]*TMath::Gaus(x, [1], [2], 1)", 38, 43.5);
    TF1* func3 = new TF1("func3", "[0]*TMath::Gaus(x, [1], [2], 1)", 38, 43.5);
    
    func1->SetParameter(0, func->GetParameter(0));
    func1->SetParameter(1, func->GetParameter(2));
    func1->SetParameter(2, func->GetParameter(3));
    
    func2->SetParameter(0, func->GetParameter(1));
    func2->SetParameter(1, func->GetParameter(5));
    func2->SetParameter(2, func->GetParameter(7));
    
    func3->SetParameter(0, func->GetParameter(4));
    func3->SetParameter(1, func->GetParameter(6));
    func3->SetParameter(2, func->GetParameter(7));
    
    func1->SetLineWidth(2);
    func1->SetLineColor(kRed+1);
    
    func2->SetLineWidth(2);
    func2->SetLineColor(kBlue+1);
    
    func3->SetLineWidth(2);
    func3->SetLineColor(kGreen+1);

    TLegend *leg1 = new TLegend(0.6, 0.6, 0.9, 0.9);

    
    int Ne = func->GetParameter(0)/h1->GetBinWidth(1);
    int Nmu = func->GetParameter(1)/h1->GetBinWidth(1);
    int Npi = func->GetParameter(4)/h1->GetBinWidth(1);
    
    int NeErr = func->GetParError(0)/h1->GetBinWidth(1);
    int NmuErr = func->GetParError(1)/h1->GetBinWidth(1);
    int NpiErr = func->GetParError(4)/h1->GetBinWidth(1);
    
    leg1->AddEntry(h3,"Data", "PE");
    leg1->AddEntry(func, "Fit", "L");
    title = "Fit - e^{+}, N = " + to_string(Ne) + " #pm " + to_string(NeErr);
    leg1->AddEntry(func1, title.c_str(), "L");
    title = "Fit - #mu^{+}, N = " + to_string(Nmu) + " #pm " + to_string(NmuErr);
    leg1->AddEntry(func2, title.c_str(), "L");
    title = "Fit - #pi^{+}, N = " + to_string(Npi) + " #pm " + to_string(NpiErr);
    leg1->AddEntry(func3, title.c_str(), "L");
    
    
    func1->Draw("Csame");
    func2->Draw("Csame");
    func3->Draw("Csame");
    
    leg1->Draw("same");
    string name = fileName.substr(0, fileName.size()-5) + "_TOFfit.png";
    cv1->Print(name.c_str());    

    cout << GetBeta(m[0], p) << " " << GetBeta(m[1], p) << " " << GetBeta(m[2], p) << " " << GetBeta(m[3], p) << endl;
    cout << l/c/GetBeta(m[0], p) << " " << l/c/GetBeta(m[1], p) << " " << l/c/GetBeta(m[2], p) << " " << l/c/GetBeta(m[3], p) << endl;
    cout << h2->GetMean() << " " << h2->GetMean()+l/c/GetBeta(m[1], p)-l/c/GetBeta(m[0], p) << " " << h2->GetMean()+l/c/GetBeta(m[2], p)-l/c/GetBeta(m[0], p) << endl;
        
    cout << h1->GetEntries() << " " << h2->GetEntries() << " " << h3->GetEntries() << endl;
    
    cout << "N_e = " <<  func->GetParameter(0)/h1->GetBinWidth(1) << endl;
    cout << "N_mu = " <<  func->GetParameter(1)/h1->GetBinWidth(1) << endl;
    cout << "N_pi = " <<  func->GetParameter(4)/h1->GetBinWidth(1) << endl;

    TString asciiname = "ascii_" + TString(inFile->GetName()).ReplaceAll(".root",".txt");
    ofstream *asciifile = new ofstream(asciiname.Data());
    (*asciifile) << "p " << p << endl;
    (*asciifile) << "N_e " << func->GetParameter(0)/h1->GetBinWidth(1) << endl;
    (*asciifile) << "N_mu " <<  func->GetParameter(1)/h1->GetBinWidth(1) << endl;
    (*asciifile) << "N_pi " <<  func->GetParameter(4)/h1->GetBinWidth(1) << endl;
    if (asciifile) asciifile->close();

    cout << "DONE!" << endl;
    cout << "See also " << asciiname.Data() << endl;
      
    
}
