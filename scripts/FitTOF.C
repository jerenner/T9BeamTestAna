#include <iostream>
#include <fstream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TString.h"

const double cc = 299792458.;
  
double GetBeta(double mass, double mom){
    double bg = mom/mass;
    double beta = sqrt(bg*bg/(1+bg*bg));
    return beta;
}

  double GetTofDiffWrtE(double mass, double mom, double me, double L = 2.9){
    return L/cc* ( sqrt(1.+pow(mass/mom,2)) - sqrt(1.+pow(me/mom,2))*1) * 1e9;
}

void FitTOF(string fileName, int p){
    double m[4]={0.511, 105.66, 139.57, 938.470};
    double c = 0.299792458;
    double l = 2.9;

    string ptag = " (Pos)";
    if (p < 0.)
      ptag = " (Neg)";
    
    TFile *inFile = new TFile(fileName.c_str(), "READ");
    
    TH1D* h1 = (TH1D*) inFile->Get("hTOFAll");
    TH1D* h2 = (TH1D*) inFile->Get("hTOFEl");
    TH1D* h3 = (TH1D*) inFile->Get("hTOFOther");

    cout << "Integrals: " 
	 << " h1: " << h1->Integral()
      	 << " h2: " << h2->Integral()
      	 << " h3: " << h3->Integral()
	 << endl;
    cout << "Means: " 
	 << " h1: " << h1->GetMean()
      	 << " h2: " << h2->GetMean()
      	 << " h3: " << h3->GetMean()
	 << endl;
    cout << "StdDev: " 
	 << " h1: " << h1->GetStdDev()
      	 << " h2: " << h2->GetStdDev()
      	 << " h3: " << h3->GetStdDev()
	 << endl;

    TCanvas *cv0 = new TCanvas("cv0", "", 800, 800);
    cv0 -> Divide(2,2);


    // pre-fit to get initial pars for the full fit
    cv0->cd(1);
    TF1* func0 = new TF1("func0", "[0]*TMath::Gaus(x, [1], [2], 1)", 38, 43.5);
    h1 -> Draw("e1");
    double a = h1 -> GetMaximum(); 
    func0 -> SetParameters(a, h1->GetMean(), h1->GetStdDev());
    cout << "The initial parameters are: " << endl;
    for (int ip = 0; ip < func0 -> GetNpar(); ++ ip) {
      cout << ip << " " << func0->GetParName(ip) << " " << func0->GetParameter(ip) << endl;
    }
    h1 -> Fit("func0");

    // draw pref-fit:
    h1 -> SetMarkerStyle(20);
    h1 -> Draw("e1 hist same");
    func0 -> Draw("same");
    cv0 -> cd(2);
    h2 -> Draw("e1");
    cv0 -> cd(3);
    h3 -> Draw("e1");
    cv0 -> cd(4);

    // the main full fit:
   
    string title = "p = " + to_string(abs(p)) + " MeV/c" + ptag ;
    h3->SetTitle(title.c_str());
    h3->GetXaxis()->SetTitle("t_{1}-t_{0} [ns]");
    h3->GetYaxis()->SetTitle("Events");
    h3->SetStats(kFALSE);
    
    TF1* func = new TF1("func", "[0]*TMath::Gaus(x, [2], [3], 1) + [1]*TMath::Gaus(x, [5], [7], 1) + [4]*TMath::Gaus(x, [6], [7], 1)", 38, 43.5);

    double A = h1 -> GetMaximum();
    func->SetParName(0, "norm1");
    func->SetParameter(0, func0 -> GetParameter(0));
    //    func->SetParLimits(0, 0.1*A, 2.*A);
    
    func->SetParName(1, "norm2");
    func->SetParameter(1, h3 -> GetMaximum()/2.);
    //    func->SetParLimits(1, 0, 0.8*A);

    double delta = 0.05;    
    func->SetParName(2, "mean_e");
    func->SetParameter(2, func0 -> GetParameter(1));//h1->GetMean());
    //    func->SetParLimits(2, (1.-delta)*h1->GetMean(), (1. + delta)*h1->GetMean());
    
    func->SetParName(3, "sigma_e");
    //func->FixParameter(2, 0.215);
    func->SetParameter(3, func0 -> GetParameter(2));//h1->GetStdDev());
    //    func->SetParLimits(3, 0.4*h1->GetStdDev(), 1.2*h1->GetStdDev()); 
    
    func->SetParName(4, "norm3");
    func->SetParameter(4, h3 -> GetMaximum()/3.);
    //    func->SetParLimits(4, 0, 0.8*A);

    double val1 = h3->GetMean()+l/c/GetBeta(m[1], p)-l/c/GetBeta(m[0], p);
    func->SetParName(5, "mean_mu");
    func->SetParameter(5, val1);
    //func->FixParameter(5, val1);
    func->SetParLimits(5, (1. - delta)*val1, (1. + delta)*val1);       
    
    /*func->SetParName(6, "sigma_mu");
    func->SetParameter(6, h3->GetStdDev());
    func->SetParLimits(6, 0.9*h3->GetStdDev(), 1.1*h3->GetStdDev()); */
    
    double val2 = h3->GetMean()+l/c/GetBeta(m[2], p)-l/c/GetBeta(m[0], p);
    func->SetParName(6, "mean_pi");
    func->SetParameter(6, val2);
    //func->FixParameter(6, val2);
    func->SetParLimits(6,  (1. - delta)*val2,  (1. - delta)*val2);       
    
    func->SetParName(7, "sigma_mupi");
    //func->FixParameter(2, 0.215);
    func->SetParameter(7, h3->GetStdDev());
    //    func->SetParLimits(7, 0.5*h3->GetStdDev(), 1.3*h3->GetStdDev()); 

    cout << "The initial parameters are: " << endl;
    for (int ip = 0; ip < func -> GetNpar(); ++ ip) {
      cout << ip << " " << func->GetParName(ip) << " " << func->GetParameter(ip) << endl;      
    }
    
    h3->Fit(func);
    
    TCanvas *cv1 = new TCanvas("cv1", "", 1200, 800);
    cv1->cd()->SetTickx(kTRUE);
    cv1->cd()->SetTicky(kTRUE);
    //cv1->cd()->SetLogy(kTRUE);
    
    h3->SetMarkerStyle(20);
    h3->SetMaximum(1.5*h3->GetMaximum());
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

    TLegend *leg1 = new TLegend(0.60, 0.6, 0.9, 0.88);
    leg1 -> SetBorderSize(0);
    int ip = int(fabs(p));
    leg1 -> SetHeader("WCTE TB July 2022");
    
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
    string namepdf = fileName.substr(0, fileName.size()-5) + "_TOFfit.pdf";
    cv1->Print(name.c_str());
    cv1->Print(namepdf.c_str());    


    // TLines
    int ls = 2;
    int lw = 2;
    double y0 = h1->GetMinimum();
    int imax = h1->GetMaximumBin();
    double y1 = 1.1*h1->GetBinContent(imax);
    double xe = h1->GetBinCenter(imax);
    TLine *le = new TLine(xe, y0, xe, y1);
    le -> SetLineColor(kRed+1);
    le -> SetLineStyle(ls);
    le -> SetLineWidth(lw);
    le -> Draw();

    double dtmu = GetTofDiffWrtE(m[1], p, m[0]);
    TLine *lmu = new TLine(xe + dtmu, y0, xe + dtmu, y1);
    lmu -> SetLineColor(kBlue+1);
    lmu -> SetLineStyle(ls);
    lmu -> SetLineWidth(lw);
    lmu -> Draw();

    double dtpi = GetTofDiffWrtE(m[2], p, m[0]);
    TLine *lpi = new TLine(xe + dtpi, y0, xe + dtpi, y1);
    lpi -> SetLineColor(kGreen+1);
    lpi -> SetLineStyle(ls);
    lpi -> SetLineWidth(lw);
    lpi -> Draw();

    double dtp = GetTofDiffWrtE(m[3], p, m[0]);
    TLine *lp = new TLine(xe + dtp, y0, xe + dtp, y1);
    lp -> SetLineColor(kBlack);
    lp -> SetLineStyle(ls);
    lp -> SetLineWidth(lw);
    lp -> Draw();
    
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
