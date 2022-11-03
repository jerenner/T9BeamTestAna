//
// Matej Pavin, modif. Jiri Kvita 2022
//

#include <iostream>
#include <fstream>

#include "TROOT.h"
#include "TStyle.h"

#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TLine.h"
#include "TString.h"
#include "TLatex.h"

// ______________________________________________________

const double cc = 299792458.;

// ______________________________________________________

double GetBeta(double mass, double mom) {
    double bg = mom/mass;
    double beta = sqrt(bg*bg/(1+bg*bg));
    return beta;
}

// ______________________________________________________

double GetTofDiffWrtE(double mass, double mom, double me, double L = 2.9) {
    return L/cc* ( sqrt(1.+pow(mass/mom,2)) - sqrt(1.+pow(me/mom,2))*1) * 1e9;
}

// ______________________________________________________

double IntegrateAlongYInGivenXWindow(TH2D *h2, double t0, double tof_reso) {
  double sum = 0.;
  int i1 = h2 -> GetXaxis() -> FindBin(t0 - tof_reso);
  int i2 = h2 -> GetXaxis() -> FindBin(t0 + tof_reso);
  cout << "Integrating "  << h2 -> GetName() << " with nx=" <<  h2 -> GetXaxis() -> GetNbins() << " in x bin range " << i1 << "," << i2 << endl;
  for (int i = i1; i <= i2; ++i) {
    for (int j = 0; j <= h2 -> GetYaxis() -> GetNbins(); ++j) {
      sum += h2 -> GetBinContent(i, j);
    }
  }
  return sum;
}

// ______________________________________________________
// ______________________________________________________


void FitTOF(string fileName, int p, bool logy = false) {
    // e, mu. pi, proton, deuteron, tritium, alpha
    double m[7] = {0.511, 105.66, 139.57, 938.470, 1876., 3.01604928*931.494102, 3727.379};
    double c = 0.299792458;
    double l = 2.9;

    bool isThreeComponentFit = fabs(p) < 299.;
    bool isProtonFit = fabs(p) > 390.;

    cout << "Configured as" << endl;
    cout << " isThreeComponentFit: " << isThreeComponentFit << endl 
	 << " isProtonFit: " << isProtonFit << endl 
	 << endl;
    
    string ptag = " (Pos)";
    if (p < 0.)
      ptag = " (Neg)";
    
    TFile *inFile = new TFile(fileName.c_str(), "READ");

    // mu, pi, protons
    double dtmu = GetTofDiffWrtE(m[1], p, m[0]);
    double dtpi = GetTofDiffWrtE(m[2], p, m[0]);
    double dtp = GetTofDiffWrtE(m[3], p, m[0]);
    // for fun also deuteron, triton, alpha
    double dtd = GetTofDiffWrtE(m[4], p, m[0]);
    double dtt = GetTofDiffWrtE(m[5], p, m[0]);
    double dta = GetTofDiffWrtE(m[6], p, m[0]);

    double tmin = 38.;
    double tmax = 43.5;
    if (isProtonFit) {
      tmax = 42 + dtp;
    }
    cout << "Will fit TOF distribution betweem " << tmin << " " << tmax << endl;

    string basehname = "hTOFOther";
    if (isProtonFit) {
      basehname = "hTOFAllWide"; 
    }
    
    TH1D* h1 = (TH1D*) inFile->Get("hTOFAll");
    TH1D* h2 = (TH1D*) inFile->Get("hTOFEl");
    TH1D* h3 = (TH1D*) inFile->Get(basehname.c_str());

    int imax = h1->GetMaximumBin();
    double xe = h1->GetBinCenter(imax);
    
    // time diffs of mu and pi to e
    //for (int ibin = 0; ibin < h1->GetNbinsX()+1; ++ibin) {
    //  cout << "hall i=" << ibin << " y=" << h1->GetBinContent(ibin) << endl;
    //}
    
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
    
    TF1* func0 = new TF1("func0", "[0]*TMath::Gaus(x, [1], [2], 1)", tmin, tmax);
    func0 -> SetNpx(1000);
    h1 -> Draw("e1");
    double a = h1 -> GetMaximum(); 
    func0 -> SetParameters(a, h1->GetMean(), h1->GetStdDev());
    cout << "The initial parameters are: " << endl;
    for (int ip = 0; ip < func0 -> GetNpar(); ++ ip) {
      cout << ip << " " << func0->GetParName(ip) << " " << func0->GetParameter(ip) << endl;
    }
    if (isProtonFit and fabs(p) > 950.) {
      func0->SetParameter(2, 39.5);
      h1 -> Fit("func0", "", "", 38, 40);
    }
    else
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
    func -> SetNpx(1000);
    double A = h3 -> GetMaximum();//func0 -> GetParameter(0); //h1 -> GetMaximum();
    cout << "A=" << A << endl;
    func->SetParName(0, "norm1");
    func->SetParameter(0, A);
    if (!isThreeComponentFit)
      func->SetParLimits(0, 0.6*A, 1.2*A);
    if (isProtonFit) {
      func->SetParLimits(0, 0, 1.2*A);
    }
    
    func->SetParName(1, "norm2");
    func->SetParameter(1, 0.3*A); // h3 -> GetMaximum()/2.);
    if (!isThreeComponentFit)
      func->SetParLimits(1, 0, 0.6*A);
    if (isProtonFit) {
      func->SetParameter(1, A); 
      func->SetParLimits(1, 0, 1.2*A);
      
      if (fabs(p) > 950) {
	func->SetParameter(0, A/3.); 
	func->SetParLimits(0, 0.1*A, 0.5*A);
	func->SetParameter(1, A);
	func->SetParLimits(1, 0.8*A, 2*A);	
      }
    }

    double delta = 0.05;    
    func->SetParName(2, "mean_e");
    func->SetParameter(2, func0 -> GetParameter(1));//h1->GetMean());
    //    func->SetParLimits(2, (1.-delta)*h1->GetMean(), (1. + delta)*h1->GetMean());
    if (isProtonFit and fabs(p) > 950) {
      func->SetParameter(2, 39.5);
    }

    func->SetParName(3, "sigma_e");
    //func->FixParameter(2, 0.215);
    func->SetParameter(3, func0 -> GetParameter(2));//h1->GetStdDev());
    //    func->SetParLimits(3, 0.4*h1->GetStdDev(), 1.2*h1->GetStdDev()); 
    
    func->SetParName(4, "norm3");
    func->SetParameter(4, 0.1*A); //h3 -> GetMaximum()/3.);
    func->SetParLimits(4, 0, 0.5*A);
    
    double mdelta = 0.15;
    double dtval_mu = xe + dtmu; // func0 -> GetParameter(1) + dtmu ;// h3 ->GetMean() +l/c/GetBeta(m[1], p)-l/c/GetBeta(m[0], p);
    func->SetParName(5, "mean_mu");
    if (!isProtonFit) {
      func->SetParameter(5, dtval_mu);
    //func->FixParameter(5, dtval_mu);
      func->SetParLimits(5, (1. - mdelta)*dtval_mu, (1. + mdelta)*dtval_mu);       
    } else {
      func->SetParName(5, "mean_p");
      double valp = xe + dtp;
      func->SetParameter(5, valp);
      //func->FixParameter(5, valp);
      double pdelta = 0.15;
      func->SetParLimits(5, (1. - pdelta)*valp, (1. + pdelta)*valp);
      if (fabs(p) > 950) {
	double valp = 39.5 + dtp;
	func->SetParameter(5, valp);
	func->SetParLimits(5, (1. - pdelta)*valp, (1. + pdelta)*valp);
      }
    }

    
    /*func->SetParName(6, "sigma_mu");
    func->SetParameter(6, h3->GetStdDev());
    func->SetParLimits(6, 0.9*h3->GetStdDev(), 1.1*h3->GetStdDev()); */
    
    double dtval_pi = func0 -> GetParameter(1) + dtpi; // h3->GetMean()  l/c/GetBeta(m[2], p)-l/c/GetBeta(m[0], p);
    func->SetParName(6, "mean_pi");
    func->SetParameter(6, dtval_pi);
    func->FixParameter(6, dtval_pi);
    func->SetParLimits(6,  (1. - mdelta)*dtval_pi,  (1. - mdelta)*dtval_pi);       
    
    func->SetParName(7, "sigma_mupi");
    //func->FixParameter(2, 0.215);
    func->SetParameter(7, h3->GetStdDev());
    //    func->SetParLimits(7, 0.5*h3->GetStdDev(), 1.3*h3->GetStdDev()); 
    if (isProtonFit) {
     func->SetParameter(7, 0.5);
     func->SetParName(7, "sigma_p");
    }
    
    // special cases

  // 300n
    if (fabs(p + 300) < 1.e-3) {
	func->SetParameter(0, A); 
	func->SetParLimits(0, 0.1*A, 2*A);
	func->SetParameter(1, A);
	func->SetParLimits(1, 0.*A, 2*A);
	func->SetParameter(5, 0.5*( dtval_mu + dtval_pi) );
	func->SetParameter(7, 0.9);
      }


    
    // 340n
    if (fabs(p + 340) < 1.e-3) {
	func->SetParameter(0, 1.5*A); 
	func->SetParLimits(0, 0.1*A, 1.5*A);
	func->SetParameter(1, A/4);
	func->SetParLimits(1, 0.*A, 1.3*A);
	func->SetParameter(6, 1.1*dtval_mu);
      }

    // 340p
    if (fabs(p - 340) < 1.e-3) {
	func->SetParameter(0, 1.*A); 
	func->SetParLimits(0, 0.1*A, 1.5*A);
	func->SetParameter(1, A/4);
	func->SetParLimits(1, 0.*A, 1.*A);
	func->SetParameter(6, 1.1*dtval_mu);
      }

    // 360n
    if (fabs(p + 360) < 1.e-3) {
	func->SetParameter(0, 1.*A); 
	func->SetParLimits(0, 0.1*A, 1.5*A);
	func->SetParameter(1, A/4);
	func->SetParLimits(1, 0.*A, 1.*A);
	func->SetParameter(6, 1.1*dtval_mu);
      }


   // 360p
    if (fabs(p - 360) < 1.e-3) {
	func->SetParameter(0, 1.*A); 
	func->SetParLimits(0, 0.1*A, 1.5*A);
	func->SetParameter(1, A/4);
	func->SetParLimits(1, 0.*A, 1.*A);
	func->SetParameter(6, 1.1*dtval_mu);
      }


  

    cout << "The initial parameters are: " << endl;
    for (int ip = 0; ip < func -> GetNpar(); ++ ip) {
      cout << ip << " " << func->GetParName(ip) << " " << func->GetParameter(ip) << endl;      
    }

    if (!isThreeComponentFit) {
      func -> FixParameter(4, 0.);
      func -> FixParameter(6, 0.);
    }
    h3->Fit(func);
    
    TCanvas *cv1 = new TCanvas("cv1", "", 1200, 800);
    cv1->cd()->SetTickx(kTRUE);
    cv1->cd()->SetTicky(kTRUE);
    if (logy) {
      cv1->cd()->SetLogy(kTRUE);
      h3->SetMaximum(20.*h3->GetMaximum());
    }
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

    if (isThreeComponentFit) {
      func3->SetParameter(0, func->GetParameter(4));
      func3->SetParameter(1, func->GetParameter(6));
      func3->SetParameter(2, func->GetParameter(7));
      func3->SetLineWidth(2);
      func3->SetLineColor(kGreen+1);
    } 
    
    func1->SetLineWidth(2);
    func1->SetLineColor(kRed+1);
    
    func2->SetLineWidth(2);
    func2->SetLineColor(kBlue+1);
    
    if (!isThreeComponentFit) {
      func2->SetLineColor(kBlack);
      func2->SetLineStyle(2);
    }

    func1->SetLineStyle(2);
    func3->SetLineStyle(2);

    
    TLegend *leg1 = new TLegend(0.60, 0.6, 0.9, 0.88);
    leg1 -> SetBorderSize(0);
    int ip = int(fabs(p));
    leg1 -> SetHeader("WCTE TB July 2022");
    
    int Ne = func->GetParameter(0)/h1->GetBinWidth(1);
    int Nmu = func->GetParameter(1)/h1->GetBinWidth(1);
    int Npi = 0.;
    int NeErr = func->GetParError(0)/h1->GetBinWidth(1);
    int NmuErr = func->GetParError(1)/h1->GetBinWidth(1);
    int NpiErr = 0.;
    if (isThreeComponentFit) {
      func->GetParameter(4)/h1->GetBinWidth(1);
      int NpiErr = func->GetParError(4)/h1->GetBinWidth(1);
    }
        
    leg1->AddEntry(h3,"Data", "PE");
    leg1->AddEntry(func, "Fit", "L");

    string ssignum = "+";
    if (p < 0)
      ssignum = "-";
    title = "Fit e^{" + ssignum + "}, N = " + to_string(Ne) + " #pm " + to_string(NeErr);
    leg1->AddEntry(func1, title.c_str(), "L");
  
    
    if (isThreeComponentFit) {
      title = "Fit #mu^{" + ssignum + "}, N = " + to_string(Nmu) + " #pm " + to_string(NmuErr);
      title = "Fit #pi^{" + ssignum + "}, N = " + to_string(Npi) + " #pm " + to_string(NpiErr);
      leg1->AddEntry(func3, title.c_str(), "L");
    } else {
      if (!isProtonFit)
	title = "Fit #mu^{" + ssignum + "}+#pi^{" + ssignum + "}, N = " + to_string(Nmu) + " #pm " + to_string(NmuErr);
      else	
	title = "Fit p^{" + ssignum + "}, N = " + to_string(Nmu) + " #pm " + to_string(NmuErr);
    }
    leg1->AddEntry(func2, title.c_str(), "L");
    
    func1->Draw("Csame");
    func2->Draw("Csame");          
    if (isThreeComponentFit)
      func3->Draw("Csame");
    
    leg1->Draw("same");
 

    // TLines
    int ls = 2;
    int lw = 2;
    double y0 = h3->GetMinimum();
    double y1 = 1.1* h3->GetBinContent(imax);
    double dy = 0.02 * h3->GetBinContent(imax);
    double dx = -0.005 * ( h3->GetXaxis()->GetXmax() - h3->GetXaxis()->GetXmin() );
    xe = func -> GetParameter(2);
    
    TLine *le = new TLine(xe, y0, xe, y1);
    le -> SetLineColor(kRed+1);
    le -> SetLineStyle(ls);
    le -> SetLineWidth(lw);
    le -> Draw();
    TLatex *char_e = new TLatex(xe + dx, y1 + dy, "e");
    char_e -> SetTextSize(0.03);
    char_e -> SetTextColor(kRed+1);
    char_e -> Draw();

    TLine *lmu = new TLine(xe + dtmu, y0, xe + dtmu, y1);
    lmu -> SetLineColor(kBlue+1);
    lmu -> SetLineStyle(ls);
    lmu -> SetLineWidth(lw);
    lmu -> Draw();
    TLatex *char_mu = new TLatex(xe + dx + dtmu, y1 + dy, "#mu");
    char_mu -> SetTextSize(0.03);
    char_mu -> SetTextColor(kBlue+1);
    char_mu -> Draw();

    TLine *lpi = new TLine(xe + dtpi, y0, xe + dtpi, y1);
    lpi -> SetLineColor(kGreen+1);
    lpi -> SetLineStyle(ls);
    lpi -> SetLineWidth(lw);
    lpi -> Draw();
    TLatex *char_pi = new TLatex(xe + dx + dtpi, y1 + dy, "#pi");
    char_pi -> SetTextSize(0.03);
    char_pi -> SetTextColor(kGreen+1);
    char_pi -> Draw();

    //double dtpp = GetTofDiffWrtE(m[3], p, m[0]);
    TLine *lp = new TLine(xe + dtp, y0, xe + dtp, y1);
    lp -> SetLineColor(kBlack);
    lp -> SetLineStyle(ls);
    lp -> SetLineWidth(lw);
    lp -> Draw();
    TLatex *char_p = new TLatex(xe + dx + dtp, y1 + dy, "p");
    char_p -> SetTextSize(0.03);
    char_p -> SetTextColor(kBlack);
    char_p -> Draw();
    
    TLine *ld = new TLine(xe + dtd, y0, xe + dtd, y1);
    ld -> SetLineColor(kBlack);
    ld -> SetLineStyle(ls);
    ld -> SetLineWidth(lw);
    ld -> Draw();
    TLatex *char_d = new TLatex(xe + dx + dtd, y1 + dy, "d");
    char_d -> SetTextSize(0.03);
    char_d -> SetTextColor(kBlack);
    char_d -> Draw();
    
    TLine *lt = new TLine(xe + dtt, y0, xe + dtt, y1);
    lt -> SetLineColor(kBlack);
    lt -> SetLineStyle(ls);
    lt -> SetLineWidth(lw);
    lt -> Draw();
    TLatex *char_t = new TLatex(xe + dx + dtt, y1 + dy, "t");
    char_t -> SetTextSize(0.03);
    char_t -> SetTextColor(kBlack);
    char_t -> Draw();
    
    TLine *la = new TLine(xe + dta, y0, xe + dta, y1);
    la -> SetLineColor(kBlack);
    la -> SetLineStyle(ls);
    la -> SetLineWidth(lw);
    la -> Draw();
    TLatex *char_a = new TLatex(xe + dx + dta, y1 + dy, "#alpha");
    char_a -> SetTextSize(0.03);
    char_a -> SetTextColor(kBlack);
    char_a -> Draw();
    

    double chi2 = func0 -> GetChisquare();
    int ndf = func0 -> GetNDF();
    TLatex *chi2text = new TLatex(0.13, 0.85, Form("#chi^{2}/ndf = %2.1f/%i = %1.1f", chi2, ndf, chi2/ndf));
    chi2text -> SetNDC();
    chi2text -> SetTextSize(0.03);
    chi2text -> Draw();
     
    cout << GetBeta(m[0], p) << " " << GetBeta(m[1], p) << " " << GetBeta(m[2], p) << " " << GetBeta(m[3], p) << endl;
    cout << l/c/GetBeta(m[0], p) << " " << l/c/GetBeta(m[1], p) << " " << l/c/GetBeta(m[2], p) << " " << l/c/GetBeta(m[3], p) << endl;
    cout << h2->GetMean() << " " << h2->GetMean()+l/c/GetBeta(m[1], p)-l/c/GetBeta(m[0], p) << " " << h2->GetMean()+l/c/GetBeta(m[2], p)-l/c/GetBeta(m[0], p) << endl;
        
    cout << h1->GetEntries() << " " << h2->GetEntries() << " " << h3->GetEntries() << endl;
    
    cout << "N_e = " <<  func->GetParameter(0)/h1->GetBinWidth(1) << endl;
    if (isThreeComponentFit) {
      cout << "N_mu = " <<  func->GetParameter(1)/h1->GetBinWidth(1) << endl;
      cout << "N_pi = " <<  func->GetParameter(4)/h1->GetBinWidth(1) << endl;
    } else {
      cout << "N_mupi = " <<  func->GetParameter(1)/h1->GetBinWidth(1) << endl;
    }
    
    TString asciiname = "ascii_" + TString(inFile->GetName()).ReplaceAll(".root",".txt");
    ofstream *asciifile = new ofstream(asciiname.Data());
    (*asciifile) << "mom " << p << endl;
    (*asciifile) << "N_e " << func->GetParameter(0)/h1->GetBinWidth(1) << " " << NeErr << endl;
    if (isProtonFit) {
	(*asciifile) << "N_proton " <<  func->GetParameter(1)/h1->GetBinWidth(1) << " " << NmuErr << endl;
    } else {

      double tof_reso = 0.300; // ns	
      TH2D* hACT2_act2cut = (TH2D*) inFile -> Get("hRef_TOFACT2V_act2cut");
      TH2D* hACT3_act3cut = (TH2D*) inFile -> Get("hRef_TOFACT3V_act3cut");
      TH2D* hACT2_all = (TH2D*) inFile -> Get("hRef_TOFACT2V");
      TH2D* hACT3_all = (TH2D*) inFile -> Get("hRef_TOFACT3V");
      hACT2_act2cut -> SetStats(0);
      hACT2_all -> SetStats(0);
      hACT3_act3cut -> SetStats(0);
      hACT3_all -> SetStats(0);
      
      if (isThreeComponentFit) {
	(*asciifile) << "N_mu " <<  func->GetParameter(1)/h1->GetBinWidth(1) << " " << NmuErr << endl;
	(*asciifile) << "N_pi " <<  func->GetParameter(4)/h1->GetBinWidth(1) << " " << NpiErr << endl;
	// integrate the ASC2 distribution below non-ele thr in the TOF time widow to get the ASC muon and pion ID efficiency

	// mu

	double nACT2_act2cut_mu = IntegrateAlongYInGivenXWindow(hACT2_act2cut, func->GetParameter(5), tof_reso);
	double nACT2_all_mu    = IntegrateAlongYInGivenXWindow(hACT2_all, func->GetParameter(5), tof_reso);
	double eff_mu          = nACT2_act2cut_mu / nACT2_all_mu;
  cout << " nACT2_act2cut_mu=" << nACT2_act2cut_mu
       << " nACT2_all_mu=" << nACT2_all_mu
       << endl;

	// pi
	double pioff = 0.35;
	double nACT3_act3cut_pi = IntegrateAlongYInGivenXWindow(hACT3_act3cut, func->GetParameter(6) + pioff, tof_reso);
	double nACT3_all_pi    = IntegrateAlongYInGivenXWindow(hACT3_all, func->GetParameter(6) + pioff, tof_reso);
	double eff_pi          = nACT3_act3cut_pi / nACT3_all_pi;
cout << " nACT3_act3cut_pi=" << nACT3_act3cut_pi
     << " nACT3_all_pi=" << nACT3_all_pi
     << endl;

	TString canname = asciiname.ReplaceAll("ascii", "ACTCuts").ReplaceAll(".txt", "").ReplaceAll("_output", "");
	TCanvas *h2can = new TCanvas(canname, canname, 200, 200, 1000, 500);
	h2can -> Divide(2,2);

	gStyle->SetPalette(kSolar);
	
	h2can -> cd(1);
	hACT2_all -> Draw("colz");
	h2can -> cd(2);
	hACT2_act2cut -> Draw("colz");
	double x1 = func->GetParameter(5) - tof_reso;
 	double x2 = func->GetParameter(5) + tof_reso;
	double y1 = 	hACT2_act2cut -> GetYaxis() -> GetXmax();
	double y2 = 	hACT2_act2cut -> GetYaxis() -> GetXmin();

	h2can -> cd(3);
	hACT3_all -> Draw("colz");
	h2can -> cd(4);
	hACT3_act3cut -> Draw("colz");
	x1 = func->GetParameter(5) - tof_reso;
 	x2 = func->GetParameter(5) + tof_reso;
	y1 = 	hACT3_act3cut -> GetYaxis() -> GetXmax();
	y2 = 	hACT3_act3cut -> GetYaxis() -> GetXmin();

	h2can -> cd(2);
	
	TLine *xlmu1 = new TLine(x1, y1, x1, y2);
	xlmu1 -> SetLineColor(kBlue+1);
	xlmu1 -> SetLineStyle(ls);
	xlmu1 -> SetLineWidth(lw);
	xlmu1 -> Draw();

	TLine *xlmu2 = new TLine(x2, y1, x2, y2);
	xlmu2 -> SetLineColor(kBlue+1);
	xlmu2 -> SetLineStyle(ls);
	xlmu2 -> SetLineWidth(lw);
	xlmu2 -> Draw();


	x1 = pioff+ func->GetParameter(6) - tof_reso;
 	x2 = pioff + func->GetParameter(6) + tof_reso;
	
	TLine *xlpi1 = new TLine(x1, y1, x1, y2);
	xlpi1 -> SetLineColor(kGreen+1);
	xlpi1 -> SetLineStyle(ls);
	xlpi1 -> SetLineWidth(lw);
	xlpi1 -> Draw();

	TLine *xlpi2 = new TLine(x2, y1, x2, y2);
	xlpi2 -> SetLineColor(kGreen+1);
	xlpi2 -> SetLineStyle(ls);
	xlpi2 -> SetLineWidth(lw);
	xlpi2 -> Draw();

	h2can -> Print(canname + ".png");
	h2can -> Print(canname + ".pdf");

	  
	// so far zero error on the efficiency...
	(*asciifile) << "eff_mu " << eff_mu  << " " << 0. << endl;
	(*asciifile) << "eff_pi " << eff_pi  << " " << 0. << endl;
	
	
      } else     {
	(*asciifile) << "N_mupi " <<  func->GetParameter(1)/h1->GetBinWidth(1) << " " << NmuErr << endl;
      
      	// mu+pi
	double nACT2_act2cut_mupi = IntegrateAlongYInGivenXWindow(hACT2_act2cut, func->GetParameter(5), tof_reso);
	double nACT2_all_mupi    = IntegrateAlongYInGivenXWindow(hACT2_all, func->GetParameter(5), tof_reso);
	double eff_mupi          = nACT2_act2cut_mupi / nACT2_all_mupi;
	// so far zero error on the efficiency...
	(*asciifile) << "eff_mupi " << eff_mupi  << " " << 0. << endl;

      }
    }
    if (asciifile) asciifile->close();

    string name = fileName.substr(0, fileName.size()-5) + "_TOFfit.png";
    string namepdf = fileName.substr(0, fileName.size()-5) + "_TOFfit.pdf";
    cv1->Print(name.c_str());
    cv1->Print(namepdf.c_str());    

    
    cout << "DONE!" << endl;
    cout << "See also " << asciiname.Data() << endl;
      
    
}


// ______________________________________________________
// ______________________________________________________
// ______________________________________________________

