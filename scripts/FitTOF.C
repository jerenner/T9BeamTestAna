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

void SetInitialFitPars(bool isThreeComponentFit, bool isProtonFit,
		       TH1D *h3, 
		       int p,
		       double A,
		       TF1 *func0,
		       TF1 *func,
		       double tof_el,
		       double delta, double mdelta,
		       double dtmu, double dtp,
		       double tof_mu, double tof_pi)
{

  // electrons:  norm [0], mean [2], sigma [3]
  // muons:      norm [1], mean [5], sigma [7]
  // pions:      norm [4], mean [6], sigma [7] -- same!
  
    func->SetParName(0, "norm1");
    func->SetParameter(0, A/2);
    if (!isThreeComponentFit)
      func->SetParLimits(0, 0.6*A, 1.2*A);
    if (isProtonFit) {
      func->SetParLimits(0, 0, 1.2*A);
    }
    
    func->SetParName(1, "norm_mu");
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

   
    func->SetParName(2, "mean_e");
    func->SetParameter(2, func0 -> GetParameter(1));//h1->GetMean());
    //    func->SetParLimits(2, (1.-delta)*h1->GetMean(), (1. + delta)*h1->GetMean());
    //if (isProtonFit && fabs(p) > 950) {
    func->SetParameter(2, 39.5);
    //}

    func->SetParName(3, "sigma_e");
    //func->FixParameter(2, 0.215);
    func->SetParameter(3, func0 -> GetParameter(2));//h1->GetStdDev());
    //    func->SetParLimits(3, 0.4*h1->GetStdDev(), 1.2*h1->GetStdDev()); 
    
    func->SetParName(4, "norm_pi");
    func->SetParameter(4, 0.1*A); //h3 -> GetMaximum()/3.);
    func->SetParLimits(4, 0, 1.*A); // was: 0.1. but not enough for high momenta

    func->SetParName(5, "mean_mu");

    if (!isProtonFit) {
      func->SetParameter(5, tof_mu);
      //func->FixParameter(5, tof_mu);
      //func->SetParLimits(5, (1. - mdelta)*tof_mu, (1. + mdelta)*tof_mu);
      func->SetParLimits(5, tof_el + dtmu*(1. - mdelta),       tof_el + dtmu*(1. + mdelta) );
      cout << "Muons mean par limits: " << tof_el + dtmu*(1. - mdelta) << ", " << tof_el + dtmu*(1. + mdelta)  << endl;
      if (!isThreeComponentFit) {
	func->SetParName(1, "normmupi");
      }
      // jk 28.11.2022
      func->SetParLimits(7, 0.05, 0.4); 
    } else {
      // proton fit
      func->SetParName(5, "mean_p");
      double valp = tof_el + dtp;
      func->SetParameter(5, valp);
      //func->FixParameter(5, valp);
      double pdelta = 0.15;
      func->SetParLimits(5, (1. - pdelta)*valp, (1. + pdelta)*valp);
      if (fabs(p) > 950) {
	double valp = 39.5 + dtp;
	func->SetParameter(5, valp);
	func->SetParLimits(5, (1. - pdelta)*valp, (1. + pdelta)*valp);
      }
    } // proton fit
    
    /*func->SetParName(6, "sigma_mu");
    func->SetParameter(6, h3->GetStdDev());
    func->SetParLimits(6, 0.9*h3->GetStdDev(), 1.1*h3->GetStdDev()); */
    
    func->SetParName(6, "mean_pi");
    func->SetParameter(6, tof_pi);
    //func->FixParameter(6, tof_pi);
    func->SetParLimits(6,  (1. - mdelta)*tof_pi,  (1. - mdelta)*tof_pi);       
    
    func->SetParName(7, "sigma_mupi");
    //func->FixParameter(2, 0.215);
    func->SetParameter(7, h3->GetStdDev());
    // jk 28.11.2022
    //func->SetParLimits(7, 0.05, 0.6); 
    //func->SetParLimits(7, 0.5*h3->GetStdDev(), 1.3*h3->GetStdDev());


    if (!isProtonFit) { 
      // adjust for the momentum bias:
      // Dec 2022
      if (p < 0) {
	// negative bias fo the negative beam
	if (fabs(p) < 255) {
	  cout << "OK, adjusting beam biases for negative low momenta..." << endl;
	  double bias_mu = 0.5 * 200. / fabs(p);
	  double bias_pi = 0.5 * 200. / fabs(p);
	  func->SetParameter(5, func -> GetParameter(5) + bias_mu);
	  func->SetParameter(6, func -> GetParameter(6) + bias_pi);
	} else {
	}
      } else  {
	// positive bias fo the positive beam
	if (fabs(p) < 255) {
	  cout << "OK, adjusting beam biases for positive low momenta..." << endl;
	  double bias_mu = -0.2 * 200. / fabs(p);
	  double bias_pi = -0.25 * 200. / fabs(p);
	  func->SetParameter(5, func -> GetParameter(5) + bias_mu);
	  func->SetParameter(6, func -> GetParameter(6) + bias_pi);
	} else {
	}
      }

    
      if (fabs(p) > 255) {
	// pions expected to have a large contribution
	double sf = 1.;
	cout << h3 -> GetName() << endl;
	if (TString(h3->GetName()).Contains("act2cut") || TString(h3->GetName()).Contains("act3cut"))
	  sf = 1.9;
	cout << "OK, adjusting mu and pi parameters for higher momenta... p=" << p << " sf=" << sf << endl;
	// pions expected to have a large contribution
	func->SetParameter(1, sf*0.5*A);
	func->SetParameter(4, sf*0.5*A);
	func->SetParameter(7, 0.15);
      }
    }
    else {
      // just a name change;-)
      func->SetParameter(7, 0.5);
      func->SetParName(7, "sigma_p");
    }
    
}


// ______________________________________________________

void GetIndividualFitComponents(TF1 *func, TF1 *&func1, TF1 *&func2, TF1 *&func3, bool isThreeComponentFit, TString tag = "")
{
  func1 = new TF1("func1" + tag, "[0]*TMath::Gaus(x, [1], [2], 1)", 38, 43.5);
  func2 = new TF1("func2" + tag, "[0]*TMath::Gaus(x, [1], [2], 1)", 38, 43.5);
  func3 = new TF1("func3" + tag, "[0]*TMath::Gaus(x, [1], [2], 1)", 38, 43.5);
  
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

  return;
}

// ______________________________________________________

void DrawExpectedTimeArrivalLinesLegAndChi2(TH1D *h1,
					    TH1D *h3,
					    TF1* func,
					    TF1* func1, TF1* func2, TF1* func3,
					    int p, int imax, 
					    bool isThreeComponentFit, bool isProtonFit,
					    int Ne, int Nmu, int Npi,
					    int NeErr, int NmuErr, int NpiErr,
					    double tof_el, double dtmu, double dtpi, double dtp, double dtd, double dtt, double dta)
{
    string ptag = " (Pos)";
    if (p < 0.)
      ptag = " (Neg)";
    
    string title = "p = " + to_string(abs(p)) + " MeV/c" + ptag ;
    if (TString(h3->GetName()).Contains("act2cut"))
      title = title + ", ACT2 cuts";
    if (TString(h3->GetName()).Contains("act3cut"))
      title = title + ", ACT3 cuts";
    h3->SetTitle(title.c_str());
    h3->GetXaxis()->SetTitle("t_{1}-t_{0} [ns]");
    h3->GetYaxis()->SetTitle("Events");
    h3->SetStats(kFALSE);
  
   // TLines
    int ls = 2;
    int lw = 2;
    
    double y0 = h3->GetMinimum();
    double y1 = 1.1* h3->GetBinContent(imax);
    double dy = 0.02 * h3->GetBinContent(imax);
    double dx = -0.005 * ( h3->GetXaxis()->GetXmax() - h3->GetXaxis()->GetXmin() );
    
    TLine *le = new TLine(tof_el, y0, tof_el, y1);
    le -> SetLineColor(kRed+1);
    le -> SetLineStyle(ls);
    le -> SetLineWidth(lw);
    le -> Draw();
    TLatex *char_e = new TLatex(tof_el + dx, y1 + dy, "e");
    char_e -> SetTextSize(0.03);
    char_e -> SetTextColor(kRed+1);
    char_e -> Draw();

    TLine *lmu = new TLine(tof_el + dtmu, y0, tof_el + dtmu, y1);
    lmu -> SetLineColor(kBlue+1);
    lmu -> SetLineStyle(ls);
    lmu -> SetLineWidth(lw);
    lmu -> Draw();
    TLatex *char_mu = new TLatex(tof_el + dx + dtmu, y1 + dy, "#mu");
    char_mu -> SetTextSize(0.03);
    char_mu -> SetTextColor(kBlue+1);
    char_mu -> Draw();

    TLine *lpi = new TLine(tof_el + dtpi, y0, tof_el + dtpi, y1);
    lpi -> SetLineColor(kGreen+1);
    lpi -> SetLineStyle(ls);
    lpi -> SetLineWidth(lw);
    lpi -> Draw();
    TLatex *char_pi = new TLatex(tof_el + dx + dtpi, y1 + dy, "#pi");
    char_pi -> SetTextSize(0.03);
    char_pi -> SetTextColor(kGreen+1);
    char_pi -> Draw();

    //double dtpp = GetTofDiffWrtE(m[3], p, m[0]);
    TLine *lp = new TLine(tof_el + dtp, y0, tof_el + dtp, y1);
    lp -> SetLineColor(kBlack);
    lp -> SetLineStyle(ls);
    lp -> SetLineWidth(lw);
    lp -> Draw();
    TLatex *char_p = new TLatex(tof_el + dx + dtp, y1 + dy, "p");
    char_p -> SetTextSize(0.03);
    char_p -> SetTextColor(kBlack);
    char_p -> Draw();
    
    TLine *ld = new TLine(tof_el + dtd, y0, tof_el + dtd, y1);
    ld -> SetLineColor(kBlack);
    ld -> SetLineStyle(ls);
    ld -> SetLineWidth(lw);
    ld -> Draw();
    TLatex *char_d = new TLatex(tof_el + dx + dtd, y1 + dy, "d");
    char_d -> SetTextSize(0.03);
    char_d -> SetTextColor(kBlack);
    char_d -> Draw();
    
    TLine *lt = new TLine(tof_el + dtt, y0, tof_el + dtt, y1);
    lt -> SetLineColor(kBlack);
    lt -> SetLineStyle(ls);
    lt -> SetLineWidth(lw);
    lt -> Draw();
    TLatex *char_t = new TLatex(tof_el + dx + dtt, y1 + dy, "t");
    char_t -> SetTextSize(0.03);
    char_t -> SetTextColor(kBlack);
    char_t -> Draw();
    
    TLine *la = new TLine(tof_el + dta, y0, tof_el + dta, y1);
    la -> SetLineColor(kBlack);
    la -> SetLineStyle(ls);
    la -> SetLineWidth(lw);
    la -> Draw();
    TLatex *char_a = new TLatex(tof_el + dx + dta, y1 + dy, "#alpha");
    char_a -> SetTextSize(0.03);
    char_a -> SetTextColor(kBlack);
    char_a -> Draw();

  TLegend *leg1 = new TLegend(0.60, 0.6, 0.9, 0.88);
    leg1 -> SetBorderSize(0);
    int ip = int(fabs(p));
    leg1 -> SetHeader("WCTE TB July 2022");

    leg1->AddEntry(h3,"Data", "PE");
    leg1->AddEntry(func, "Fit", "L");

    string ssignum = "+";
    if (p < 0)
      ssignum = "-";
    title = "Fit e^{" + ssignum + "}, N = " + to_string(Ne) + " #pm " + to_string(NeErr);
    leg1->AddEntry(func1, title.c_str(), "L");
    
    if (isThreeComponentFit) {
      title = "Fit #mu^{" + ssignum + "}, N = " + to_string(Nmu) + " #pm " + to_string(NmuErr);
      leg1->AddEntry(func2, title.c_str(), "L");
      title = "Fit #pi^{" + ssignum + "}, N = " + to_string(Npi) + " #pm " + to_string(NpiErr);
      leg1->AddEntry(func3, title.c_str(), "L");
    } else {
      if (!isProtonFit)
	title = "Fit #mu^{" + ssignum + "}+#pi^{" + ssignum + "}, N = " + to_string(Nmu) + " #pm " + to_string(NmuErr); 
      else	
	title = "Fit p^{" + ssignum + "}, N = " + to_string(Nmu) + " #pm " + to_string(NmuErr);
    }
    leg1->Draw("same");

    double chi2 = func -> GetChisquare(); // was func0!
    int ndf = func -> GetNDF(); // was func0!
    TLatex *chi2text = new TLatex(0.13, 0.85, Form("#chi^{2}/ndf = %2.1f/%i = %1.1f", chi2, ndf, chi2/ndf));
    chi2text -> SetNDC();
    chi2text -> SetTextSize(0.03);
    chi2text -> Draw();

    return;
}

// ______________________________________________________

void AssignFitResults(TF1 *func, TH1D *h1, int &Ne, int &Nmu, int &Npi, 
                      int &NeErr, int &NmuErr, int &NpiErr, bool isThreeComponentFit)
{
 Ne  = (int) func->GetParameter(0) / h1->GetBinWidth(1);
 Nmu = (int) func->GetParameter(1) / h1->GetBinWidth(1);
 Npi = 0;
 NeErr  = (int) func->GetParError(0) / h1->GetBinWidth(1);
 NmuErr = (int) func->GetParError(1) / h1->GetBinWidth(1);
 NpiErr = 0;
 if (isThreeComponentFit) {
      Npi = func->GetParameter(4) / h1->GetBinWidth(1);
      NpiErr = (int) func->GetParError(4) / h1->GetBinWidth(1);
  }
 return;   
}

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

void PrintIntegrals(TH1D *h1, TH1D *h2, TH1D *h3)
{
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
    
    // time diffs of mu and pi to e
    //for (int ibin = 0; ibin < h1->GetNbinsX()+1; ++ibin) {
    //  cout << "hall i=" << ibin << " y=" << h1->GetBinContent(ibin) << endl;
    //}

    
}

// ______________________________________________________
// ______________________________________________________
// ______________________________________________________


void FitTOF(string fileName, int p, bool logy = true) {

    // masses: [MeV]
    //             e, mu, pi, proton, deuteron, tritium, alpha
    double m[7] = {0.511, 105.66, 139.57, 938.470, 1876., 3.01604928*931.494102, 3727.379};
    double c = 0.299792458; // [m/ns]
    double l = 2.9; // [m]

    bool isThreeComponentFit = fabs(p) < 301.;
    bool isProtonFit = fabs(p) > 390.;

    
    cout << "Configured as" << endl;
    cout << " isThreeComponentFit: " << isThreeComponentFit << endl 
	 << " isProtonFit: " << isProtonFit << endl 
	 << endl;
    
    // dt = delta_t = delay times of mu, pi, protons w.r.t. electrons
    double dtmu = GetTofDiffWrtE(m[1], p, m[0]);
    double dtpi = GetTofDiffWrtE(m[2], p, m[0]);
    double dtp = GetTofDiffWrtE(m[3], p, m[0]);
    // for fun also deuteron, triton, alpha
    double dtd = GetTofDiffWrtE(m[4], p, m[0]);
    double dtt = GetTofDiffWrtE(m[5], p, m[0]);
    double dta = GetTofDiffWrtE(m[6], p, m[0]);

    // HARD-CODED TOF BOUNDRIES!
    double tmin = 38.;
    double tmax = 43.5;
    if (isProtonFit) {
      tmax = 42. + dtp;
    }
    cout << "Will fit TOF distribution betweem " << tmin << " " << tmax << endl;
   
    // ORIGINAL: string basehname = "hTOFOther"; original, for standard TOF fit
    // hack by Jiri, to extract mu and pi efficiencies ont he full sample and after dedicated ACT2 cuts
    string basehname = "hTOFAll"; // still, this is after custom ACT1 electron removal cut
    if (isProtonFit) {
      basehname = "hTOFAllWide"; 
    }

    // read histograms
    TFile *inFile = new TFile(fileName.c_str(), "READ");
    // TODO
    // carefully: to rename h1 to hTOFAll, h2 to hTOFEl, h3 to hTOF!!
    TH1D* h1 = (TH1D*) inFile->Get("hTOFAll");
    TH1D* h2 = (TH1D*) inFile->Get("hTOFEl");
    TH1D* h3 = (TH1D*) inFile->Get(basehname.c_str());

    // TOF distr. after ACT cuts:
    TH1D* h3_act2cut = (TH1D*) inFile->Get("hTOF_act2cut");
    TH1D* h3_act3cut = (TH1D*) inFile->Get("hTOF_act3cut");

    int imax = h1->GetMaximumBin();
    double tof_el = h1->GetBinCenter(imax);
    
    PrintIntegrals(h1, h2, h3);

    // pre-fit to get initial pars for the full fit
    TCanvas *can_aux = new TCanvas("can_aux", "", 800, 800);
    can_aux -> Divide(2,2);
    can_aux->cd(1);
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
    can_aux -> cd(2);
    h2 -> Draw("e1");
    can_aux -> cd(3);
    h3 -> Draw("e1");
    can_aux -> cd(4);

    // the main full fit:
    TF1* func = new TF1("func", "[0]*TMath::Gaus(x, [2], [3], 1) + [1]*TMath::Gaus(x, [5], [7], 1) + [4]*TMath::Gaus(x, [6], [7], 1)", 38, 43.5);
    func -> SetNpx(1000);
    func->SetLineWidth(2);
    func->SetLineColor(kBlack);

    double A = h3 -> GetMaximum();//func0 -> GetParameter(0); //h1 -> GetMaximum();
    cout << "A=" << A << endl;

    // for fit parameters initial values and limits:
    double delta = 0.05;    
    double mdelta = 0.8;
    double tof_mu = tof_el + dtmu; // func0 -> GetParameter(1) + dtmu ;// h3 ->GetMean() +l/c/GetBeta(m[1], p)-l/c/GetBeta(m[0], p);
    double tof_pi = func0 -> GetParameter(1) + dtpi; // h3->GetMean()  l/c/GetBeta(m[2], p)-l/c/GetBeta(m[0], p);
    
    // set initial fit parameters
    SetInitialFitPars(isThreeComponentFit, isProtonFit, h3, p, A, func0, func, tof_el, delta, mdelta, dtmu, dtp, tof_mu, tof_pi);
    // special cases
    // old: for the case of w/o any el cuts SetParsSpecialMomenta(p, A, func,  tof_el, tof_mu, tof_pi);

    // print fit initial parameters
    cout << "The initial parameters are: " << endl;
    for (int ip = 0; ip < func -> GetNpar(); ++ ip) {
      cout << ip << " " << func->GetParName(ip) << " " << func->GetParameter(ip) << endl;      
    }

    if (!isThreeComponentFit) {
      func -> FixParameter(4, 0.);
      func -> FixParameter(6, 0.);
    }

    // MAIN FIT!
    h3->Fit(func);

    // fits for TOF after ACT2 nad ACT3 cuts:
    TF1 *func_act2cut = (TF1*)func -> Clone(TString(func -> GetName()) + "_act2cut");
    SetInitialFitPars(isThreeComponentFit, isProtonFit, h3_act2cut, p, A, func0, func_act2cut, tof_el, delta, mdelta, dtmu, dtp, tof_mu, tof_pi);
    h3_act2cut->Fit(func_act2cut);

    TF1 *func_act3cut = (TF1*)func -> Clone(TString(func -> GetName()) + "_act3cut");
    SetInitialFitPars(isThreeComponentFit, isProtonFit, h3_act3cut, p, A, func0, func_act3cut, tof_el, delta, mdelta, dtmu, dtp, tof_mu, tof_pi);
    h3_act3cut->Fit(func_act3cut);

    // plot main fit result
    TCanvas *can_main_fit = new TCanvas("can_main_fit", "", 1200, 800);
    can_main_fit->cd()->SetTickx(kTRUE);
    can_main_fit->cd()->SetTicky(kTRUE);
    if (logy) {
      can_main_fit->cd()->SetLogy(kTRUE);
      h3->SetMaximum(20.*h3->GetMaximum());
    }
    h3->SetMarkerStyle(20);
    h3->SetMaximum(1.5*h3->GetMaximum());
    h3->Draw("ep");
    func->Draw("csame");

    // draw add also the total fits after ACT cuts:
    // can be removed
    /*
      func_act2cut->SetLineStyle(2);
      func_act2cut->Draw("csame");
      func_act3cut->SetLineStyle(3);
      func_act3cut->Draw("csame");
    */
    
    // now get the the individual components of the main fit
    TF1* func1 = 0;
    TF1* func2 = 0;
    TF1* func3 = 0;
    GetIndividualFitComponents(func, func1, func2, func3, isThreeComponentFit, "");
    int Ne, Nmu, Npi, NeErr, NmuErr, NpiErr;
    AssignFitResults(func, h1, Ne, Nmu, Npi, NeErr, NmuErr, NpiErr, isThreeComponentFit);

    // and the individual components of the fit after ACT2 amplitude cuts
    TF1* func1_act2cut = 0;
    TF1* func2_act2cut = 0;
    TF1* func3_act2cut = 0;
    GetIndividualFitComponents(func_act2cut, func1_act2cut, func2_act2cut, func3_act2cut, isThreeComponentFit, "_act2cut");
    int Ne_act2cut, Nmu_act2cut, Npi_act2cut, NeErr_act2cut, NmuErr_act2cut, NpiErr_act2cut;
    AssignFitResults(func_act2cut, h1, Ne_act2cut, Nmu_act2cut, Npi_act2cut, NeErr_act2cut, NmuErr_act2cut, NpiErr_act2cut, isThreeComponentFit);

    // could get the individual components of the fit after ACT3 cuts here as well
    // and the individual components of the fit after ACT2 amplitude cuts
    TF1* func1_act3cut = 0;
    TF1* func2_act3cut = 0;
    TF1* func3_act3cut = 0;
    GetIndividualFitComponents(func_act3cut, func1_act3cut, func2_act3cut, func3_act3cut, isThreeComponentFit, "_act3cut");
    int Ne_act3cut, Nmu_act3cut, Npi_act3cut, NeErr_act3cut, NmuErr_act3cut, NpiErr_act3cut;
    AssignFitResults(func_act3cut, h1, Ne_act3cut, Nmu_act3cut, Npi_act3cut, NeErr_act3cut, NmuErr_act3cut, NpiErr_act3cut, isThreeComponentFit);

    // finish drawing of the main fit
    can_main_fit->cd();
    func1->Draw("Csame");
    func2->Draw("Csame");          
    if (isThreeComponentFit)
      func3->Draw("Csame");
    
    tof_el = func -> GetParameter(2);
    DrawExpectedTimeArrivalLinesLegAndChi2(h1, h3, func, func1, func2, func3,
					   p, imax, 
					   isThreeComponentFit, isProtonFit,
					   Ne, Nmu, Npi,
					   NeErr, NmuErr, NpiErr,
					   tof_el, dtmu, dtpi, dtp, dtd, dtt, dta);


    // ---------------------------------------------------------------------------------------------------
    // and now draw also the histogram, main fit and indivdual components also for the ACT2cut results
    h3_act2cut->SetStats(kFALSE);
    TCanvas *can_main_fit_act2cut = new TCanvas("can_main_fit_act2cut", "", 1200, 800);
    can_main_fit_act2cut->cd()->SetTickx(kTRUE);
    can_main_fit_act2cut->cd()->SetTicky(kTRUE);
    if (logy) {
      can_main_fit_act2cut->cd()->SetLogy(kTRUE);
    }
    h3_act2cut->SetMaximum(h3->GetMaximum()); // yes, set max same as in the case of the fit w/o cuts;)
    h3_act2cut->SetMarkerStyle(20);
    h3_act2cut->Draw("ep");
    func_act2cut->Draw("csame");

    func1_act2cut->Draw("Csame");
    func2_act2cut->Draw("Csame");          
    if (isThreeComponentFit)
      func3_act2cut->Draw("Csame");
    
    tof_el = func_act2cut -> GetParameter(2);
    DrawExpectedTimeArrivalLinesLegAndChi2(h1, h3_act2cut, func_act2cut, func1_act2cut, func2_act2cut, func3_act2cut,
					   p, imax, 
					   isThreeComponentFit, isProtonFit,
					   Ne_act2cut, Nmu_act2cut, Npi_act2cut,
					   NeErr_act2cut, NmuErr_act2cut, NpiErr_act2cut,
					   tof_el, dtmu, dtpi, dtp, dtd, dtt, dta);


    // ---------------------------------------------------------------------------------------------------
    // sorry, the same code just act2cut --> act3cut
    // and now draw also the histogram, main fit and indivdual components also for the ACT3cut results
    h3_act3cut->SetStats(kFALSE);
    TCanvas *can_main_fit_act3cut = new TCanvas("can_main_fit_act3cut", "", 1200, 800);
    can_main_fit_act3cut->cd()->SetTickx(kTRUE);
    can_main_fit_act3cut->cd()->SetTicky(kTRUE);
    if (logy) {
      can_main_fit_act3cut->cd()->SetLogy(kTRUE);
    }
    h3_act3cut->SetMaximum(h3->GetMaximum()); // yes, set max same as in the case of the fit w/o cuts;)
    h3_act3cut->SetMarkerStyle(20);
    h3_act3cut->Draw("ep");
    func_act3cut->Draw("csame");

    func1_act3cut->Draw("Csame");
    func2_act3cut->Draw("Csame");          
    if (isThreeComponentFit)
      func3_act3cut->Draw("Csame");
    
    tof_el = func_act3cut -> GetParameter(2);
    DrawExpectedTimeArrivalLinesLegAndChi2(h1, h3_act3cut, func_act3cut, func1_act3cut, func2_act3cut, func3_act3cut,
					   p, imax, 
					   isThreeComponentFit, isProtonFit,
					   Ne_act3cut, Nmu_act3cut, Npi_act3cut,
					   NeErr_act3cut, NmuErr_act3cut, NpiErr_act3cut,
					   tof_el, dtmu, dtpi, dtp, dtd, dtt, dta);



    
    // some printouts    
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

    // NOW SOME NUMBERS POSTPROCESSING ;-)
    if (isProtonFit) {
      (*asciifile) << "N_proton " <<  func->GetParameter(1)/h1->GetBinWidth(1) << " " << NmuErr << endl;
    } else {
      // not a proton fit
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

	// a by hand poion shift by Jiri!!!
	//    double pioff = 0.35;
	double pioff = 0.;

	// try the fit components results after the act cut and w/o the cut
	double eff_mu = 0.;
	if (Nmu > 0.)
	  eff_mu = Nmu_act3cut /  (1.*Nmu);
	double eff_pi = 0.;
	if (Npi > 0.)
	  eff_pi = Npi_act3cut / (1.*Npi);

	// Draw the ACT2 amplitude vs TOF arrival time, before and after cuts:
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

	h2can -> Print(canname + ".png");
	h2can -> Print(canname + ".pdf");
	  
	// so far zero error on the efficiency...
	(*asciifile) << "eff_mu " << eff_mu  << " " << 0. << endl;
	(*asciifile) << "eff_pi " << eff_pi  << " " << 0. << endl;
	
      } else     {
        // two-component fit
	(*asciifile) << "N_mupi " <<  func->GetParameter(1)/h1->GetBinWidth(1) << " " << NmuErr << endl;      
    
	double eff_mupi = Nmu_act2cut /  (1.*Nmu);
	// so far zero error on the efficiency...
	(*asciifile) << "eff_mupi " << eff_mupi  << " " << 0. << endl;
	
      } // two-component fit
    } // not a proton fit

    if (asciifile) asciifile->close();

    // print:
    string name = fileName.substr(0, fileName.size()-5) + "_TOFfit.png";
    string namepdf = fileName.substr(0, fileName.size()-5) + "_TOFfit.pdf";
    can_main_fit->Print(name.c_str());
    can_main_fit->Print(namepdf.c_str());    

    /*
    string name_act2cut = fileName.substr(0, fileName.size()-5) + "_TOFfit_act2cut.png";
    string namepdf_act2cut = fileName.substr(0, fileName.size()-5) + "_TOFfit_act2cut.pdf";
    can_main_fit_act2cut->Print(name_act2cut.c_str());
    can_main_fit_act2cut->Print(namepdf_act2cut.c_str());    
    */
    
    string name_act3cut = fileName.substr(0, fileName.size()-5) + "_TOFfit_act3cut.png";
    string namepdf_act3cut = fileName.substr(0, fileName.size()-5) + "_TOFfit_act3cut.pdf";
    can_main_fit_act3cut->Print(name_act3cut.c_str());
    can_main_fit_act3cut->Print(namepdf_act3cut.c_str());    

    cout << "DONE!" << endl;
    cout << "See also " << asciiname.Data() << endl;
}


// ______________________________________________________
// ______________________________________________________
// ______________________________________________________
