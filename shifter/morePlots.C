// jk
// 22.7.2023
// then developed by others


// ----------------------------------------------------------
/*
void makeLine(x1, x2, y1, y2)
{
  TLine line = new TLine(x1, y1, x2, y2);
  line -> SetLineColor(ROOT.kGreen);
  line -> SetLineWidth(2);
  line -> Draw();
  return line;
}
*/

// ----------------------------------------------------------

void morePlots(TString filename = "") {
  gStyle->SetPaintTextFormat("1.2f");
  TString pngdir = "png_results/";
  TString pdfdir = "pdf_results/";
  gSystem->Exec("mkdir -p " + pngdir);
  gSystem->Exec("mkdir -p " + pdfdir);

  TString opt2d = "colz";
    
  gStyle -> SetOptFit(111);
  gStyle -> SetPalette(kSolar);

  vector<TCanvas*> cans;
  
  TFile *rfile = new TFile(filename, "read");
  TString hnames[] = {
			"hTOFAll",
			"hTOFAllLow",
			"hTOFAllWide",
			"hRef_pbA_act23A",
			"hRef_TOFACT23A",
			"hRef_TOFPbA",
			"hnHitsHodoscope",
			"LeadGlassPhotonAVsPositronHodoOcc",
			"LeadGlassPhotonAVsPositronMaxHodoOcc",
			"HodoOccScatter",
			"HodoOccScatterFrac"
			

  };
  int nhs = sizeof(hnames) / sizeof(TString);
  
  TString ftag = filename.ReplaceAll("output_","").ReplaceAll("_plots.root","").ReplaceAll("/","_more");
  gSystem->Exec("mkdir -p pdf png");
  for (int ih = 0; ih < nhs; ih++ ) {
    TString hname = hnames[ih];
    TH1D *h = (TH1D*) rfile -> Get(hname);
    TString canname = Form("WCTEJuly2023_Quick1D_%s_%s",ftag.Data(), hname.Data());
    canname = canname.ReplaceAll("_list_root","").ReplaceAll("_ntuple","");
    TCanvas *can = new TCanvas(canname, canname, ih*350, ih*100, 800, 800);
    if (h)  {
      h -> SetStats(111111);
      
      TString opt = "colz";
      if (hname.Contains("Frac")) {
	opt = "colz text";
      }
      if (hname.Contains("hTOF")) {
	opt = "hist";	
	h -> SetFillColor(kPink+9);
	h -> SetFillStyle(1111);
	gPad -> SetLogy(1);
      } else {
	gPad -> SetLogz(1);
      }
      h -> Draw(opt);
      if (hname.Contains("Frac")) {
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	gPad->Update();
      }
    }
  }

  
  //##################################
  //#       plots all the canvas     #
  //##################################
  
  //srun = "";
  //tokens = filename.split("_");
  //for token in tokens:
  //if (token.Contains("00") {
  //	srun = token.replace("000","");
  //  }
  /*
  for (auto can : cans) {
    can -> cd();
    can -> Update();
    can -> Print(pngdir + can -> GetName() + ".png");
    can -> Print(pdfdir + can -> GetName() + ".pdf");
  }
  */
  //if (!gBatch)
  //  gApplication -> Run();
  
} // macro
  
