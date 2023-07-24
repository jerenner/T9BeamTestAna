// jk
// 22.7.2023




map<int,TString> ChNames =  { {0, "ACT-00"}, {1, "ACT-01"}, {2, "ACT-10"}, {3, "ACT-11"}, {4, "ACT-20"}, {5, "ACT-21"}, {6, "ACT-30"}, {7, "ACT-31"},
			      {8, "TOF-00"}, {9, "TOF-01"}, {10, "TOF-10"}, {11, "TOF-11"}, {12, "TOF-20"}, {13, "TOF-21"}, {14, "TOF-30"}, {15, "TOF-31"},
			      {16, "HC-00"}, {17, "HC-10"}, {18, "LG-00"},
			      {19, "X"}, {20, "X"}, {21, "X"}, {22, "X"}, {23, "X"},
			      {24, "X"}, {25, "X"}, {26, "X"}, {27, "X"}, {28, "X"}, {29, "X"}, {30, "X"}, {31, "X"} };


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

void quickPlots1d(TString filename = "") {

  TString pngdir = "png_results/";
  TString pdfdir = "pdf_results/";
  gSystem->Exec("mkdir " + pngdir);
  gSystem->Exec("mkdir " + pdfdir);

  TString opt2d = "colz";
    
  gStyle -> SetOptFit(111);
  gStyle -> SetPalette(kSolar);

  vector<TCanvas*> cans;
  
  TFile *rfile = new TFile(filename, "read");
  map<TString,int> hbasenames = {
				 {"hRef_nPeaks", kYellow},
				 {"hRef_Time", kGreen},
				 {"hRef_Charge", kCyan},
				 {"hRef_Voltage", kMagenta}

  };
  
  int nChannels = 19; // 32
  
  TString ftag = filename.ReplaceAll("output_","").ReplaceAll("_plots.root","").ReplaceAll("/","_");
  gSystem->Exec("mkdir -p pdf png");
  int ih = -1;
  for ( const auto &[hbasename, col]: hbasenames ) {
    ih++;
    vector<TH1D*> hs;

    for (int ich = 0; ich <  nChannels; ich++) {
      TString hname = hbasename;
      hname += ich;
      //cout << "Will try to get histo names " << hname.Data() << endl;
      TH1D *h = (TH1D*) rfile -> Get(hname);
      //cout << "Pushing " << ich << " " <<  hname << endl;
      hs.push_back(h);
    } // channels
    
    TString canname = Form("WCTEJuly2023_Quick1D_%s_%s",ftag.Data(), hbasename.Data());
    canname = canname.ReplaceAll("_list_root","").ReplaceAll("_ntuple","");
    int ich = -1;
    int idigi = 0;
    for (auto h : hs) {
      ich++;
      if (ich % 8 == 0) {
	idigi = ich / 8;
      
	int off = 60;
	int cw = 4*400 + 4*off;
	int ch = 2*400;
	if (idigi > 1) {
	  cw = 2*400 + 2*off;
	  ch = 400;
	}
	TCanvas *can = new TCanvas(canname + Form("_digi%i", idigi), canname + Form("_digi%i", idigi), 0, 0, cw, ch);
	if (idigi < 2) {
	  can -> Divide(4,2);
	} else {
	  can -> Divide(3,1);
	}
	cans.push_back(can);
      }
      TCanvas *can = cans[cans.size()-1];
      can -> cd(ich % 8 + 1);
      if (h) {
	h -> SetStats(0);
	// if not "Time" in h ->GetName():
	if (TString(h ->GetName()).Contains("nPeaks")) {
	  gPad -> SetLogy(1);
	  //h ->GetYaxis().SetRangeUser(1.e-4, h ->GetYaxis().GetXmax());
	} else {
	  cout << "got null pointer histo!" << endl;
	}
      }
      h -> SetFillColor(col);
      h -> SetFillStyle(1111);
      h -> SetTitle(ChNames[ich]);
      h -> Draw("hist");
    } // hs
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
  
