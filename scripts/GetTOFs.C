#include "midas_data1.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <unistd.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TF1.h>
#include <TMarker.h>

double GetBaseline(TH1D *Ch){
	double baseline = 0;
	for (int j = 0; j < 30; j++) baseline += Ch->GetBinContent(j);
	return baseline/30;
}
int GetMin(TH1D *Ch, int w_s, int w_e, int &A, double baseline){
	int i_min;
	int A_min = 99999;
	//int PtsAfterPeak = 0;
	for (int i = w_e; i > w_s; --i){
		if (Ch->GetBinContent(i) < A_min){
			i_min = i;
			A_min = Ch->GetBinContent(i);
			//PtsAfterPeak = 0;
		}
		/*else {
			++PtsAfterPeak;
			int q = 0;
			if ( PtsAfterPeak == 6) for ( int j = - 7; j < 0; ++j){
				q += baseline //proste nakej algoritmus co po peti dalsi iteracich porovna integraci 3 vs dalsi 4 .. treba... nebo 3 vs 2
			}
		}*/
	}
	A = A_min;
	return i_min;
}

int GetMin2(TH1D *Ch, int w_s, int w_e, int &A, double baseline ){ 
	
	int i_min = 0;
	int A_min = 99999;
	int i = w_e;
	
	while ( i > w_s ){
		if (Ch->GetBinContent(i) < A_min){
			i_min = i;
			A_min = Ch->GetBinContent(i);
		}
	--i;		
	}
	A = A_min;
	if (Ch->GetBinContent(i) > Ch->GetBinContent(i-1)) {
		while (Ch->GetBinContent(i) >= Ch->GetBinContent(i-1)) i++;
		return GetMin2(Ch, i, w_e, A , baseline);			
	}
	else return i_min;
}

int GetRiseEdgeTime(TH1D *Ch, int w_s, double baseline, int A, int i_min){
	double thr = (0.2 * ((double) A - baseline)) + baseline;
	int i = i_min;
	while ((Ch->GetBinContent(i+1) < thr) && (i > 1)){
		--i;
	}
	return i;
	
}

void GetTOFs(){
	char * fpath = (char*) "root_run_000143_0000_clean.root";// change later to some input file via args
	TFile * f = TFile::Open(fpath);	
	
	TTree * tr = (TTree*) f->Get("midas_data2");
	midas_data1 *m = new midas_data1(tr);
	
	int w_s[] = {250, 180, 180, 200, 290, 290, 290, 270};
	int w_e[] = {400, 400, 400, 400, 400, 400, 400, 400};
	TH1D *Chs[] = {m->Channel0, m->Channel1, m->Channel2, m->Channel3, m->Channel4, m->Channel5, m->Channel6, m->Channel7};
	Long64_t nentries = m->fChain->GetEntriesFast();
	
	const bool visual = false; 
	
	TCanvas *c = new TCanvas("c","TOF estimation",1860,1000);
	c->Divide(2,4);
	TMarker *p_edge[8];
	TMarker *p_min[8];
	for (int i = 0; i < 8; ++i){
		p_edge[i] = new TMarker(0,0,i+1);
		p_edge[i]->SetMarkerStyle(8);
		p_edge[i]->SetMarkerColor(3);
		p_edge[i]->SetMarkerSize(1);
		p_min[i] = new TMarker(0,0,i+9);
		p_min[i]->SetMarkerStyle(8);
		p_min[i]->SetMarkerColor(2);
		p_min[i]->SetMarkerSize(1);
			
		c->cd(i+1);
		Chs[i]->Draw();
		p_edge[i]->Draw();
		p_min[i]->Draw();
	}
	c->Modified();
	c->Update();
	if (not visual) c->Close();
	
	TCanvas *c_tof = new TCanvas("c_tof"," ");
	TH1D *h_tof = new TH1D("h_tof","TOF [ns]",36, 30,60);
	for ( Long64_t jentry=0; jentry<nentries; jentry++) {
		Long64_t ientry = m->LoadTree(jentry);
		if (ientry < 0) break;
		m->fChain->GetEntry(ientry);
		//TOF0
		double TOF0 = 0;
		double TOF1 = 0;
		
		//cout<<endl<<endl<<"Event #"<<ientry<<":"<<endl<<"__________________________________________________________"<<endl;
		for (int j = 0; j < 8; ++j){
			TH1D * Ch = Chs[j];
			double baseline = GetBaseline( Ch);									
			int A;
			int i_min = GetMin2(Ch, w_s[j], w_e[j], A, baseline);
			int i_edge = GetRiseEdgeTime(Ch, w_s[j], baseline, A, i_min);
			double t_edge = Ch->GetBinCenter(i_edge);
			
			
			if (j < 4 )TOF0 += t_edge;
			else TOF1 += t_edge;
			
			if (visual){
			cout<<"   Baseline = "<<baseline<<", A = "<<A<<", i_min = "<<i_min<<", TOF0 = "<<t_edge<<endl;
				c->cd(j+1);
				Ch->Draw();
				cout<<"   X = "<<t_edge<<", Y = "<<Ch->GetBinContent(i_edge)<< endl;
				p_edge[j]->SetX(t_edge);
				p_edge[j]->SetY(Ch->GetBinContent(i_edge));
				p_edge[j]->Draw();
				p_min[j]->SetX(Ch->GetBinCenter(i_min));
				p_min[j]->SetY(Ch->GetBinContent(i_min));
				p_min[j]->Draw();

			
			
				c->Modified();
				c->Update();
			}
		}
		TOF0 = TOF0/4;
		TOF1 = TOF1/4;
		h_tof->Fill(TOF1 - TOF0);
		//cout<<"   --->TOF0 = "<<TOF0<<endl;
		if (visual) {
			cout<<"Event #"<<ientry<<endl;
			cin.get();
		}
//events with bad peak identification 3-w14, 4-w8w13
		if(ientry%1000 == 0) cout<<ientry<<"/"<<nentries<<endl;
			
	}
	c_tof->cd();
	h_tof->Draw();	
		

}
