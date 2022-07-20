const int nChannels = 8;
struct MidasEvent{
    int midasEvent;
    int timestamp;
    int freqsetting;
    uint64_t triggerTime;
    TH1D* waveforms[nChannels];
    
    void Initialize(){
        midasEvent = 0;
        timestamp = 0;
        freqsetting = 0;
        triggerTime = 0;
        for(int i = 0; i < nChannels; i++) waveforms[i] = NULL;
    }
};

    vector<TLine*> pedestalLines;
    vector<TLine*> timeLines;
void DrawEvent(string rawFileName, string procFileName, int nEvent){
    
    //TPad::SetTickx(kTRUE);
    MidasEvent event1;
    event1.Initialize();
    
    MidasEvent event2;
    event2.Initialize();
    
    TFile rawFile(rawFileName.c_str(), "READ");
    
    
    TTree *tree1 = (TTree*) rawFile.Get("midas_data1");

    tree1->SetBranchAddress("midasEvent",&event1.midasEvent);
    tree1->SetBranchAddress("timestamp",&event1.timestamp);
    tree1->SetBranchAddress("freqsetting",&event1.freqsetting);
    tree1->SetBranchAddress("triggerTime",&event1.triggerTime);
    for(int i = 0; i < nChannels; i++){
        tree1->SetBranchAddress(Form("Channel%d",i),&(event1.waveforms[i]));
    }

    TTree *tree2 = (TTree*) rawFile.Get("midas_data2");

    tree2->SetBranchAddress("midasEvent",&event2.midasEvent);
    tree2->SetBranchAddress("timestamp",&event2.timestamp);
    tree2->SetBranchAddress("freqsetting",&event2.freqsetting);
    tree2->SetBranchAddress("triggerTime",&event2.triggerTime);
    for(int i = 0; i < nChannels; i++){
        tree2->SetBranchAddress(Form("Channel%d",i),&(event2.waveforms[i]));
    }



    
    vector<vector<double>> *peakVoltage = NULL;
    vector<vector<double>> *peakTime = NULL;
    vector<vector<double>> *signalTime = NULL;
    vector<vector<double>> *intCharge = NULL;
    double pedestal[16];
    double pedestalSigma[16];

    
    TFile inFile(procFileName.c_str(), "READ");
    TTree *tree3 = (TTree*) inFile.Get("anaTree");

    tree3->SetBranchAddress("PeakVoltage",&peakVoltage);
    tree3->SetBranchAddress("PeakTime",&peakTime);
    tree3->SetBranchAddress("SignalTime",&signalTime);
    tree3->SetBranchAddress("IntCharge",&intCharge);
    tree3->SetBranchAddress("Pedestal",&pedestal);
    tree3->SetBranchAddress("PedestalSigma",&pedestalSigma);
    
    tree1->GetEntry(nEvent);
    tree2->GetEntry(nEvent);
    tree3->GetEntry(nEvent);
    
    TCanvas* c = new TCanvas("c", "", 1600, 1000);
    c->Divide(4, 4);
    

    for(int i = 0; i < nChannels; i++){
        c->cd(i+1)->SetMargin(0.1, 0.1, 0.1, 0.1);
        c->cd(i+1)->SetTickx(kTRUE);
        c->cd(i+1)->SetTicky(kTRUE);
        string title = "Channel " + to_string(i);// + "; t = " + Form("%.2f", signalTime[i]) + " ns;" + "A = " + Form("%.2f", peakVoltage[i]) + " V";
        event1.waveforms[i]->SetTitle(title.c_str());
        event1.waveforms[i]->SetLineColor(kBlack);
        event1.waveforms[i]->SetStats(kFALSE);
        event1.waveforms[i]->Draw("hist");
        double low = event1.waveforms[i]->GetBinLowEdge(1);
        double high = event1.waveforms[i]->GetBinLowEdge(event1.waveforms[i]->GetNbinsX()) + event1.waveforms[i]->GetBinWidth(1);
        TLine* linePed = new TLine(low, pedestal[i]*4096/2.5, high, pedestal[i]*4096/2.5);
        
        string name1 = "linePed" + to_string(i);
        string name2 = "lineTime" + to_string(i);
        linePed->SetLineStyle(2);
        linePed->SetLineWidth(2);
        linePed->SetLineColor(kRed+1);
        //linePed->SetName(name1);
        
        pedestalLines.push_back(linePed);
        pedestalLines.back()->Draw();
        
        for(int j = 0; j< signalTime->at(i).size(); j++){
             TLine* lineTime = new TLine(signalTime->at(i).at(j), event1.waveforms[i]->GetMinimum(), signalTime->at(i).at(j), event1.waveforms[i]->GetMaximum());
            //lineTime->SetName(name2.c_str());
            lineTime->SetLineStyle(2);
            lineTime->SetLineWidth(2);
            lineTime->SetLineColor(kBlue+1);
            
            timeLines.push_back(lineTime);
            timeLines.back()->Draw();       
        }

        
        
    }
    cout << "radi" << endl;
    
    //c->Update();
   for(int i = 0; i < nChannels; i++){
        c->cd(i+9)->SetMargin(0.1, 0.1, 0.1, 0.1);
        c->cd(i+9)->SetTickx(kTRUE);
        c->cd(i+9)->SetTicky(kTRUE);
        
         
        string title = "Channel " + to_string(8+i);// + "; t = " + Form("%.2f", signalTime[8+i]) + " ns;" + "A = " + Form("%.2f", peakVoltage[8+i]) + " V";
        event2.waveforms[i]->SetTitle(title.c_str());        
        event2.waveforms[i]->SetLineColor(kBlack);
        event2.waveforms[i]->SetStats(kFALSE);
        event2.waveforms[i]->Draw("hist");
        double low = event2.waveforms[i]->GetBinLowEdge(1);
        double high = event2.waveforms[i]->GetBinLowEdge(event2.waveforms[i]->GetNbinsX()) + event2.waveforms[i]->GetBinWidth(1);
        TLine* linePed = new TLine(low, pedestal[i+8]*4096/2.5, high, pedestal[i+8]*4096/2.5);
        string name1 = "linePed" + to_string(i);
        string name2 = "lineTime" + to_string(i);
        linePed->SetLineStyle(2);
        linePed->SetLineWidth(2);
        linePed->SetLineColor(kRed+1);
       // linePed->SetName(name1);
        
        pedestalLines.push_back(linePed);
        pedestalLines.back()->Draw();
        
       for(int j = 0; j< signalTime->at(i+8).size(); j++){
             TLine* lineTime = new TLine(signalTime->at(i+8).at(j), pedestal[i+8]*4096/2.5, signalTime->at(i+8).at(j), pedestal[i+8]*4096/2.5-200);
            //lineTime->SetName(name2.c_str());
            lineTime->SetLineStyle(2);
            lineTime->SetLineWidth(3);
            lineTime->SetLineColor(kBlue);
            
            timeLines.push_back(lineTime);
            timeLines.back()->Draw();       
        }
        
    }
    
    c->Draw();
    c->SaveAs("event.png");
    
}
