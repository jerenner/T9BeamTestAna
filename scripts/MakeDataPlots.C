double GetBeta(double mass, double mom){
    double bg = mom/mass;
    double beta = sqrt(bg*bg/(1+bg*bg));
    return beta;
}


void MakeDataPlots(string fileName){
    vector<vector<double>> *peakVoltage = NULL;
    vector<vector<double>> *peakTime = NULL;
    vector<vector<double>> *signalTime = NULL;
    vector<vector<double>> *intCharge = NULL;
    double pedestal[16];
    double pedestalSigma[16];

    
    TFile inFile(fileName.c_str(), "READ");
    TTree *tree = (TTree*) inFile.Get("anaTree");

    tree->SetBranchAddress("PeakVoltage",&peakVoltage);
    tree->SetBranchAddress("PeakTime",&peakTime);
    tree->SetBranchAddress("SignalTime",&signalTime);
    tree->SetBranchAddress("IntCharge",&intCharge);
    tree->SetBranchAddress("Pedestal",&pedestal);
    tree->SetBranchAddress("PedestalSigma",&pedestalSigma);
    //tree->SetBranchAddress("PassThreshold",&passThreshold);
    
    int ent = tree->GetEntries();
    
    
    TH1D hTOFAll("hTOFAll", "", 100, 37.5, 42.5);
    TH1D hTOFEl("hTOFEl", "", 100, 37.5, 42.5);
    TH1D hTOFOther("hTOFOther", "", 100, 37.5, 42.5);
    
    vector<TH1D> hCharge;
    vector<TH1D> hVoltage;
    vector<TH1D> hHit;
    
    TH2D hTOFACT1V("hTOFACT1V", "", 200, 37.5, 42.5, 200, 0., 1.6);
    TH2D hTOFACT2V("hTOFACT2V", "", 200, 37.5, 42.5, 200, 0., 0.1);
    TH2D hTOFACT3V("hTOFACT3V", "", 200, 37.5, 42.5, 200, 0., 0.1);
    
    TH2D hTOFACT1C("hTOFACT1C", "", 200, 37.5, 42.5, 200, 0., 1.6);
    TH2D hTOFACT2C("hTOFACT2C", "", 200, 37.5, 42.5, 200, 0., 0.1);
    TH2D hTOFACT3C("hTOFACT3C", "", 200, 37.5, 42.5, 200, 0., 0.1);
    
    
    
    for(int i = 0; i < 16; i++){
        string name1 = "hCharge" + to_string(i);
        string name2 = "hVoltage" + to_string(i);
        string name3 = "hHits" + to_string(i);
        TH1D temp1(name1.c_str(), "", 200, 0., 0.8);
        TH1D temp2(name2.c_str(), "", 200, 0., 0.8);
        TH1D temp3(name3.c_str(), "", 5, -0.5, 4.5);
        
        hCharge.push_back(temp1);
        hVoltage.push_back(temp2);
        hHit.push_back(temp3);
    }
    
    for(int i = 0; i < ent; i++){
        tree->GetEntry(i);
        
        bool pedestalOK = true;
        for(int j = 8; j < 16; j++){
            if(pedestalSigma[j] > 0.007){
                pedestalOK = false;
                break;
            }
        }
        
        if(!pedestalOK) continue;

        
        vector<int> indices;
        

        for(int j = 0; j < 16; j++){
        
            int ind = std::max_element(peakVoltage->at(j).begin(),peakVoltage->at(j).end()) - peakVoltage->at(j).begin();
            indices.push_back(ind);
            hCharge.at(j).Fill(intCharge->at(j).at(ind));
            hVoltage.at(j).Fill(peakVoltage->at(j).at(ind));
            hHit.at(j).Fill(peakVoltage->at(j).size());
            
            

        }


        double t0 = (signalTime->at(8).at(indices.at(8)) + signalTime->at(9).at(indices.at(9)) + signalTime->at(10).at(indices.at(10)) + signalTime->at(11).at(indices.at(11)))/4.;
        double t1 = (signalTime->at(12).at(indices.at(12)) + signalTime->at(13).at(indices.at(13)) + signalTime->at(14).at(indices.at(14)) + signalTime->at(15).at(indices.at(15)))/4.;

        hTOFACT1V.Fill(t1-t0, peakVoltage->at(0).at(indices.at(0)) + peakVoltage->at(1).at(indices.at(1)));
        hTOFACT2V.Fill(t1-t0, peakVoltage->at(2).at(indices.at(2)) + peakVoltage->at(3).at(indices.at(3)));
        hTOFACT3V.Fill(t1-t0, peakVoltage->at(4).at(indices.at(4)) + peakVoltage->at(5).at(indices.at(5)));
        
        hTOFACT1C.Fill(t1-t0, intCharge->at(0).at(indices.at(0)) + intCharge->at(1).at(indices.at(1)));
        hTOFACT2C.Fill(t1-t0, intCharge->at(2).at(indices.at(2)) + intCharge->at(3).at(indices.at(3)));
        hTOFACT3C.Fill(t1-t0, intCharge->at(4).at(indices.at(4)) + intCharge->at(5).at(indices.at(5)));
        
        
       /* if(intCharge->at(8).at(indices.at(8)) < 0.06) continue;
        if(intCharge->at(9).at(indices.at(9)) < 0.05) continue;
        if(intCharge->at(10).at(indices.at(10)) < 0.07) continue;
        if(intCharge->at(11).at(indices.at(11)) < 0.06) continue;
        if(intCharge->at(12).at(indices.at(12)) < 0.05) continue;
        if(intCharge->at(13).at(indices.at(13)) < 0.06) continue;
        if(intCharge->at(14).at(indices.at(14)) < 0.035) continue;
        if(intCharge->at(15).at(indices.at(15)) < 0.045) continue;*/

        hTOFAll.Fill(t1-t0);
        /*if(peakVoltage->at(0).at(indices.at(0))+peakVoltage->at(1).at(indices.at(1)) > 0.1 || 
            peakVoltage->at(2).at(indices.at(2))+peakVoltage->at(3).at(indices.at(3))>0.012 || 
            peakVoltage->at(4).at(indices.at(4))+peakVoltage->at(5).at(indices.at(5))>0.04){*/

        if(peakVoltage->at(0).at(indices.at(0))+peakVoltage->at(1).at(indices.at(1)) > 0.1 || 
            peakVoltage->at(2).at(indices.at(2))+peakVoltage->at(3).at(indices.at(3))>0.025 || 
            peakVoltage->at(4).at(indices.at(4))+peakVoltage->at(5).at(indices.at(5))>0.05){
            hTOFEl.Fill(t1-t0);
        }
        else{
            hTOFOther.Fill(t1-t0);
        }

    }
    
    string outFileName = fileName.substr(0, fileName.size()-5) + "_plots.root";
    
    TFile outFile(outFileName.c_str(), "RECREATE");
    outFile.cd();
    
    hTOFAll.Write();
    hTOFEl.Write();
    hTOFOther.Write();
    hTOFACT1V.Write();
    hTOFACT2V.Write();
    hTOFACT3V.Write();
    hTOFACT1C.Write();
    hTOFACT2C.Write();
    hTOFACT3C.Write();
    
    for(int j = 0; j < 16; j++){
        hVoltage.at(j).Write();
        hCharge.at(j).Write();
        hHit.at(j).Write();
    }

}


