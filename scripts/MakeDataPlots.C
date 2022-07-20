double GetBeta(double mass, double mom){
    double bg = mom/mass;
    double beta = sqrt(bg*bg/(1+bg*bg));
    return beta;
}


void MakeDataPlots(string fileName, int mom){
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
    
    TH1D hTOF("hRef_TOFAll", "", 250, -100, 150);
    TH1D hTOFAll("hTOFAll", "", 120, 37.5, 43.5);
    TH1D hTOFAllWide("hTOFAllWide", "", 720, 37.5, 73.5);
    TH1D hTOFEl("hTOFEl", "", 120, 37.5, 43.5);
    TH1D hTOFOther("hTOFOther", "", 120, 37.5, 43.5);
    
    vector<TH1D> hCharge;
    vector<TH1D> hVoltage;
    vector<TH1D> hHit;
    vector<TH1D> hPedestalSigma;
    
    TH2D hTOFACT1V("hRef_TOFACT1V", "; t_{1}-t_{0} [ns]; ACT1 Amplitude", 200, 37.5, 43.5, 200, 0., 1.6);
    TH2D hTOFACT2V("hRef_TOFACT2V", "; t_{1}-t_{0} [ns]; ACT2 Amplitude", 200, 37.5, 43.5, 200, 0., 0.1);
    TH2D hTOFACT3V("hRef_TOFACT3V", "; t_{1}-t_{0} [ns]; ACT3 Amplitude", 200, 37.5, 42.5, 200, 0., 0.1);
    
    TH2D hTOFACT1C("hRef_TOFACT1C", "; t_{1}-t_{0} [ns]; ACT1 Charge", 200, 37.5, 43.5, 200, 0., 0.016);
    TH2D hTOFACT2C("hRef_TOFACT2C", "; t_{1}-t_{0} [ns]; ACT2 Charge", 200, 37.5, 43.5, 200, 0., 0.002);
    TH2D hTOFACT3C("hRef_TOFACT3C", "; t_{1}-t_{0} [ns]; ACT3 Charge", 200, 37.5, 43.5, 200, 0., 0.002);
    
    
    
    for(int i = 0; i < 16; i++){
        string name1 = "hRef_Charge" + to_string(i);
        string name2 = "hRef_Voltage" + to_string(i);
        string name3 = "hRef_Hits" + to_string(i);
        string name4 = "hRef_PedestalSigma" + to_string(i);
        
        string title1 = "Channel " + to_string(i) + "; Charge [nC]; Triggers";
        string title2 = "Channel " + to_string(i) + "; Total Amplitude [V]; Triggers";
        string title3 = "Channel " + to_string(i) + "; Hits per trigger; Triggers";
        string title4 = "Channel " + to_string(i) + "; #sigma_{ped} [V]; Triggers";
        TH1D temp1(name1.c_str(), title1.c_str(), 200, 0., 0.08);
        TH1D temp2(name2.c_str(), title2.c_str(), 200, 0., 0.8);
        TH1D temp3(name3.c_str(), title3.c_str(), 5, -0.5, 4.5);
        TH1D temp4(name4.c_str(), title4.c_str(), 200, 0., 0.01);
        
        hCharge.push_back(temp1);
        hVoltage.push_back(temp2);
        hHit.push_back(temp3);
        hPedestalSigma.push_back(temp4);
    }
    
    for(int i = 0; i < ent; i++){
        tree->GetEntry(i);
        

        vector<int> indices;
        

        for(int j = 0; j < 16; j++){
        
            int ind = std::max_element(peakVoltage->at(j).begin(),peakVoltage->at(j).end()) - peakVoltage->at(j).begin();
            indices.push_back(ind);
            hCharge.at(j).Fill(intCharge->at(j).at(ind));
            hVoltage.at(j).Fill(peakVoltage->at(j).at(ind));
            hHit.at(j).Fill(peakVoltage->at(j).size());
            hPedestalSigma.at(j).Fill(pedestalSigma[j]);

        }
  
  
        double t0 = (signalTime->at(8).at(indices.at(8)) + signalTime->at(9).at(indices.at(9)) + signalTime->at(10).at(indices.at(10)) + signalTime->at(11).at(indices.at(11)))/4.;
        double t1 = (signalTime->at(12).at(indices.at(12)) + signalTime->at(13).at(indices.at(13)) + signalTime->at(14).at(indices.at(14)) + signalTime->at(15).at(indices.at(15)))/4.;

        hTOFACT1V.Fill(t1-t0, peakVoltage->at(0).at(indices.at(0)) + peakVoltage->at(1).at(indices.at(1)));
        hTOFACT2V.Fill(t1-t0, peakVoltage->at(2).at(indices.at(2)) + peakVoltage->at(3).at(indices.at(3)));
        hTOFACT3V.Fill(t1-t0, peakVoltage->at(4).at(indices.at(4)) + peakVoltage->at(5).at(indices.at(5)));
        
        hTOFACT1C.Fill(t1-t0, intCharge->at(0).at(indices.at(0)) + intCharge->at(1).at(indices.at(1)));
        hTOFACT2C.Fill(t1-t0, intCharge->at(2).at(indices.at(2)) + intCharge->at(3).at(indices.at(3)));
        hTOFACT3C.Fill(t1-t0, intCharge->at(4).at(indices.at(4)) + intCharge->at(5).at(indices.at(5)));
        
        hTOF.Fill(t1-t0);
                      
        bool pass = true;
        bool isEl = false;
        switch(mom){
            case 200: {
                if(peakVoltage->at(0).at(indices.at(0))+peakVoltage->at(1).at(indices.at(1)) > 0.10)
                    isEl = true;
                if(peakVoltage->at(2).at(indices.at(2))+peakVoltage->at(3).at(indices.at(3))>0.010)
                    isEl = true; 
                if(peakVoltage->at(4).at(indices.at(4))+peakVoltage->at(5).at(indices.at(5))>0.012)
                    isEl = true;
                                                
                break;
            }
            
            case -200: {
                if(peakVoltage->at(0).at(indices.at(0))+peakVoltage->at(1).at(indices.at(1)) > 0.10)
                    isEl = true;
                if(peakVoltage->at(2).at(indices.at(2))+peakVoltage->at(3).at(indices.at(3))>0.010)
                    isEl = true; 
                if(peakVoltage->at(4).at(indices.at(4))+peakVoltage->at(5).at(indices.at(5))>0.012)
                    isEl = true;
                                                
                break;
            }
            case 220: {

                /*double pedestalSigmaCut[16] = {0.0014, 0.0018, 0.0014, 0.0012, 0.0012, 0.0012, 0.0028, 0.0014, 0.0018, 0.0014, 0.0014, 0.0014, 0.0014, 0.0014, 0.0016, 0.0014};
                for(int j = 0; j < 16; j++){
                
                    if(pedestalSigma[j] > pedestalSigmaCut[j]){
                        pass = false;
                        break;
                    }
                }*/

                if(peakVoltage->at(0).at(indices.at(0))+peakVoltage->at(1).at(indices.at(1)) > 0.12)
                    isEl = true;
                if(peakVoltage->at(2).at(indices.at(2))+peakVoltage->at(3).at(indices.at(3))>0.018)
                    isEl = true; 
                if(peakVoltage->at(4).at(indices.at(4))+peakVoltage->at(5).at(indices.at(5))>0.035)
                    isEl = true;
                                    
                break;
            }
            
            case -220: {
                /*double pedestalSigmaCut[16] = {0.0014, 0.0018, 0.0014, 0.0012, 0.0012, 0.0012, 0.0028, 0.0014, 0.0018, 0.0014, 0.0014, 0.0014, 0.005, 0.009, 0.006, 0.006};
                for(int j = 0; j < 16; j++){
                
                    if(pedestalSigma[j] > pedestalSigmaCut[j]){
                        pass = false;
                        break;
                    }
                }*/
                   
                if(peakVoltage->at(0).at(indices.at(0))+peakVoltage->at(1).at(indices.at(1)) > 0.18)
                    isEl = true;
                if(peakVoltage->at(2).at(indices.at(2))+peakVoltage->at(3).at(indices.at(3))>0.025)
                    isEl = true; 
                if(peakVoltage->at(4).at(indices.at(4))+peakVoltage->at(5).at(indices.at(5))>0.022)
                    isEl = true;
                
                
                break;
            }
            
            
	        case -300: { //; jk
                /*double pedestalSigmaCut[16] = {0.0014, 0.0018, 0.0014, 0.0012, 0.0012, 0.0012, 0.0028, 0.0014, 0.0018, 0.0014, 0.0014, 0.0014, 0.005, 0.009, 0.006, 0.006};
                for(int j = 0; j < 16; j++){
                
                    if(pedestalSigma[j] > pedestalSigmaCut[j]){
                        pass = false;
                        break;
                    }
                }*/
                   
                if(peakVoltage->at(0).at(indices.at(0))+peakVoltage->at(1).at(indices.at(1)) > 0.18)
                    isEl = true;
                if(peakVoltage->at(2).at(indices.at(2))+peakVoltage->at(3).at(indices.at(3))>0.025)
                    isEl = true; 
                if(peakVoltage->at(4).at(indices.at(4))+peakVoltage->at(5).at(indices.at(5))>0.022)
                    isEl = true;
                
                
                break;
            }




	        case 400: {
	        
                double voltageCut[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025, 0.025, 0.025, 0.025, 0.025, 0.040, 0.025, 0.025};
                for(int j = 0; j < 16; j++){
                
                    if(peakVoltage->at(j).at(indices.at(j)) < voltageCut[j]){
                        pass = false;
                        break;
                    }
                }
	            break;
	            
	        }

            default: {
                cout << "Settings not implemented for " << mom << " MeV/c beam" << endl;
                break;
            }
        }
        





        

        if(!pass) continue;
        hTOFAll.Fill(t1-t0);
        hTOFAllWide.Fill(t1-t0);

        if(isEl){
            hTOFEl.Fill(t1-t0);
        }
        else{
            hTOFOther.Fill(t1-t0);
        }

    }
    
    string outFileName = fileName.substr(0, fileName.size()-5) + "_plots.root";
    
    TFile outFile(outFileName.c_str(), "RECREATE");
    outFile.cd();
    

    hTOFACT1V.Write();
    hTOFACT2V.Write();
    hTOFACT3V.Write();
    hTOFACT1C.Write();
    hTOFACT2C.Write();
    hTOFACT3C.Write();
    hTOF.Write();
    for (auto hist: hVoltage) hist.Write();
    for (auto hist: hCharge) hist.Write();
    for (auto hist: hHit) hist.Write();
    for (auto hist: hPedestalSigma) hist.Write();


    hTOFAll.Write();
    hTOFAllWide.Write();
    hTOFEl.Write();
    hTOFOther.Write();

}


