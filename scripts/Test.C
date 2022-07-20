double GetBeta(double mass, double mom){
    double bg = mom/mass;
    double beta = sqrt(bg*bg/(1+bg*bg));
    return beta;
}


void Test(string fileName, int mom){
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
    
    vector<TH1D> hTime;
    
     TH1D hTOF("hRef_TOFAll", "", 250, -100, 150);
    int ent = tree->GetEntries();
    //double timeCut[8][2] = {{185, 215}, {125, 160}, {125, 160}, {135, 170}, {170, 205}, {185, 215}, {185, 215}, {185, 215}};
    double timeCut[8][2] = {{185, 300}, {125, 300}, {125, 300}, {135, 300}, {170, 300}, {185, 300}, {185, 300}, {185, 300}};
    for(int i = 0; i < ent; i++){
        tree->GetEntry(i);   
        
        vector<int> indices(16, -1);
        
        
        int nFound = 0;
        for(int j = 8; j <16; j++){
        
            bool found = false;
            int ind = 0;
            double max = -1;
            for(int k = 0; k < signalTime->at(j).size(); k++){
                if(signalTime->at(j).at(k) < timeCut[j-8][0] || signalTime->at(j).at(k) > timeCut[j-8][1]) continue;
                
                if(peakVoltage->at(j).at(k) < max) continue;
                
                    
                ind = k;
                found = true;           

               
            }
            
            if(found){
                nFound = nFound | (1 << (j-8));
                indices.at(j) = ind;
                //nFound ++;         
            }
        }
        std::bitset<8> c(nFound & (1));
        cout << "Event: " << i << ", N found = " << c << endl;
        if(nFound != 255) continue;

        double t0 = (signalTime->at(8).at(indices.at(8)) + signalTime->at(9).at(indices.at(9)) + signalTime->at(10).at(indices.at(10)) + signalTime->at(11).at(indices.at(11)))/4.;
        double t1 = (signalTime->at(12).at(indices.at(12)) + signalTime->at(13).at(indices.at(13)) + signalTime->at(14).at(indices.at(14)) + signalTime->at(15).at(indices.at(15)))/4.;
        hTOF.Fill(t1-t0);
           
    }
    
    TFile outFile("test.root", "RECREATE");
    outFile.cd();
    hTOF.Write(); 
}
