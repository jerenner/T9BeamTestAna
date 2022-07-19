#ifndef ANALYSISCONFIG_H
#define ANALYSISCONFIG_H

#include <string>
#include <vector>
#include <iostream>
//#define DEBUG 0
class AnalysisConfig {

	public:

		AnalysisConfig (std::string configFile);
		~AnalysisConfig (){}

		void Print();

        std::string GetConfigFile(){return fConfigFile;}
        int GetNumberOfChannels(){return fNChannels;}
        double GetVoltageScale(){return fVoltageScale;}

        bool IsActive(int i);
        bool MeasureCharge(int i);
        bool MeasureTime(int i);

        double GetPedestalBinLow(int i);
        double GetPedestalBinHigh(int i);

        double GetAnalysisBinLow(int i);
        double GetAnalysisBinHigh(int i);

        int GetPolarity(int i);
        double GetThreshold(int i);

        void SetConfigFile(std::string configFile);
 	private:
        void ReadConfig();
        void CheckBin(int i);

        std::string fConfigFile;
        int fNChannels;
        double fVoltageScale;

        std::vector<double> fThresholds;

        std::vector<bool> fActiveChannels;
        std::vector<bool> fChargeMeasurements;
        std::vector<bool> fTimeMeasurements;

        std::vector<double> fPedestalWindowLow;
        std::vector<double> fPedestalWindowHigh;

        std::vector<double> fAnalysisWindowLow;
        std::vector<double> fAnalysisWindowHigh;

        std::vector<int> fPolarity;
};

#endif
