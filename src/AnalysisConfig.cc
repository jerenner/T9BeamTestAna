#include "AnalysisConfig.h"
#include <iostream>
#include <stdlib.h>
#include "Utility.h"
#include <algorithm>
#include <fstream>
#include "jsmn.h"
#include <sstream>
#include "JSONError.h"

using namespace std;
using namespace utl;

AnalysisConfig::AnalysisConfig (std::string configFile){
    fConfigFile = configFile;
    CheckFile(fConfigFile);
    ReadConfig();
}

void AnalysisConfig::SetConfigFile(std::string configFile){
    fConfigFile = configFile;
    CheckFile(fConfigFile);
    ReadConfig();
}
void AnalysisConfig::ReadConfig(){

    fActiveChannels.clear();
    fChargeMeasurements.clear();
    fTimeMeasurements.clear();

    fPedestalWindowLow.clear();
    fPedestalWindowHigh.clear();
    fAnalysisWindowLow.clear();
    fAnalysisWindowHigh.clear();

	ifstream in(fConfigFile);
	string contents((istreambuf_iterator<char>(in)),
		istreambuf_iterator<char>());


	contents.erase(remove(contents.begin(), contents.end(), '\t'), contents.end());
	contents.erase(remove(contents.begin(), contents.end(), '\n'), contents.end());
	contents.erase(remove(contents.begin(), contents.end(), ' '), contents.end());
	jsmn_parser p;
	jsmntok_t t[256*2]; // modified by Jiri 2023 from 256 to 256*2


	jsmn_init(&p);
	int r = jsmn_parse(&p, contents.c_str(), contents.size(), t, sizeof(t)/sizeof(t[0]));
	if (r < 0) {
		cerr << "Error in " << __func__ << ", line " << __LINE__ << " in " << __FILE__ << endl;
		cerr << "Failed to parse JSON file. Exiting..." << endl;

		exit(EXIT_FAILURE);
	}

	if (r < 3) {
		cerr << "Error in " << __func__ << ", line " << __LINE__ << " in " << __FILE__ << endl;
		cerr << "Number of tokens too small. Exiting..." << endl;
		exit(EXIT_FAILURE);
	}
	CheckJSMNType(t[1], JSMN_STRING, __LINE__, __func__, __FILE__);
	CheckJSMNType(t[2], JSMN_OBJECT, __LINE__, __func__, __FILE__);


	string configName = contents.substr(t[1].start, t[1].end-t[1].start);
	if(configName != "WaveAnalysis"){
		cerr << "Error in " << __func__ << ", line " << __LINE__ << " in " << __FILE__ << endl;
		cerr << "First object in the JSON file must be WaveAnalysis. Provided name: " << configName << " Exiting..." << endl;
		exit(EXIT_FAILURE);
	}
	for(int i = 3; i < r; i++){
		CheckJSMNType(t[i], JSMN_STRING, __LINE__, __func__, __FILE__);
		string optName = contents.substr(t[i].start, t[i].end-t[i].start);


		if(i < r-1) i++;
		else continue;

		if(optName == "NumberOfChannels"){
			CheckJSMNType(t[i], JSMN_PRIMITIVE, __LINE__, __func__, __FILE__);
            fNChannels = stoi(contents.substr(t[i].start, t[i].end-t[i].start));
		}
		else if(optName == "VoltageScale"){
			CheckJSMNType(t[i], JSMN_PRIMITIVE, __LINE__, __func__, __FILE__);
            fVoltageScale = stod(contents.substr(t[i].start, t[i].end-t[i].start));
		}
		else if(optName == "Thresholds"){
			CheckJSMNType(t[i], JSMN_ARRAY, __LINE__, __func__, __FILE__);

			for(int j = 0; j < t[i].size; j++){
				CheckJSMNType(t[i+j+1], JSMN_PRIMITIVE, __LINE__, __func__, __FILE__);
                fThresholds.push_back(stod(contents.substr(t[i+j+1].start, t[i+j+1].end-t[i+j+1].start)));
			}
			i += t[i].size;
		}
		else if(optName == "ActiveChannels"){
			CheckJSMNType(t[i], JSMN_ARRAY, __LINE__, __func__, __FILE__);

			for(int j = 0; j < t[i].size; j++){
				CheckJSMNType(t[i+j+1], JSMN_PRIMITIVE, __LINE__, __func__, __FILE__);
                bool flag;
                string val = contents.substr(t[i+j+1].start, t[i+j+1].end-t[i+j+1].start);
                istringstream is(val);
				is >> boolalpha >> flag;
                fActiveChannels.push_back(flag);
			}
			i += t[i].size;
		}
		else if(optName == "ChargeMeasurements"){
			CheckJSMNType(t[i], JSMN_ARRAY, __LINE__, __func__, __FILE__);

			for(int j = 0; j < t[i].size; j++){
				CheckJSMNType(t[i+j+1], JSMN_PRIMITIVE, __LINE__, __func__, __FILE__);
                bool flag;
                string val = contents.substr(t[i+j+1].start, t[i+j+1].end-t[i+j+1].start);
                istringstream is(val);
				is >> boolalpha >> flag;
                fChargeMeasurements.push_back(flag);
			}
			i += t[i].size;
		}
		else if(optName == "TimeMeasurements"){
			CheckJSMNType(t[i], JSMN_ARRAY, __LINE__, __func__, __FILE__);

			for(int j = 0; j < t[i].size; j++){
				CheckJSMNType(t[i+j+1], JSMN_PRIMITIVE, __LINE__, __func__, __FILE__);
                bool flag;
                string val = contents.substr(t[i+j+1].start, t[i+j+1].end-t[i+j+1].start);
                istringstream is(val);
				is >> boolalpha >> flag;
                fTimeMeasurements.push_back(flag);
			}
			i += t[i].size;
		}
		else if(optName == "PedestalWindowLow"){
			CheckJSMNType(t[i], JSMN_ARRAY, __LINE__, __func__, __FILE__);

			for(int j = 0; j < t[i].size; j++){
				CheckJSMNType(t[i+j+1], JSMN_PRIMITIVE, __LINE__, __func__, __FILE__);
                fPedestalWindowLow.push_back(stod(contents.substr(t[i+j+1].start, t[i+j+1].end-t[i+j+1].start)));
			}
			i += t[i].size;
		}
		else if(optName == "PedestalWindowHigh"){
			CheckJSMNType(t[i], JSMN_ARRAY, __LINE__, __func__, __FILE__);

			for(int j = 0; j < t[i].size; j++){
				CheckJSMNType(t[i+j+1], JSMN_PRIMITIVE, __LINE__, __func__, __FILE__);
                fPedestalWindowHigh.push_back(stod(contents.substr(t[i+j+1].start, t[i+j+1].end-t[i+j+1].start)));
			}
			i += t[i].size;
		}
		else if(optName == "AnalysisWindowLow"){
			CheckJSMNType(t[i], JSMN_ARRAY, __LINE__, __func__, __FILE__);

			for(int j = 0; j < t[i].size; j++){
				CheckJSMNType(t[i+j+1], JSMN_PRIMITIVE, __LINE__, __func__, __FILE__);
                fAnalysisWindowLow.push_back(stod(contents.substr(t[i+j+1].start, t[i+j+1].end-t[i+j+1].start)));
			}
			i += t[i].size;
		}
		else if(optName == "AnalysisWindowHigh"){
			CheckJSMNType(t[i], JSMN_ARRAY, __LINE__, __func__, __FILE__);

			for(int j = 0; j < t[i].size; j++){
				CheckJSMNType(t[i+j+1], JSMN_PRIMITIVE, __LINE__, __func__, __FILE__);
                fAnalysisWindowHigh.push_back(stod(contents.substr(t[i+j+1].start, t[i+j+1].end-t[i+j+1].start)));
			}
			i += t[i].size;
		}
		else if(optName == "Polarity"){
			CheckJSMNType(t[i], JSMN_ARRAY, __LINE__, __func__, __FILE__);

			for(int j = 0; j < t[i].size; j++){
				CheckJSMNType(t[i+j+1], JSMN_PRIMITIVE, __LINE__, __func__, __FILE__);
                fPolarity.push_back(stoi(contents.substr(t[i+j+1].start, t[i+j+1].end-t[i+j+1].start)));
			}
			i += t[i].size;
		}
		else if(optName == "comment"){
			//skip
		}

	}


	if((fActiveChannels.size() != fNChannels) || (fChargeMeasurements.size() != fNChannels) || (fTimeMeasurements.size() != fNChannels) ||
        (fPedestalWindowLow.size() != fNChannels) || (fPedestalWindowHigh.size() != fNChannels) || (fAnalysisWindowLow.size() != fNChannels) ||
        (fAnalysisWindowHigh.size() != fNChannels) || (fPolarity.size() != fNChannels)
    ){
        std::cerr << "Error in " << __func__ << ", line " << __LINE__ << " in " << __FILE__ << std::endl;
        std::cerr << "Inconsistent vector size!" << std::endl;
        exit(-1);
    }
}

void AnalysisConfig::Print(){
    cout << "**********************************************************************" << endl;
    cout << "********************** WAVEFORM ANALYSIS CONFIG **********************" << endl;
    cout << "**********************************************************************" << endl;
    cout << "* Number of channels:   " << fNChannels << endl;
    cout << "* Voltage scale:   " << fVoltageScale << endl;
    cout << "* Active channels:      ";
    for(auto i : fActiveChannels)
        printf("%5d, ", (int)i);
    cout << endl;

    cout << "* Thresholds:           ";
    for(auto i : fThresholds)
        printf("%.3f, ", i);
    cout << endl;

    cout << "* Polarity:             ";
    for(auto i : fPolarity)
        printf("%5d, ", i);
    cout << endl;

    cout << "* Charge measurements:  ";
    for(auto i : fChargeMeasurements)
        printf("%5d, ", (int)i);
    cout << endl;

    cout << "* Time measurements:    ";
    for(auto i : fTimeMeasurements)
        printf("%5d, ", (int)i);
    cout << endl;

    cout << "* Pedestal window low:  ";
    for(auto i : fPedestalWindowLow)
        printf("%f, ", i);
    cout << endl;
    cout << "* Pedestal window high: ";
    for(auto i : fPedestalWindowHigh)
        printf("%f, ", i);
    cout << endl;
    cout << endl;
    cout << "* Analysis window low:  ";
    for(auto i : fAnalysisWindowLow)
        printf("%f, ", i);
    cout << endl;
    cout << "* Analysis window high: ";
    for(auto i : fAnalysisWindowHigh)
        printf("%f, ", i);
    cout << endl;
    cout << "**********************************************************************" << endl;
}

void AnalysisConfig::CheckBin(int i){
    if(i < 0 || i >= fNChannels){
        std::cerr << "Error in " << __func__ << ", line " << __LINE__ << " in " << __FILE__ << std::endl;
        std::cerr << "Inconsistent vector size!" << std::endl;
        exit(-1);
    }

}
bool AnalysisConfig::IsActive(int i){
    CheckBin(i);
    return fActiveChannels.at(i);
}

bool AnalysisConfig::MeasureCharge(int i){
    CheckBin(i);
    return fChargeMeasurements.at(i);
}

bool AnalysisConfig::MeasureTime(int i){
    CheckBin(i);
    return fTimeMeasurements.at(i);
}

double AnalysisConfig::GetPedestalBinLow(int i){
    CheckBin(i);
    return fPedestalWindowLow.at(i);
}

double AnalysisConfig::GetPedestalBinHigh(int i){
    CheckBin(i);
    return fPedestalWindowHigh.at(i);
}

double AnalysisConfig::GetAnalysisBinLow(int i){
    CheckBin(i);
    return fAnalysisWindowLow.at(i);
}

double AnalysisConfig::GetAnalysisBinHigh(int i){
    CheckBin(i);
    return fAnalysisWindowHigh.at(i);
}

int AnalysisConfig::GetPolarity(int i){
    CheckBin(i);
    return fPolarity.at(i);
}

double AnalysisConfig::GetThreshold(int i){
    CheckBin(i);
    return fThresholds.at(i);
}
