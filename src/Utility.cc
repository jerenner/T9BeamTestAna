#include "Utility.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

void utl::SplitFilename (const std::string& str, std::string &file, std::string &path)
{
	std::cout << "Splitting: " << str << '\n';
	std::size_t found = str.find_last_of("/\\");
	path = str.substr(0,found);
	file = str.substr(found+1);
	std::cout << " path: " << path << '\n';
	std::cout << " file: " << file << '\n';
}

std::string utl::ConvertToString(double x, int prec){
	std::stringstream ss;
	std::string val;

	ss << std::fixed << std::setprecision(prec);
	ss << x;
	ss >> val;

	return val;
}


bool utl::CheckFile(std::string name){
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } 
    else {
        return false;
    }   
}

void utl::ReadInputList(std::string fileName, std::vector<std::string> &job){

        std::ifstream ifs;
        std::string line;
        ifs.open(fileName.c_str());


        if (!ifs.is_open()){
                std::cerr << "Error in " << __func__ << ", line " << __LINE__ << " in " << __FILE__ << std::endl;
                std::cerr << "Provided list file does not exist!" << std::endl;
                std::cerr << "Exiting..." << std::endl;
                exit (EXIT_FAILURE);
        }
        while (getline(ifs, line)){
                if(line.at(0)=='#')
                        continue;

                if(!utl::CheckFile(line)){
                        std::cerr << "Warning in " << __func__ << ", line " << __LINE__ << " in " << __FILE__ << std::endl;
                        std::cerr << "File " << line << " does not exist! Skipping..." << std::endl;
                        continue;
                }
                job.push_back(line);


        }
        ifs.close();

}
