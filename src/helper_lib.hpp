#ifndef HELPER_LIB_HPP
#define HELPER_LIB_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <stdexcept>

// config structure
struct Config {
    std::string mode;
    std::string input_file;
    int start_ind = 0;
    int end_ind = 0;
    int nPoints = 0;
    std::string output_file = "out.txt";
    std::string calibration_output_tag = "tmp";
    bool plot = false;
    bool write_chains = false;
    bool write_log = false;
    bool print_summary = false;
    double c1 = 1;
    double c2 = 1;
    double c3 = 1;
    double c4 = 1;
};

//  single point result structure
struct SinglePointResult {
    double L_mean;
    double L_std;
    double x_mean;
    double x_std;
    double y_mean;
    double y_std;
};

// function to read config file
Config readConfigFile(const std::string& filename);

// class for the reading of the input file
class DataReader {
public:
    DataReader(const std::string& input_file, const std::string& mode) {
        readFile(input_file, mode);
    }

    const std::vector<int>& getRun(){return run;};
    const std::vector<int>& getEvent(){return event;};
    const std::vector<int>& getTrigger(){return trigger;};
    const std::vector<int>& getIndx(){return indx;};

    const std::vector<double>& getXtrue(){return xtrue;};
    const std::vector<double>& getYtrue(){return ytrue;};
    const std::vector<double>& getL1(){return L1;};
    const std::vector<double>& getL2(){return L2;};
    const std::vector<double>& getL3(){return L3;};
    const std::vector<double>& getL4(){return L4;};

    const std::vector<double>& getSc_integral(){return sc_integral;};

private:
    std::vector<int> run;
    std::vector<int> event;
    std::vector<int> trigger;
    std::vector<int> indx;
    
    std::vector<double> L1;
    std::vector<double> L2;
    std::vector<double> L3;
    std::vector<double> L4;

    std::vector<double> xtrue;
    std::vector<double> ytrue;

    std::vector<double> sc_integral;

    void readFile(const std::string& filename, const std::string& mode);
};
#endif /* HELPER_LIB_HPP */
