#ifndef HELPER_LIB_HPP
#define HELPER_LIB_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

void how_to_use();

// class for the reading of the input file
class DataReader {
public:
    DataReader(const std::string& input_file, const std::string& mode) {
        readFile(input_file, mode);
    }

    const std::vector<int>& getRun();
    const std::vector<int>& getEvent();
    const std::vector<int>& getTrigger();
    const std::vector<int>& getIndx();

    const std::vector<double>& getXtrue();
    const std::vector<double>& getYtrue();
    const std::vector<double>& getL1();
    const std::vector<double>& getL2();
    const std::vector<double>& getL3();
    const std::vector<double>& getL4();

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

    void readFile(const std::string& filename, const std::string& mode);
};
#endif /* HELPER_LIB_HPP */
