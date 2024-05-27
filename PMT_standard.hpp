//
//  PMT_standard.hpp
//  LIMEPMTfits
//
//  Created by Stefano Piacentini on 23/09/22.
//  Modified by Francesco Borra on 28/06/23
//
#ifndef PMT_standard_hpp
#define PMT_standard_hpp

#include <BAT/BCModel.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <omp.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

#include "Math/ProbFunc.h"

class PMTfit : public BCModel
{
public:

    PMTfit(const std::string& mode, int nth, int index, double *L, double x, double y);

    ~PMTfit(){};

    double LogLikelihood(const std::vector<double>& pars);

    double LogAPrioriProbability(const std::vector<double>& pars);

    double D2(double x, double y, int i);
    
private:
    double Lmax;
    double cmax;
    std::string mode_;
    std::int index_;


    // prior parameters
    double L_mean = 0;
    double L_std = 0;
    
    double c_mean[4] = {0.};
    double c_std[4] = {0.};
    
    //PMT positions (in cm)
    double x1 = 2.3;
    double y1 = 30.7;

    double x2 = 30.7;
    double y2 = 30.7;

    double x3 = 30.7;
    double y3 = 2.3;

    double x4 = 2.3;
    double y4 = 2.3;

    double zGEM = 19;

    double data[4] = {0.};

    double xTrue = 0.;
    double yTrue = 0.;

};

void print_usage() {
    std::cout << "Usage: program -m mode -i input_file -s start_ind -e end_ind -o output_file [-p] [-c] [-l] [-h]" << std::endl;
}


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



#endif /* PMT_standard_hpp */
