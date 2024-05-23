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

    PMTfit(const std::string& name, std::string model, std::string datafile,
           int nth, int Np, int Iprec);

    ~PMTfit(){};

    double LogLikelihood(const std::vector<double>& pars);

    double LogAPrioriProbability(const std::vector<double>& pars);

    double D2(double x, double y, int i);
    
//    std::vector<int> run;
    int getRun(int i);
//    std::vector<int> event;
    int getEvent(int i);
//    std::vector<int> trigger;
    int getTrigger(int i);
//
//    std::vector<int> mj2;
    int getMj2(int i);
//    std::vector<int> GE;
    int getGE(int i);
//    std::vector<int> indx;
    int getIndx(int i);
//
//    std::vector<double> xtrue;
    double getXtrue(int i);
//    std::vector<double> ytrue;
    double getYtrue(int i);
//    std::vector<double> sX;
    double getSigX(int i);
//    std::vector<double> integral;
    double getIntegral(int i);
//     std::vector<int> time;

private:
    double Lmax;
    unsigned int Npoints;
    std::string r_type_;
    std::string res_dir_;

    // prior parameters
    double L_mean = 0;
    double L_std = 0;
    
    double c_mean[4] = {0.};
    double c_std[4] = {0.};
    
    int iprec_;

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

    std::vector<std::vector<double>> data {4};

    // Imported values
    std::vector<int> run;
    std::vector<int> event;
    std::vector<int> trigger;
    
    std::vector<int> mj2;
    std::vector<int> GE;
    std::vector<int> indx;
    
    std::vector<double> xtrue;
    std::vector<double> ytrue;
    std::vector<double> sX;
    std::vector<double> integral;

    
    int ReadInputFile(std::string filename);
        
};

void print_usage() {
    std::cout << "Usage: program -m mode -i input_file -r start_row -o output_file [-p] [-c] [-l] [-h]" << std::endl;
}


#endif /* PMT_standard_hpp */
