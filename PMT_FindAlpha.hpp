//
//  PMT_FindAlpha.cpp
//  LIMEPMTfits
//
//  Created by Stefano Piacentini on 23/09/22.
//  Modified by Francesco Borra on 28/06/23
//
#ifndef PMT_FindAlpha_hpp
#define PMT_FindAlpha_hpp

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

class PMTfindalpha : public BCModel
{
public:

    PMTfindalpha(const std::string& mode, int nth, int nPoints, 
                   const std::vector<double>& L1_inp, const std::vector<double>& L2_inp, const std::vector<double>& L3_inp, 
                   const std::vector<double>& L4_inp, const std::vector<double>&x, const std::vector<double>& y);

    ~PMTfindalpha(){};

    double LogLikelihood(const std::vector<double>& pars);

    double LogAPrioriProbability(const std::vector<double>& pars);

    double D2(double x, double y, int i);
    
private:
    double Lmax;
    double cmax;
    std::string mode_;
    unsigned int nPoints;

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

    std::vector<std::vector<double>> data {4};

    std::vector<double> xTrue;
    std::vector<double> yTrue;

};

#endif /* PMT_FindAlpha_hpp */
