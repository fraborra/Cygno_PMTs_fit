//
//  PMT_calibration.hpp
//  LIMEPMTfits
//
//  Created by Stefano Piacentini on 23/09/22.
//  Modified by Francesco Borra on 28/06/23
//
#ifndef PMT_calibration_hpp
#define PMT_calibration_hpp

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

class PMTcalibration : public BCModel
{
public:

    PMTcalibration(const std::string& mode, int nth, int index, 
                   double *L1_inp, double *L2_inp, double *L3_inp, double *L4_inp, 
                   double x, double y);

    ~PMTcalibration(){};

    double LogLikelihood(const std::vector<double>& pars);

    double LogAPrioriProbability(const std::vector<double>& pars);

    double D2(double x, double y, int i);
    
private:
    double Lmax;
    double cmax;
    std::string mode_;
    int index_;


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

    double L1[4] = {0.};
    double L2[4] = {0.};
    double L3[4] = {0.};
    double L4[4] = {0.};

    double xTrue = 0.;
    double yTrue = 0.;

};

#endif /* PMT_calibration_hpp */