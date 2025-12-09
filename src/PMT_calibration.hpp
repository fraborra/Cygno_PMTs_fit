//
//  PMT_calibration.hpp
//
//  Created by Stefano Piacentini on 23/09/22.
//  Modified by Francesco Borra on 28/06/23
//
#ifndef PMT_calibration_hpp
#define PMT_calibration_hpp

#include <BAT/BCModel.h>

#include <fstream>
#include <iostream>
#include <omp.h>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

#include "Math/ProbFunc.h"

#include "PMT_association.hpp"
#include "helper_lib.hpp"

class PMTcalibration : public PMTassociation {
  public:
    PMTcalibration() {};
    PMTcalibration(const std::string &mode, int nth, int nPoints,
                   const std::vector<double> &L1_inp,
                   const std::vector<double> &L2_inp,
                   const std::vector<double> &L3_inp,
                   const std::vector<double> &L4_inp,
                   const std::vector<double> &x, const std::vector<double> &y);

    ~PMTcalibration() {};

    double LogLikelihood(const std::vector<double> &pars) override;

  private:
    double Lmax = 40000; // The smaller the smaller the parameter space --> can
                         // be changed if needed if trying to fit higher energy
                         // spot/longer integrals must be modified!

    double cmax =
        20; // The prior for the c_i can be tweaked to reduce parameter space
    std::string mode_;
    unsigned int nPoints;

    // PMT positions (in cm)
    double x1 = 2.3;
    double y1 = 30.7;

    double x2 = 30.7;
    double y2 = 30.7;

    double x3 = 30.7;
    double y3 = 2.3;

    double x4 = 2.3;
    double y4 = 2.3;

    double zGEM = 19;

    std::vector<std::vector<double>> data{4};

    std::vector<double> xTrue;
    std::vector<double> yTrue;
};

void runCalibrationFit(const Config config);

#endif /* PMT_calibration_hpp */
