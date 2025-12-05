//
//  PMT_singleEvent.hpp
//  Francesco Borra, Dic 2025
//

#ifndef PMT_singleEvent_hpp
#define PMT_singleEvent_hpp

#include <BAT/BCModel.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <omp.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include "Math/ProbFunc.h"
#include "PMT_association.hpp"
#include "helper_lib.hpp"


class PMTSingleEvent : public PMTassociation
{
public:

    PMTSingleEvent(const std::string& mode, int nth, double *L, double *c_tmp);
    
    ~PMTSingleEvent(){};

//     double LogLikelihood(const std::vector<double>& pars);

    // double D2(double x, double y, int i);
    
    // double EvaluateSigma(double mu);
    
private:
    double Lmax = 40000; // The smaller the smaller the parameter space --> can be changed if needed
                         // if trying to fit higher energy spot/longer integrals must be modified!
    double c[4] = {1.};

    std::string mode_;
        
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
};

SinglePointResult runSingleEvent(const double q_values[4],
                                 const double calib[4],
                                 int NIterPrerun,
                                 int NIter,
                                 bool second_round);

void saveSingleEventResults(const SinglePointResult results);

#endif /* PMT_singleEvent_hpp */
