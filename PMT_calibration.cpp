//
//  PMT_calibration.cpp
//  LIMEPMTfits
//
//  Created by Stefano Piacentini on 23/09/22.
//  Modified by Francesco Borra on 28/06/23.
//

#include "PMT_calibration.hpp"
#include <TMath.h>
#include <BAT/BCMath.h>
#include <cmath>


// PMTfit class
PMTcalibration::PMTcalibration(const std::string& mode, int nth, int nP, 
                const std::vector<double>& L1_inp, const std::vector<double>& L2_inp, const std::vector<double>& L3_inp, 
                const std::vector<double>& L4_inp, const std::vector<double>&x, const std::vector<double>& y): BCModel(mode)
{
    std::cout<<"Starting fit for '"<<mode<<" reconstruction'"<<std::endl;

    mode_ = mode;
    Lmax = 40000; // if trying to fit higher energy spot/longer integrals must be modified!
    // The prior for the c_i can be tweaked to reduce parameter space
    cmax = 10;
    nPoints = nP;
    // std::cout << "Dentro il fit" << std::endl << std::endl;

    for (int i = 0; i < nPoints; ++i) {
        data[0].push_back(L1_inp[i]);
        data[1].push_back(L2_inp[i]);
        data[2].push_back(L3_inp[i]);
        data[3].push_back(L4_inp[i]);

        xTrue.push_back(x[i]);
        yTrue.push_back(y[i]);
    }
        
    //DEFINING parameters
     if (mode_.compare("PMTcalibration") == 0){
        AddParameter("L", 0, Lmax, "L", "[a.u.]");
        GetParameter("L").Fix(4000.0); // just to have c_i values smaller, can put any value, 
                                       // we are only interested in the c_i ratios
        
        for(unsigned int i=0; i < nPoints; i++) {
            AddParameter("x_"+std::to_string(i), 0, 33, "x_"+std::to_string(i), "[cm]");
            AddParameter("y_"+std::to_string(i), 0, 33, "y_"+std::to_string(i), "[cm]");
            //  FIXING x and y COORDINATES         
            GetParameter("x_"+std::to_string(i)).Fix(xTrue[i]);
            GetParameter("y_"+std::to_string(i)).Fix(yTrue[i]);
        }

        // The prior for the c_i can be tweaked to reduce parameter space
        AddParameter("c1", 0., cmax, "c1", "[counts]");
        AddParameter("c2", 0., cmax, "c2", "[counts]");
        AddParameter("c3", 0., cmax, "c3", "[counts]");
        AddParameter("c4", 0., cmax, "c4", "[counts]");

    } else {
        throw std::runtime_error("Unknown model '"+mode_+"'.\n");
    }
    
    omp_set_dynamic(0);
    omp_set_num_threads(nth);
}


// Compute Likelihood
double PMTcalibration::LogLikelihood(const std::vector<double>& pars) {

    double LL = 0.;
    
    // NEED TO RESTORE A LOOP FOR NPOINTS
    for(unsigned int i=0; i<nPoints; i++) { // i == nPoint index

        for(unsigned int j=0; j<4; j++) { // j == PMT index
            double Lij = data[j][i];  // here the data

            double sLij = 0.1*Lij; // for now set to 10% of the integral
            
            // int k = 1+ 2*nPoints +j;  // index for c_i (pars[3] is c_1 if nPoints=1 and so on)
            int k = pars.size()-4 + j;

            double tmp = D2(pars[1 + 2*i], pars[2 + 2*i], j);     // compute r_i**2
            LL += BCMath::LogGaus(Lij,                            // x, namely Lj
                                 (pars[0]*pars[k])/(pow(tmp, 2)), // mu, namely the light computed in the step (c_i * Lj / r_i^4)
                                 sLij,                            // sigma
                                 true                             // norm factor
                                 );
        }
    }

    return LL;
}

// Define prior
double PMTcalibration::LogAPrioriProbability(const std::vector<double>& pars) {
    double LL = 0.;
        //flat priors everywhere

    // L prior
    if(pars[0]<0 || pars[0]>Lmax) {
        LL += log(0.0);
    } else {
        LL += log(1.0/Lmax);
    }
    // x_i and y_i prior
    for(unsigned int i=1; i < 1 + nPoints*2; i++) {
        if(pars[i]<0 || pars[i]>33.) {
            LL += log(0.0);
        } else {
            LL += log(1.0/33.);
        }
    }
    // c_i prior
    for(unsigned int j=(pars.size()-4); j<pars.size(); j++){
            if(pars[j]<0. || pars[j]>cmax) {
                LL += log(0.0);
            } else {
                LL += log(1.0/cmax);
            }
        }
    return LL;
}

// Function to calculate distance between the PMT and the chosen position
double PMTcalibration::D2(double x, double y, int i) {
    if (i == 0) {
        return (x - x1)*(x - x1) + (y - y1)*(y - y1) + zGEM*zGEM;
    } else if (i == 1) {
        return (x - x2)*(x - x2) + (y - y2)*(y - y2) + zGEM*zGEM;
    } else if (i == 2) {
        return (x - x3)*(x - x3) + (y - y3)*(y - y3) + zGEM*zGEM;
    } else if (i == 3) {
        return (x - x4)*(x - x4) + (y - y4)*(y - y4) + zGEM*zGEM;
    } else {
        throw std::runtime_error("Uknown value of PMT index.\n");
    }
    return 0.;
}