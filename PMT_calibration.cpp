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
PMTcalibration::PMTcalibration(const std::string& mode, int nth, 
               int index, double *L1_inp, double *L2_inp, double *L3_inp, double *L4_inp, 
               double x, double y) : BCModel(mode)
{
    std::cout<<"Starting fit for '"<<mode<<" reconstruction'"<<std::endl;

    mode_ = mode;
    Lmax = 20000;
    cmax = 2;
    index_ = index;
    xTrue = x;
    yTrue = y;

    // std::cout << "Dentro il fit" << std::endl << std::endl;

    for (int i = 0; i < 4; ++i) {

        L1[i] = L1_inp[i];
        L2[i] = L2_inp[i];
        L3[i] = L3_inp[i];
        L4[i] = L4_inp[i];
        // std::cout << "L"<<i<<":\t" << L[i] << "\t";
    }
        // std::cout << std::endl;
        
    //DEFINING parameters
     if (mode_.compare("PMTcalibration") == 0){
        AddParameter("L", 0, Lmax, "L", "[a.u.]");
        GetParameter("L").Fix(4000.0); // just to have c_i values smaller, can put any value, 
                                       // we are only interested in the c_i ratios
        AddParameter("x", 0, 33, "x", "[cm]");
        AddParameter("y", 0, 33, "y", "[cm]"); 
        //  FIXING x and y COORDINATES         
        GetParameter("x").Fix(xTrue);
        GetParameter("y").Fix(yTrue);

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
    
    for(unsigned int j=0; j<4; j++) {
        double Lj = data[j];  // here the data
        // std::cout<< "i: " << j << "\t" << "L: "<< Lj << std::endl;

        double sLj = 0.1*Lj; // for now set to 10% of the integral
        
        int k = 3 +j;  // index for c_i (pars[3] is c_1 and so on)
        // std::cout<< "k: " << k << "\t" << "cj: "<< pars[k] << std::endl;
        std::cout.flush();
        double tmp = D2(pars[1], pars[2], j);                // compute r_i**2
        LL += BCMath::LogGaus(Lj,                             // x, namely Lj
                             (pars[0]*pars[k])/(pow(tmp, 2)), // mu, namely the light computed in the step (c_i * Lj / r_i^4)
                             sLj,                             // sigma
                             true                             // norm factor
                             );
    }

    return LL;
}

// Define prior
double PMTcalibration::LogAPrioriProbability(const std::vector<double>& pars) {
    double LL = 0.;
        //flat priors everywhere
        if(pars[0]<0 || pars[0]>Lmax) {
            LL += log(0.0);
        } else {
            LL += log(1.0/Lmax);
        }
    
    for(unsigned int i=1; i<3; i++) {
        if(pars[i]<0 || pars[i]>33.) {
            LL += log(0.0);
        } else {
            LL += log(1.0/33.);
        }
    }
    
    for(unsigned int j=3; j<7; j++){
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