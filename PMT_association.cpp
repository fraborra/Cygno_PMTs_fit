//
//  PMT_association.cpp
//  LIMEPMTfits
//
//  Created by Stefano Piacentini on 23/09/22.
//  Modified by Francesco Borra on 28/06/23.
//

#include "PMT_association.hpp"
#include <TMath.h>
#include <BAT/BCMath.h>
#include <cmath>


// PMTfit class
PMTassociation::PMTassociation(const std::string& mode, int nth, 
                               double *L, double *c_tmp) : BCModel(mode)
{
    std::cout<<"Starting fit for '"<<mode<<" reconstruction'"<<std::endl;

    mode_ = mode;
    Lmax = 40000; //The smaller the smaller the parameter space

    for (int i = 0; i < 4; ++i) {
        data[i] = L[i];
        c[i] = c_tmp[i];
    }

    //DEFINING parameters
    if (mode_.compare("association") == 0) {
        AddParameter("L", 0, Lmax, "L", "[a.u.]");
        GetParameter("L").SetPriorConstant();

        AddParameter("x", 0, 33, "x", "[cm]");
        GetParameter("x").SetPriorConstant();

        AddParameter("y", 0, 33, "y", "[cm]");
        GetParameter("y").SetPriorConstant();
        
    } else {
        throw std::runtime_error("Unknown model '"+mode_+"'.\n");
    }
    
    omp_set_dynamic(0);
    omp_set_num_threads(nth);
}


// Compute Likelihood
double PMTassociation::LogLikelihood(const std::vector<double>& pars) {

    double LL = 0.;
    
    for(unsigned int j=0; j<4; j++) {
        double Lj = data[j];  // here the data

        double sLj = 0.1*Lj; // for now set to 10% of the integral
        
        double rij = D2(pars[1], pars[2], j);                 // compute r_i**2
        LL += BCMath::LogGaus(Lj,                             // x, namely Lj
                             (pars[0]*c[j])/(pow(rij, 2)),    // mu, namely the light computed in the step (c_i * Lj / r_i^4)
                             sLj,                             // sigma
                             true                             // norm factor
                             );
    }

    return LL;
}


// Function to calculate distance between the PMT and the chosen position
double PMTassociation::D2(double x, double y, int i) {
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