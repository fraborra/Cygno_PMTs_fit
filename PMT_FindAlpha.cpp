//
//  PMT_FindAlpha.cpp
//  LIMEPMTfits
//
//  Created by Stefano Piacentini on 23/09/22.
//  Modified by Francesco Borra on 28/06/23.
//

#include "PMT_FindAlpha.hpp"
#include <TMath.h>
#include <BAT/BCMath.h>
#include <cmath>


// PMTfit class
PMTfindalpha::PMTfindalpha(const std::string& mode, int nth, int nP, 
                const std::vector<double>& L1_inp, const std::vector<double>& L2_inp, const std::vector<double>& L3_inp, 
                const std::vector<double>& L4_inp, const std::vector<double>&x, const std::vector<double>& y): BCModel(mode)
{
    std::cout<<"Starting fit for '"<<mode<<" reconstruction'"<<std::endl;

    mode_ = mode;
    Lmax = 500000;
    double Lmin = 0.;// if trying to fit higher energy spot/longer integrals must be modified!
    // The prior for the c_i can be tweaked to reduce parameter space
    cmax = 10;
    nPoints = nP;
    
    for (int i = 0; i < nPoints; ++i) {
        data[0].push_back(L1_inp[i]);
        data[1].push_back(L2_inp[i]);
        data[2].push_back(L3_inp[i]);
        data[3].push_back(L4_inp[i]);

        xTrue.push_back(x[i]);
        yTrue.push_back(y[i]);
	std::cout << L1_inp[i] << std::endl;
    }
        
    //DEFINING parameters
     if (mode_.compare("PMTfindalpha") == 0){
        AddParameter("L", Lmin, Lmax, "L", "[a.u.]");
	GetParameter("L").SetPriorConstant();
        
	AddParameter("alpha",2.,6.,"#alpha","");
	GetParameter("alpha").SetPriorConstant();
	
    } else {
        throw std::runtime_error("Unknown model '"+mode_+"'.\n");
    }
    
     omp_set_dynamic(0);
     omp_set_num_threads(nth);
}


// Compute Likelihood
double PMTfindalpha::LogLikelihood(const std::vector<double>& pars) {

    double LL = 0.;
    double c[4] = {1.0,0.965,0.860,0.827};

    // NEED TO RESTORE A LOOP FOR NPOINTS
    for(unsigned int i=0; i<nPoints; i++) { // i == nPoint index

        for(unsigned int j=0; j<4; j++) { // j == PMT index
            double Lij = data[j][i];  // here the data

            double sLij = 0.1*Lij; // for now set to 10% of the integral
            
	    double rij = D2(xTrue[i], yTrue[i], j);     // compute r_i**2
            LL += BCMath::LogGaus(Lij,                            // x, namely Lj
                                 (pars[0]*c[j])/(pow(rij, 0.5*pars[1])), // mu, namely the light computed in the step (c_i * Lj / r_i^4)
                                 sLij,                            // sigma
                                 true                             // norm factor
                                 );
        }
    }

    return LL;
}


// Function to calculate distance between the PMT and the chosen position
double PMTfindalpha::D2(double x, double y, int i) {
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
