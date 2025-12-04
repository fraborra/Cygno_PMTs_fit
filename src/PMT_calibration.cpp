//
//  PMT_calibration.cpp
//
//  Created by Stefano Piacentini on 23/09/22.
//  Modified by Francesco Borra on 28/06/23.
//

#include <TMath.h>
#include <BAT/BCMath.h>
#include <cmath>
#include "PMT_calibration.hpp"
#include "PMT_association.hpp"
#include "helper_lib.hpp"

// PMTfit class
PMTcalibration::PMTcalibration(const std::string& mode, int nth, int nP, 
                const std::vector<double>& L1_inp, const std::vector<double>& L2_inp, const std::vector<double>& L3_inp, 
                const std::vector<double>& L4_inp, const std::vector<double>&x, const std::vector<double>& y) : PMTassociation()
{
    std::cout<<"Starting fit for '"<<mode<<" reconstruction'"<<std::endl;

    mode_ = mode;
    Lmax = 40000; // if trying to fit higher energy spot/longer integrals must be modified!
    // The prior for the c_i can be tweaked to reduce parameter space
    cmax = 20;
    nPoints = nP;
    // std::cout << "Dentro il fit" << std::endl << std::endl;

    for (unsigned int i = 0; i < nPoints; ++i) {
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

        // The prior for the c_i can be tweaked to reduce parameter space
        AddParameter("c1", 0., cmax, "c1", "[counts]");
        GetParameter("c1").SetPriorConstant();
        
        AddParameter("c2", 0., cmax, "c2", "[counts]");
        GetParameter("c2").SetPriorConstant();
        
        AddParameter("c3", 0., cmax, "c3", "[counts]");
        GetParameter("c3").SetPriorConstant();
        
        AddParameter("c4", 0., cmax, "c4", "[counts]");
        GetParameter("c4").SetPriorConstant();

    } else {
        throw std::runtime_error("Unknown model '"+mode_+"'.\n");
    }
    
    omp_set_dynamic(0);
    omp_set_num_threads(nth);
}


// Compute Likelihood
double PMTcalibration::LogLikelihood(const std::vector<double>& pars) {

    double LL = 0.;
    
    for(unsigned int i=0; i<nPoints; i++) { // i == nPoint index
        
        for(unsigned int j=0; j<4; j++) { // j == PMT index

            double Lij = data[j][i];  // here the data

            // int k = 1+j;  // index for c_i (pars[1] is c_1)
            int k = 1 + j;
            
            double rij = D2(xTrue[i], yTrue[i], j);               // compute r_ij**2
            double mu = (pars[0]*pars[k])/(pow(rij, 2));

            double sLij = EvaluateSigma(mu); //

            LL += BCMath::LogGaus(Lij,      // x, namely Lj
                                 mu,        // mu, namely the light computed in the step (c_i * Lij / r_ij^4)
                                 sLij,      // sigma
                                 true       // norm factor
                                 );
        }
    }

    return LL;
}

// // Evaluate sigma
// double PMTcalibration::EvaluateSigma(double mu) {
//     //func: sigma = [0]*sqrt(x)+x*[1]
//     double par0 = 0.02;
//     double par1 = 0.06;

//     return par0*sqrt(mu)+par1*mu;
// }



// // Function to calculate distance between the PMT and the chosen position
// double PMTcalibration::D2(double x, double y, int i) {
//     if (i == 0) {
//         return (x - x1)*(x - x1) + (y - y1)*(y - y1) + zGEM*zGEM;
//     } else if (i == 1) {
//         return (x - x2)*(x - x2) + (y - y2)*(y - y2) + zGEM*zGEM;
//     } else if (i == 2) {
//         return (x - x3)*(x - x3) + (y - y3)*(y - y3) + zGEM*zGEM;
//     } else if (i == 3) {
//         return (x - x4)*(x - x4) + (y - y4)*(y - y4) + zGEM*zGEM;
//     } else {
//         throw std::runtime_error("Uknown value of PMT index.\n");
//     }
//     return 0.;
// }

void runCalibrationFit(const Config config){

    // Setting chains parameters
    int Nch = 12;          //number of parallel MCMC chains
    int NIter = 100000;    //number of step per chain

    // prepare helper variables
    std::vector<double> L1;
    std::vector<double> L2;
    std::vector<double> L3;
    std::vector<double> L4;
    std::vector<double> x;
    std::vector<double> y;

    std::string mode = config.mode;
    std::string input_file = config.input_file;
    int start_ind = config.start_ind;
    int end_ind = config.end_ind;
    int nPoints = config.nPoints;
    std::string output_tag = config.calibration_output_tag;
    bool plot = config.plot;

    // read data
    DataReader data(input_file, mode);
    
    int index_max = static_cast<int>(data.getRun().size());
    if (end_ind == -1 || index_max<end_ind) {end_ind = index_max;}


    //  Create results folder if not existing
    std::string res_dir;
    res_dir = "./data/output_"+mode;
    int com = std::system(("mkdir -p "+res_dir).c_str());
    if(com == 0) {
    std::cerr << "Failed to create directory: " << res_dir << std::endl;
    }
    res_dir += "/";

    // loop over the nPoints points to store data
    for(int point = 0; point<end_ind; point++){
        // HERE I NORMALIZE THE Li USING SC_INTEGRAL, SO THAT IF THE LY IS NOT CONSTANT IS ALL GOOD
        L1.push_back(data.getL1()[point]/data.getSc_integral()[point]*10000);
        L2.push_back(data.getL2()[point]/data.getSc_integral()[point]*10000);
        L3.push_back(data.getL3()[point]/data.getSc_integral()[point]*10000);
        L4.push_back(data.getL4()[point]/data.getSc_integral()[point]*10000);

        x.push_back(data.getXtrue()[point]);
        y.push_back(data.getYtrue()[point]);
    }
        
    // Create log
    BCLog::OpenLog("./logs/calibration_log.txt", BCLog::detail, BCLog::detail);

    // INITIALIZE THE MODEL
    PMTcalibration cal(mode, Nch, end_ind, L1, L2, L3, L4, x, y);

    // Setting MCMC algorithm and precision
    cal.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
    cal.SetPrecision(BCEngineMCMC::kMedium);

    BCLog::OutSummary("Model created");

    // Setting prerun iterations to 10^6
    cal.SetNIterationsPreRunMax(1000000);
    
    // Setting MC run iterations and number of parallel chains
    cal.SetNIterationsRun(NIter);
    cal.SetNChains(Nch);
    
    // Prefix for BAT outputs
    std::string BAT_out_prefix_cal = res_dir+cal.GetSafeName() + "_" + output_tag;

    cal.WriteMarkovChain(BAT_out_prefix_cal + "_mcmc.root", "RECREATE");

// ===============================================================
    // Run MCMC, marginalizing posterior
    cal.MarginalizeAll();
// ===============================================================

    // Run mode finding; by default using Minuit
    cal.FindMode(cal.GetBestFitParameters());

    if (plot) {
        // Draw all marginalized distributions into a PDF file
        cal.PrintAllMarginalized(BAT_out_prefix_cal + "_plots.pdf");
    
        // Print summary plots
        cal.PrintParameterPlot(BAT_out_prefix_cal + "_parameters.pdf");
        cal.PrintCorrelationPlot(BAT_out_prefix_cal + "_correlation.pdf");
        cal.PrintCorrelationMatrix(BAT_out_prefix_cal + "_correlationMatrix.pdf");
        cal.PrintKnowledgeUpdatePlots(BAT_out_prefix_cal + "_update.pdf");
    }
    
    // Print results of the analysis
    cal.PrintSummary();

    // Close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

}