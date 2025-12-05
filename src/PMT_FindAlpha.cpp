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

#include "PMT_calibration.hpp"
#include "PMT_association.hpp"
#include "helper_lib.hpp"


// PMTfindalpha class
PMTfindalpha::PMTfindalpha(const std::string& mode, int nth, int nP, 
                const std::vector<double>& L1_inp, const std::vector<double>& L2_inp, const std::vector<double>& L3_inp, 
                const std::vector<double>& L4_inp, const std::vector<double>&x, const std::vector<double>& y, double *c_tmp) : PMTcalibration()
{
    std::cout<<"Starting fit for '"<<mode<<" reconstruction'"<<std::endl;

    mode_ = mode;
    double Lmin = 0.;
    for (int i = 0; i < 4; ++i) {
        c[i] = c_tmp[i];
    }

    nPoints = nP;
    
    for (unsigned int i = 0; i < nPoints; ++i) {
        data[0].push_back(L1_inp[i]);
        data[1].push_back(L2_inp[i]);
        data[2].push_back(L3_inp[i]);
        data[3].push_back(L4_inp[i]);

        xTrue.push_back(x[i]);
        yTrue.push_back(y[i]);
    }
        
    //DEFINING parameters
    if (mode_ == "PMTfindalpha"){
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

    for(unsigned int i=0; i<nPoints; i++) { // i == nPoint index

        for(int j=0; j<4; j++) { // j == PMT index

            double Lij = data[j][i];  // here the data

            
            double rij = D2(xTrue[i], yTrue[i], j);             // compute r_i**2
            double mu = (pars[0]*c[j])/(pow(rij, 0.5*pars[1])); // compute mu

            double sLij = EvaluateSigma(mu);        // compute uncertainty on Lj

            LL += BCMath::LogGaus(Lij,              // x, namely Lj
                                    mu,             // mu, namely the light computed in the step (c_i * Lj / r_i^4)
                                    sLij,           // sigma
                                    true            // norm factor
                                    );
        }
    }

    return LL;
}

void runFindAlphaFit(const Config config){
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
    std::string output_tag = config.calibration_output_tag;
    bool plot = config.plot;
    double c1 = 1;
    double c2 = config.c2/config.c1;
    double c3 = config.c3/config.c1;
    double c4 = config.c4/config.c1;
    double c_list[4] = {c1, c2, c3, c4};


    //  Create results folder if not existing
    std::string res_dir;
    res_dir = "./output_"+mode;
    int com = std::system(("mkdir -p "+res_dir).c_str());
    if(com == 0) {
    std::cerr << "Failed to create directory: " << res_dir << std::endl;
    }
    res_dir += "/";
    
    // Read data
    DataReader data(input_file, mode);
    
    int index_max = static_cast<int>(data.getRun().size());
    if (end_ind == -1 || index_max<end_ind) {end_ind = index_max;}

    // store data
    for(int point = start_ind; point<end_ind; point++){
        // HERE I NORMALIZE THE Li USING SC_INTEGRAL, SO THAT IF THE LY IS NOT CONSTANT IS ALL GOOD
        L1.push_back(data.getL1()[point]/data.getSc_integral()[point]*10000);
        L2.push_back(data.getL2()[point]/data.getSc_integral()[point]*10000);
        L3.push_back(data.getL3()[point]/data.getSc_integral()[point]*10000);
        L4.push_back(data.getL4()[point]/data.getSc_integral()[point]*10000);
        
        x.push_back(data.getXtrue()[point]);
        y.push_back(data.getYtrue()[point]);
    }

    // Create log
    BCLog::OpenLog("./logs/findalpha_log.txt", BCLog::detail, BCLog::detail);

    // INITIALIZE THE MODEL
    PMTfindalpha alpha(mode, Nch, end_ind, L1, L2, L3, L4, x, y, c_list);

    // Setting MCMC algorithm and precision
    alpha.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
    alpha.SetPrecision(BCEngineMCMC::kMedium);

    BCLog::OutSummary("Model created");

    // Setting prerun iterations to 10^6
    alpha.SetNIterationsPreRunMax(1000000);
    
    // Setting MC run iterations and number of parallel chains
    alpha.SetNIterationsRun(NIter);
    alpha.SetNChains(Nch);
    
    // Prefix for BAT outputs
    std::string BAT_out_prefix_alpha = res_dir+alpha.GetSafeName() + "_" + output_tag;

    alpha.WriteMarkovChain(BAT_out_prefix_alpha + "_mcmc.root", "RECREATE");

    // ===============================================================
    // Run MCMC, marginalizing posterior
    alpha.MarginalizeAll();
    // ===============================================================

    // Run mode finding; by default using Minuit
    alpha.FindMode(alpha.GetBestFitParameters());
    
    if (plot) {
        // Draw all marginalized distributions into a PDF file
        alpha.PrintAllMarginalized(BAT_out_prefix_alpha + "_plots.pdf");
            
        // Print summary plots
        alpha.PrintParameterPlot(BAT_out_prefix_alpha + "_parameters.pdf");
        alpha.PrintCorrelationPlot(BAT_out_prefix_alpha + "_correlation.pdf");
        alpha.PrintCorrelationMatrix(BAT_out_prefix_alpha + "_correlationMatrix.pdf");
        alpha.PrintKnowledgeUpdatePlots(BAT_out_prefix_alpha + "_update.pdf");
    }
    
    // Print results of the analysis
    alpha.PrintSummary();
    
    // Close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

}

