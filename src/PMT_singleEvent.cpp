// 
//  PMTSingleEvent.cpp
//  Francesco Borra, Dic 2025
//  
#include <TMath.h>
#include <BAT/BCMath.h>
#include <cmath>
#include "PMT_association.hpp"
#include "PMT_singleEvent.hpp"
#include "helper_lib.hpp"
// #include <nlohmann/json.hpp>


// PMTSingleEvent class
PMTSingleEvent::PMTSingleEvent(const std::string& mode, int nth, 
                               double *L, double *c_tmp) : PMTassociation()
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


// // Compute Likelihood
// double PMTSingleEvent::LogLikelihood(const std::vector<double>& pars) {

//     double LL = 0.;
    
//     for(unsigned int j=0; j<4; j++) {
//         double Lj = data[j];  // here the data
        
//         double rij = D2(pars[1], pars[2], j);          // compute r_i**2
//         double mu = (pars[0]*c[j])/(pow(rij, 2));     // compute mu

// //         double sLj = 0.1*Lj; // for now set to 10% of the integral        
//         double sLj = EvaluateSigma(mu);        // compute uncertainty on Lj
        
//         LL += BCMath::LogGaus(Lj,               // x, namely Lj
//                              mu,                // mu, namely the light computed in the step (c_i * Lj / r_i^4)
//                              sLj,               // sigma
//                              true               // norm factor
//                              );
//     }

//     return LL;
// }

// // Evaluate sigma
// double PMTSingleEvent::EvaluateSigma(double mu) {
//     //func: sigma = [0]*sqrt(x)+x*[1]
//     double par0 = 0.02;
//     double par1 = 0.06;

//     return par0*sqrt(mu)+par1*mu;
// }

// // Function to calculate distance between the PMT and the chosen position
// double PMTSingleEvent::D2(double x, double y, int i) {
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

// run fit function
SinglePointResult runSingleEvent(
                    const double q_values[4],
                    const double calib[4],
                    int NIterPrerun = 100000,
                    int NIter = 10000,
                    bool second_round = false
                    ) 
{

    std::string mode = "association";
    double c1 = calib[0];
    double c2 = calib[1]/c1;
    double c3 = calib[2]/c1;
    double c4 = calib[3]/c1;
    
    double c_list[4] = {c1, c2, c3, c4};


    // ==================================== START =========================//
    // Setting chains parameters
    int Nch = 6;           //number of parallel MCMC chains
    // int NIter = 10000;   //number of step per chain

    double L[4] = {0.};
    L[0] = q_values[0];
    L[1] = q_values[1];
    L[2] = q_values[2];
    L[3] = q_values[3];
    
    // INITIALIZE THE MODEL
    PMTSingleEvent m(mode, Nch, L, c_list);

    // Setting MCMC algorithm and precision
    m.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
    m.SetPrecision(BCEngineMCMC::kMedium);

    // Setting prerun iterations to 10^5 (for fast integration, if it does not converge it is saved as not converged)
    m.SetNIterationsPreRunMax(NIterPrerun);

    // Setting MC run iterations and number of parallel chains
    m.SetNIterationsRun(NIter);
    m.SetNChains(Nch);

// ===============================================================
    // Run MCMC, marginalizing posterior
    m.MarginalizeAll();
// ===============================================================

    // Run mode finding; by default using Minuit
    m.FindMode(m.GetBestFitParameters());

    std::vector<unsigned> H1Indices = m.GetH1DPrintOrder();

    // Check if the pre run has converged:
    int status = m.GetNIterationsConvergenceGlobal();

    SinglePointResult results;
    // start results storing
    if (status>0){ // If prerun converged then store the results
        BCH1D posteriorL = m.GetMarginalized(H1Indices[0]);
        BCH1D posteriorx = m.GetMarginalized(H1Indices[1]);
        BCH1D posteriory = m.GetMarginalized(H1Indices[2]);
        
        double L_mean = posteriorL.GetHistogram()->GetMean();
        double L_std  = posteriorL.GetHistogram()->GetRMS();
        
        double x_mean = posteriorx.GetHistogram()->GetMean();
        double x_std  = posteriorx.GetHistogram()->GetRMS();
        
        double y_mean = posteriory.GetHistogram()->GetMean();
        double y_std  = posteriory.GetHistogram()->GetRMS();
        
        BCH2D postLx = m.GetMarginalized(H1Indices[0], H1Indices[1]);
        BCH2D postLy = m.GetMarginalized(H1Indices[0], H1Indices[2]);
        BCH2D postxy = m.GetMarginalized(H1Indices[1], H1Indices[2]);
        
//         double corrLx = postLx.GetHistogram()->GetCorrelationFactor();
//         double corrLy = postLy.GetHistogram()->GetCorrelationFactor();
//         double corrxy = postxy.GetHistogram()->GetCorrelationFactor();

        results = {L_mean, L_std, x_mean, x_std, y_mean, y_std};

    } else if (!second_round) { // If prerun not converged try again with more iterations and more prerun iterations

        NIter = 50000;   //number of step per chain
        NIterPrerun = 500000; //number of prerun iterations

        results = runSingleEvent(
                    q_values,
                    calib,
                    NIterPrerun,
                    NIter,
                    second_round = true
                    );

    } else {    // If prerun not converged again then store -1 to all the parameters
    
//         double L_mean = -1.;   
//         double L_std  = -1.;
        
//         double x_mean = -1.;
//         double x_std  = -1.;
        
//         double y_mean = -1.;
//         double y_std  = -1.;
        
//         double corrLx = -1.;
//         double corrLx = -1.;
//         double corrLx = -1.;
        std::cout << "association, status < 0" << std::endl;
    
        results = {-1, -1, -1, -1, -1, -1};

    }

    return results;

}

// void saveSingleEventResults(const SinglePointResult results) 
// {
//     std::string res_dir = "./output_singleEvent";
//     int com = std::system(("mkdir -p "+res_dir).c_str());
        
//     if(com == 0) {
//     std::cerr << "Failed to create directory: " << res_dir << std::endl;
//     }
//     res_dir += "/";

//     const std::string& output_file = res_dir+"/results.json";

//     // save results to json
//     nlohmann::json j;

//     j["L_mean"] = results.L_mean;
//     j["L_std"]  = results.L_std;
//     j["x_mean"] = results.x_mean;
//     j["x_std"]  = results.x_std;
//     j["y_mean"] = results.y_mean;
//     j["y_std"]  = results.y_std;

//     std::ofstream(output_file) << j.dump(4); // 4 = indentation
// }
