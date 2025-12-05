//
//  PMT_singleEvent.cpp
//  Francesco Borra, Dic 2025
//
#include "PMT_singleEvent.hpp"
#include "PMT_association.hpp"
#include "helper_lib.hpp"
#include <BAT/BCMath.h>
#include <TMath.h>
#include <cmath>

// PMTSingleEvent class
// PMTSingleEvent::PMTSingleEvent(const std::string &mode, int nth, double *L,
//                                double *c_tmp)
//     : PMTassociation() {
//   std::cout << "Starting fit for '" << mode << " reconstruction'" <<
//   std::endl;

//   mode_ = mode;

//   for (int i = 0; i < 4; ++i) {
//     data[i] = L[i];
//     c[i] = c_tmp[i];
//   }

//   // DEFINING parameters
//   if (mode_ == "association") {
//     AddParameter("L", 0, Lmax, "L", "[a.u.]");
//     GetParameter("L").SetPriorConstant();

//     AddParameter("x", 0, 33, "x", "[cm]");
//     GetParameter("x").SetPriorConstant();

//     AddParameter("y", 0, 33, "y", "[cm]");
//     GetParameter("y").SetPriorConstant();

//   } else {
//     throw std::runtime_error("Unknown model '" + mode_ + "'.\n");
//   }

//   omp_set_dynamic(0);
//   omp_set_num_threads(nth);
// }

// run fit function
SinglePointResult runSingleEvent(const double q_values[4],
                                 const double calib[4],
                                 int NIterPrerun = 100000, int NIter = 10000,
                                 bool second_round = false) {

  std::string mode = "association";
  double c1 = calib[0];
  double c2 = calib[1] / c1;
  double c3 = calib[2] / c1;
  double c4 = calib[3] / c1;

  double c_list[4] = {c1, c2, c3, c4};

  // ==================================== START =========================//
  // Setting chains parameters
  int Nch = 6; // number of parallel MCMC chains
  // int NIter = 10000;   //number of step per chain

  double L[4] = {0.};
  L[0] = q_values[0];
  L[1] = q_values[1];
  L[2] = q_values[2];
  L[3] = q_values[3];

  // INITIALIZE THE MODEL
  // PMTSingleEvent m(mode, Nch, L, c_list);
  // testing association class
  PMTassociation m(mode, Nch, L, c_list);

  // Setting MCMC algorithm and precision
  m.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
  m.SetPrecision(BCEngineMCMC::kMedium);

  // Setting prerun iterations to 10^5 (for fast integration, if it does not
  // converge it is saved as not converged)
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
  if (status > 0) { // If prerun converged then store the results
    BCH1D posteriorL = m.GetMarginalized(H1Indices[0]);
    BCH1D posteriorx = m.GetMarginalized(H1Indices[1]);
    BCH1D posteriory = m.GetMarginalized(H1Indices[2]);

    double L_mean = posteriorL.GetHistogram()->GetMean();
    double L_std = posteriorL.GetHistogram()->GetRMS();

    double x_mean = posteriorx.GetHistogram()->GetMean();
    double x_std = posteriorx.GetHistogram()->GetRMS();

    double y_mean = posteriory.GetHistogram()->GetMean();
    double y_std = posteriory.GetHistogram()->GetRMS();

    BCH2D postLx = m.GetMarginalized(H1Indices[0], H1Indices[1]);
    BCH2D postLy = m.GetMarginalized(H1Indices[0], H1Indices[2]);
    BCH2D postxy = m.GetMarginalized(H1Indices[1], H1Indices[2]);

    results = {L_mean, L_std, x_mean, x_std, y_mean, y_std};

  } else if (!second_round) { // If prerun not converged try again with more
                              // iterations and more prerun iterations

    NIter = 50000;        // number of step per chain
    NIterPrerun = 500000; // number of prerun iterations

    results = runSingleEvent(q_values, calib, NIterPrerun, NIter,
                             second_round = true);

  } else { // If prerun not converged again then store -1 to all the parameters

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
