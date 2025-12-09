//
//  PMT_association.cpp
//
//  Created by Stefano Piacentini on 23/09/22.
//  Modified by Francesco Borra on 28/06/23.
//

#include "PMT_association.hpp"
#include "helper_lib.hpp"
#include <BAT/BCMath.h>
#include <TMath.h>
#include <cmath>

// PMTfit class
PMTassociation::PMTassociation(const std::string &mode, int nth, double *L,
                               double *c_tmp)
    : BCModel(mode) {
    std::cout << "Starting fit for '" << mode << " reconstruction'"
              << std::endl;

    mode_ = mode;

    for (int i = 0; i < 4; ++i) {
        data[i] = L[i];
        c[i] = c_tmp[i];
    }

    // DEFINING parameters
    if (mode_ == "association") {
        AddParameter("L", 0, Lmax, "L", "[a.u.]");
        GetParameter("L").SetPriorConstant();

        AddParameter("x", 0, 33, "x", "[cm]");
        GetParameter("x").SetPriorConstant();

        AddParameter("y", 0, 33, "y", "[cm]");
        GetParameter("y").SetPriorConstant();

    } else {
        throw std::runtime_error("Unknown model '" + mode_ + "'.\n");
    }

    omp_set_dynamic(0);
    omp_set_num_threads(nth);
}

// Compute Likelihood
double PMTassociation::LogLikelihood(const std::vector<double> &pars) {

    double LL = 0.;

    for (unsigned int j = 0; j < 4; j++) {
        double Lj = data[j]; // here the data

        double rij = D2(pars[1], pars[2], j);         // compute r_i**2
        double mu = (pars[0] * c[j]) / (pow(rij, 2)); // compute mu

        double sLj = EvaluateSigma(mu); // compute uncertainty on Lj

        LL += BCMath::LogGaus(
            Lj,  // x, namely Lj
            mu,  // mu, namely the light computed in the step (c_i * Lj / r_i^4)
            sLj, // sigma
            true // norm factor
        );
    }

    return LL;
}

// Evaluate sigma
double PMTassociation::EvaluateSigma(double mu) {
    // func: sigma = [0]*sqrt(x)+x*[1]
    double par0 = 0.02;
    double par1 = 0.06;

    return par0 * sqrt(mu) + par1 * mu;
}

// Function to calculate distance between the PMT and the chosen position
double PMTassociation::D2(double x, double y, int i) {
    if (i == 0) {
        return (x - x1) * (x - x1) + (y - y1) * (y - y1) + zGEM * zGEM;
    } else if (i == 1) {
        return (x - x2) * (x - x2) + (y - y2) * (y - y2) + zGEM * zGEM;
    } else if (i == 2) {
        return (x - x3) * (x - x3) + (y - y3) * (y - y3) + zGEM * zGEM;
    } else if (i == 3) {
        return (x - x4) * (x - x4) + (y - y4) * (y - y4) + zGEM * zGEM;
    } else {
        throw std::runtime_error("Uknown value of PMT index.\n");
    }
    return 0.;
}

// run fit function
void runAssociationFit(const Config config) {

    std::string mode = config.mode;
    std::string input_file = config.input_file;
    int start_ind = config.start_ind;
    int end_ind = config.end_ind;
    std::string output_file = config.output_file;
    bool plot = config.plot;
    bool write_chains = config.write_chains;
    bool write_log = config.write_log;
    bool print_summary = config.print_summary;
    int nPoints = config.nPoints;
    double c1 = 1;
    double c2 = config.c2 / config.c1;
    double c3 = config.c3 / config.c1;
    double c4 = config.c4 / config.c1;

    double c_list[4] = {c1, c2, c3, c4};

    //  Create results folder if not existing and plot/log/chains requested
    std::string res_dir;
    res_dir = "./data/output_" + mode;
    if (plot || write_chains || write_log || print_summary) {
        int com = std::system(("mkdir -p " + res_dir).c_str());

        if (com == 0) {
            std::cerr << "Failed to create directory: " << res_dir << std::endl;
        }
        res_dir += "/";
    }

    DataReader data(input_file, mode);

    int index_max = static_cast<int>(data.getRun().size());
    if (end_ind == -1 || index_max < end_ind) {
        end_ind = index_max;
    }
    // If write chains it will fit max 5 point to avoid too much data
    if (write_chains) {
        end_ind = std::min(end_ind, 5);
    }

    // prepare output variables
    std::vector<double> L_mean;
    std::vector<double> L_std;
    std::vector<double> x_mean;
    std::vector<double> x_std;
    std::vector<double> y_mean;
    std::vector<double> y_std;
    std::vector<double> corrLx;
    std::vector<double> corrLy;
    std::vector<double> corrxy;

    // ==================================== START =========================//
    // Setting chains parameters
    int Nch = 6;                 // number of parallel MCMC chains
    int NIter = nPoints * 10000; // number of step per chain

    // BEGIN OF THE FIT LOOP
    for (int index = start_ind; index < end_ind; index++) {
        // import the L from the array
        double L[4] = {0.};
        L[0] = data.getL1()[index];
        L[1] = data.getL2()[index];
        L[2] = data.getL3()[index];
        L[3] = data.getL4()[index];

        // Create log
        if (write_log) {
            BCLog::OpenLog("./logs/association_log.txt", BCLog::detail,
                           BCLog::detail);
        }

        // INITIALIZE THE MODEL
        PMTassociation m(mode, Nch, L, c_list);

        // Setting MCMC algorithm and precision
        m.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
        m.SetPrecision(BCEngineMCMC::kMedium);
        if (write_log) {
            BCLog::OutSummary("Model created");
        }
        // Setting prerun iterations to 10^5 (for fast integration, if it does
        // not converge it is saved as not converged)
        m.SetNIterationsPreRunMax(100000);

        // Setting MC run iterations and number of parallel chains
        m.SetNIterationsRun(NIter);
        m.SetNChains(Nch);

        // Name for BAT outputs
        std::string BAT_out_prefix =
            res_dir + m.GetSafeName() + "_" + std::to_string(index);

        // Write MCMC on root file (The full chains are not needed for the
        // position reconstruction)
        if (write_chains) {
            m.WriteMarkovChain(BAT_out_prefix + "_mcmc.root", "RECREATE");
        }

        // ===============================================================
        // Run MCMC, marginalizing posterior
        m.MarginalizeAll();
        // ===============================================================

        // Run mode finding; by default using Minuit
        m.FindMode(m.GetBestFitParameters());

        if (plot) {
            // Draw all marginalized distributions into a PDF file
            m.PrintAllMarginalized(BAT_out_prefix + "_plots.pdf");

            // Print summary plots
            m.PrintParameterPlot(BAT_out_prefix + "_parameters.pdf");
            m.PrintCorrelationPlot(BAT_out_prefix + "_correlation.pdf");
            m.PrintCorrelationMatrix(BAT_out_prefix + "_correlationMatrix.pdf");
            m.PrintKnowledgeUpdatePlots(BAT_out_prefix + "_update.pdf");
        }

        // Print results of the analysis
        if (print_summary) {
            m.PrintSummary();
        }

        //==================
        // RESULTS
        std::vector<unsigned> H1Indices = m.GetH1DPrintOrder();

        // Check if the pre run has converged:
        int status = m.GetNIterationsConvergenceGlobal();

        // start results storing
        if (status > 0) { // If prerun converged then store the results
            BCH1D posteriorL = m.GetMarginalized(H1Indices[0]);
            BCH1D posteriorx = m.GetMarginalized(H1Indices[1]);
            BCH1D posteriory = m.GetMarginalized(H1Indices[2]);

            L_mean.push_back(posteriorL.GetHistogram()->GetMean());
            L_std.push_back(posteriorL.GetHistogram()->GetRMS());

            x_mean.push_back(posteriorx.GetHistogram()->GetMean());
            x_std.push_back(posteriorx.GetHistogram()->GetRMS());

            y_mean.push_back(posteriory.GetHistogram()->GetMean());
            y_std.push_back(posteriory.GetHistogram()->GetRMS());

            BCH2D postLx = m.GetMarginalized(H1Indices[0], H1Indices[1]);
            BCH2D postLy = m.GetMarginalized(H1Indices[0], H1Indices[2]);
            BCH2D postxy = m.GetMarginalized(H1Indices[1], H1Indices[2]);

            corrLx.push_back(postLx.GetHistogram()->GetCorrelationFactor());
            corrLy.push_back(postLy.GetHistogram()->GetCorrelationFactor());
            corrxy.push_back(postxy.GetHistogram()->GetCorrelationFactor());

        } else { // If prerun not converged then store -1 to all the parameters
            L_mean.push_back(-1);
            L_std.push_back(-1);

            x_mean.push_back(-1);
            x_std.push_back(-1);

            y_mean.push_back(-1);
            y_std.push_back(-1);

            corrLx.push_back(-1);
            corrLx.push_back(-1);
            corrLx.push_back(-1);

            std::cout << "association, status < 0" << std::endl;
        }

        if (index % 100 == 0) { // print every 100 index the iteration number
            std::cout << "Iteration number: " << index << std::endl;
        }
    } // end for loop over row indices

    // PRINT RESULTS ON FILE
    std::ofstream outfile;
    outfile.open(output_file, std::ios_base::trunc);

    std::vector<int> run = data.getRun();
    std::vector<int> event = data.getEvent();
    std::vector<int> trigger = data.getTrigger();
    std::vector<int> indx = data.getIndx();

    int control = static_cast<int>(L_mean.size());
    int i = 0;
    for (int index = start_ind; index < end_ind; index++) {
        outfile << run[index] << "\t" << event[index] << "\t" << trigger[index]
                << "\t" << indx[index] << "\t" << L_mean[i] << "\t" << L_std[i]
                << "\t"                                  // L and L_std
                << x_mean[i] << "\t" << x_std[i] << "\t" // x and x_std
                << y_mean[i] << "\t" << y_std[i] << "\t" // y and y_std
                << corrLx[i] << "\t" << corrLy[i] << "\t"
                << corrxy[i] // Correlations

                << std::endl;
        i++;
        if (i >= control) {
            break;
        }
    }

    outfile.close();

    if (write_log) {
        // Close log file
        BCLog::OutSummary("Exiting");
        BCLog::CloseLog();
    }
}
