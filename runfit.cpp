//
//  runfit.cpp
//  LIMEPMTfits
//
//  Created by Stefano Piacentini on 23/09/22.
//  Modified by Francesco Borra on 28/06/23.
//

#include <BAT/BCLog.h>

#include <iostream>
#include <string>
#include <stdexcept>
#include <omp.h>

#include "PMT_association.hpp"
#include "PMT_calibration.hpp"
#include "helper_lib.hpp"

int main(int argc, char *argv[]) {

    if (argc < 2) {
        std::cerr << "Error: you need to specify the config file." << std::endl;
        std::cerr << "Usage: " << argv[0] << " <config_file>" << std::endl;
        return 1;
    }

    std::string configFile = argv[1];
    Config config = readConfigFile(configFile);

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
    double c2 = config.c2/config.c1;
    double c3 = config.c3/config.c1;
    double c4 = config.c4/config.c1;

    double c_list[4] = {c1, c2, c3, c4};


    std::string res_dir;

    const std::string known_mode[2] = {"association", "PMTcalibration"};
                                        // MAYBE ADD ONE THAT ONLY GIVES THE CHAINS AND POSTERIORS
    
    bool compare_r_type = false;
    for(int i=0; i<2;i++){
        if(mode.compare(known_mode[i]) == 0) {
            compare_r_type = true;
        }
    }
    if(compare_r_type) {
        std::cout<<"Initialization of '"<<mode<<" reconstruction'..."<<std::endl;
    } else {
        throw std::invalid_argument("Unknown reconstruction type '"+mode+"'");
    }
    
//  Create results folder if not existing and plot/log/chains requested
    res_dir = "./output_"+mode;

    int com = std::system(("mkdir "+res_dir).c_str());
    if(com ==0) {
        std::cerr << "Failed to create directory: " << res_dir << std::endl;
    }
    res_dir += "/";
    
    // Reading input file
    DataReader data(input_file, mode);

    int index_max = static_cast<int>(data.getRun().size());
    if (end_ind == -1 || index_max<end_ind) {end_ind = index_max;}

    // prepare output variables
    std::vector<double> L_mean;
    std::vector<double> L_std;
    std::vector<double> x_mean;
    std::vector<double> x_std;
    std::vector<double> y_mean;
    std::vector<double> y_std;

// ==================================== START =========================//
    if(mode.compare("association") == 0) { // START OF PMT ASSOCIATION if
        // Setting chains parameters
        int Nch = 6;           //number of parallel MCMC chains
        int NIter = nPoints*10000;   //number of step per chain

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
                BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);
            }

            // INITIALIZE THE MODEL
            PMTassociation m(mode, Nch, L, c_list);

            // Setting MCMC algorithm and precision
            m.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
            m.SetPrecision(BCEngineMCMC::kMedium);
            if (write_log) {
                BCLog::OutSummary("Model created");
            }
            // Setting prerun iterations to 10^5 (for fast integration, if it does not converge it is saved as not converged)
            m.SetNIterationsPreRunMax(100000);

            // Setting MC run iterations and number of parallel chains
            m.SetNIterationsRun(NIter);
            m.SetNChains(Nch);

            // Name for BAT outputs
            std::string BAT_out_prefix = res_dir+m.GetSafeName() + "_" + input_file + "_" +std::to_string(index);

            // Write MCMC on root file (The full chains are not needed for the position reconstruction)
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
            if (print_summary) {m.PrintSummary();}

            //==================
            // RESULTS
            std::vector<unsigned> H1Indices = m.GetH1DPrintOrder();

            // Check if the pre run has converged:
            int status = m.GetNIterationsConvergenceGlobal();

            // start results storing
            if (status>0){ // If prerun converged then store the results
                BCH1D posteriorL = m.GetMarginalized(H1Indices[0]);
                BCH1D posteriorx = m.GetMarginalized(H1Indices[1]);
                BCH1D posteriory = m.GetMarginalized(H1Indices[2]);

                L_mean.push_back(posteriorL.GetHistogram()->GetMean());
                L_std.push_back(posteriorL.GetHistogram()->GetRMS());
                
                x_mean.push_back(posteriorx.GetHistogram()->GetMean());
                x_std.push_back(posteriorx.GetHistogram()->GetRMS());
                
                y_mean.push_back(posteriory.GetHistogram()->GetMean());
                y_std.push_back(posteriory.GetHistogram()->GetRMS());

            } else { // If prerun not converged then store -1 to all the parameters
                L_mean.push_back(-1);
                L_std.push_back(-1);
                
                x_mean.push_back(-1);
                x_std.push_back(-1);
                
                y_mean.push_back(-1);
                y_std.push_back(-1);

                std::cout << "association, status < 0" << std::endl;
            }

            if (index % 100 == 0) {// print every 100 index the iteration number
                std::cout << "Iteration number: " << index << std::endl;
            }
        } // end for loop over row indices
    } // END PMT ASSOCIATION if



    // START PMT CALIBRATION IF
    else if (mode.compare("PMTcalibration") == 0) { 

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
        BCLog::OpenLog("calibration_log.txt", BCLog::detail, BCLog::detail);

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
        std::string BAT_out_prefix_cal = res_dir+cal.GetSafeName() + "_" + input_file;

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
        
        return 0;


        // }// end for loop over indices (calibration)
    } // END PMT CALIBRATION if
    else { throw std::invalid_argument("Unknown reconstruction type '"+mode+"'");}



    // PRINT RESULTS ON FILE
    std::ofstream outfile;
    outfile.open(output_file, std::ios_base::trunc);

    std::vector<int> run = data.getRun();
    std::vector<int> event = data.getEvent();
    std::vector<int> trigger = data.getTrigger();
    std::vector<int> indx = data.getIndx();

    if(mode.compare("association") == 0) {
        int control = static_cast<int>(L_mean.size());
        int i = 0;
        for (int index = start_ind; index < end_ind; index++)
        {        
            outfile << run[index] <<"\t"<< event[index] <<"\t"<< trigger[index] <<"\t"<< indx[index] <<"\t"
            << L_mean[i] <<"\t"<< L_std[i] <<"\t"  // L and L_std
            << x_mean[i] <<"\t"<< x_std[i] <<"\t"  // x and x_std
            << y_mean[i] <<"\t"<< y_std[i]         // y and y_std
            << std::endl;
            i++;
            if (i >= control) {
                break;
            }

        }     
    }

    outfile.close();

    if (write_log) {
        // Close log file
        BCLog::OutSummary("Exiting");
        BCLog::CloseLog();
    }

    return 0;
}
