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
    if (write_chains || write_log || plot) {
        res_dir = "./output_"+mode;

        int com = std::system(("mkdir "+res_dir).c_str());
        if(com != 0) {
            std::cout<<com<<std::endl;
        } else {
            std::cerr << "Failed to create directory: " << res_dir << std::endl;
        }
        res_dir += "/";
        std::cout << res_dir << std::endl;
    }
    
    // Reading input file
    DataReader data(input_file, mode);

    std::vector<int>::size_type index_max = data.getRun().size();
    if (end_ind == -1)
    {
        if(mode.compare("association") == 0) {
            end_ind = index_max;
        }
        else if(mode.compare("PMTcalibration") == 0) {
            end_ind = index_max - index_max%nPoints;
        }

    } else if (index_max<end_ind)
    {
        if(mode.compare("association") == 0) {
            end_ind = index_max;
        }
        else if(mode.compare("PMTcalibration") == 0) {
            end_ind = index_max - index_max%nPoints;
        }
    }

    // prepare output variables
    std::vector<double> L_mean;
    std::vector<double> L_std;
    std::vector<double> x_mean;
    std::vector<double> x_std;
    std::vector<double> y_mean;
    std::vector<double> y_std;

    std::vector<double> c1_mean;
    std::vector<double> c1_std;
    std::vector<double> c2_mean;
    std::vector<double> c2_std;
    std::vector<double> c3_mean;
    std::vector<double> c3_std;
    std::vector<double> c4_mean;
    std::vector<double> c4_std;

// ==================================== START =========================//
    if(mode.compare("association") == 0) { // START OF PMT ASSOCIATION if
        // Setting chains parameters
        int Nch = 6;           //number of parallel MCMC chains
        int NIter = 1*10000;   //number of step per chain

        // BEGIN OF THE FIT LOOP
        for (int index = start_ind; index < end_ind; index++) {
            // import the L from the array
            double L[4] = {0.};
            L[0] = data.getL1()[index];
            L[1] = data.getL2()[index];
            L[2] = data.getL3()[index];
            L[3] = data.getL4()[index];

            // INITIALIZE THE MODEL
            PMTassociation m(mode, Nch, L);

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

// ===============================================================
            // Run MCMC, marginalizing posterior
            m.MarginalizeAll();
// ===============================================================

            // Run mode finding; by default using Minuit
            m.FindMode(m.GetBestFitParameters());

            // Write MCMC on root file (The full chains are not needed for the position reconstruction)
            if (write_chains) {
                m.WriteMarkovChain("prova_mcmc.root", "RECREATE");//, true, true);
                // m.WriteMarkovChain(res_dir+m.GetSafeName()+std::to_string(index) + "_mcmc.root", "RECREATE");
            }

            if (plot) {
                // Draw all marginalized distributions into a PDF file
                m.PrintAllMarginalized(res_dir+m.GetSafeName() + "_" + std::to_string(index) + "_plots.pdf");
            
                // Print summary plots
                m.PrintParameterPlot(res_dir+m.GetSafeName() + "_" + std::to_string(index) + "_parameters.pdf");
                m.PrintCorrelationPlot(res_dir+m.GetSafeName() + "_" + std::to_string(index) + "_correlation.pdf");
                m.PrintCorrelationMatrix(res_dir+m.GetSafeName() + "_" + std::to_string(index) + "_correlationMatrix.pdf");
                m.PrintKnowledgeUpdatePlots(res_dir+m.GetSafeName() + "_" + std::to_string(index) + "_update.pdf");
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
        int NIter = nPoints*40000;   //number of step per chain

        // BEGIN OF THE FIT LOOP
        for (int index = start_ind; index < end_ind; index += nPoints) {

            // prepare helper variables
            std::vector<double> L1;
            std::vector<double> L2;
            std::vector<double> L3;
            std::vector<double> L4;
            std::vector<double> x;
            std::vector<double> y;

            // loop over the nPoints points to fit
            for(int point = 0; point<nPoints; point++){
                L1.push_back(data.getL1()[index+point]);
                L2.push_back(data.getL2()[index+point]);
                L3.push_back(data.getL3()[index+point]);
                L4.push_back(data.getL4()[index+point]);

                x.push_back(data.getXtrue()[index+point]);
                y.push_back(data.getYtrue()[index+point]);
            }

            // INITIALIZE THE MODEL
            PMTcalibration m(mode, Nch, nPoints, L1, L2, L3, L4, x, y);

            // Setting MCMC algorithm and precision
            m.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
            m.SetPrecision(BCEngineMCMC::kMedium);
            if (write_log) {
                BCLog::OutSummary("Model created");
            }

            // Setting prerun iterations to 10^6
            m.SetNIterationsPreRunMax(1000000);
            
            // Setting initial position for the parameters ## NEED TO BE checked
            std::vector<double> x0;
            x0.push_back(1000.);    // L
            for(int i=0; i<nPoints; i++) { // x and y
                x0.push_back(x[i]);
                x0.push_back(y[i]);
            }
            for (int i =0; i<4; i++){// The 4 PMT "calibrations"
                x0.push_back(1.0);
            }

            std::cout << "lunghezza di x0: " << x0.size() <<std::endl;
            std::cout << "Number of parameters: " << m.GetNParameters() << std::endl;
            // m.SetInitialPositions(x0);

            // Setting MC run iterations and number of parallel chains
            m.SetNIterationsRun(NIter);
            m.SetNChains(Nch);
// ===============================================================
            // Run MCMC, marginalizing posterior
            m.MarginalizeAll();
// ===============================================================

            // Run mode finding; by default using Minuit
            m.FindMode(m.GetBestFitParameters());

            // Write MCMC on root file (The full chains are not needed for the position reconstruction)
            if (write_chains) {
                m.WriteMarkovChain("prova_mcmc_cal.root", "RECREATE");//, true, true);
                // m.WriteMarkovChain(res_dir+m.GetSafeName()+std::to_string(index) + "_mcmc.root", "RECREATE");
            }

            if (plot) {
                // Draw all marginalized distributions into a PDF file
                m.PrintAllMarginalized(res_dir+m.GetSafeName() + "_" + std::to_string(index) + "_plots.pdf");
            
                // Print summary plots
                m.PrintParameterPlot(res_dir+m.GetSafeName() + "_" + std::to_string(index) + "_parameters.pdf");
                m.PrintCorrelationPlot(res_dir+m.GetSafeName() + "_" + std::to_string(index) + "_correlation.pdf");
                m.PrintCorrelationMatrix(res_dir+m.GetSafeName() + "_" + std::to_string(index) + "_correlationMatrix.pdf");
                m.PrintKnowledgeUpdatePlots(res_dir+m.GetSafeName() + "_" + std::to_string(index) + "_update.pdf");
            }
            
            // Print results of the analysis
            if (print_summary) {m.PrintSummary();}

            // ========================================================================================================
            // STORE RESULTS
            std::vector<unsigned> H1Indices = m.GetH1DPrintOrder();

            // Check if the pre run has converged:
            int status = m.GetNIterationsConvergenceGlobal();

            if (status>0){ // If prerun converged then store the results
                BCH1D posteriorc1 = m.GetMarginalized(H1Indices[0]);
                BCH1D posteriorc2 = m.GetMarginalized(H1Indices[1]);
                BCH1D posteriorc3 = m.GetMarginalized(H1Indices[2]);
                BCH1D posteriorc4 = m.GetMarginalized(H1Indices[3]);

                c1_mean.push_back(posteriorc1.GetHistogram()->GetMean());
                c1_std.push_back(posteriorc1.GetHistogram()->GetRMS());

                c2_mean.push_back(posteriorc2.GetHistogram()->GetMean());
                c2_std.push_back(posteriorc2.GetHistogram()->GetRMS());

                c3_mean.push_back(posteriorc3.GetHistogram()->GetMean());
                c3_std.push_back(posteriorc3.GetHistogram()->GetRMS());

                c4_mean.push_back(posteriorc4.GetHistogram()->GetMean());
                c4_std.push_back(posteriorc4.GetHistogram()->GetRMS());

                std::cout << "calibration, status > 0" << std::endl;

            } else { // If prerun not converged then store -1 to all the parameters
                c1_mean.push_back(-1);
                c1_std.push_back(-1);
                
                c2_mean.push_back(-1);
                c2_std.push_back(-1);

                c3_mean.push_back(-1);
                c3_std.push_back(-1);
                
                c4_mean.push_back(-1);
                c4_std.push_back(-1);
                std::cout << "calibration, status < 0" << std::endl;
            } // end store results

            if (index % 25 == 0) {// print every 25 index the iteration number
                std::cout << "Iteration number: " << index << std::endl;
            }

        }// end for loop over indices (calibration)
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
        std::vector<int>::size_type control = L_mean.size();
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

    } else if (mode.compare("PMTcalibration") == 0) {
        std::vector<int>::size_type control = c1_mean.size();
        int i = 0;
        for (int index = start_ind; index < end_ind; index += nPoints) // NEED TO CHOOSE THE OUTPUT
        {        
            outfile << run[index] <<"\t"<< event[index] <<"\t"<< trigger[index] <<"\t"<< indx[index] <<"\t"
            << c1_mean[i] <<"\t"<< c1_std[i] <<"\t"  // c1 and std
            << c2_mean[i] <<"\t"<< c2_std[i] <<"\t"  // c2 and std
            << c3_mean[i] <<"\t"<< c3_std[i] <<"\t"  // c3 and std
            << c4_mean[i] <<"\t"<< c4_std[i]         // c4 and std
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
