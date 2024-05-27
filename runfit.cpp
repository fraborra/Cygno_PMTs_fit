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
// #include <filesystem>

#include "PMT_standard.hpp"

namespace fs = std::filesystem;


int main(int argc, char *argv[]) {

    int option;
    std::string mode;
    std::string input_file;
    std::string res_dir;
    int start_ind = -1;
    int end_ind = -1;
    std::string output_file;
    bool plot = false;
    bool write_chains = false;
    bool write_log = false;
    
    // Parsing arguments
    while ((option = getopt(argc, argv, "hm:i:s:e:o:pcl")) != -1) {
        switch (option) {
            case 'h':
                print_usage();
                return 0;
            case 'm':
                mode = optarg;
                break;
            case 'i':
                input_file = optarg;
                break;
            case 's':
                start_ind = std::stoi(optarg);
                break;
            case 'e':
                end_ind = std::stoi(optarg);
                break;
            case 'o':
                output_file = optarg;
                break;
            case 'p':
                plot = true;
                break;
            case 'c':
                write_chains = true;
                break;
            case 'l':
                write_log = true;
                break;
            default:
                print_usage();
                return 1;
        }
    }
    
    if (mode.empty() || input_file.empty() || start_ind == -1 || output_file.empty() || end_ind == -1) {
        std::cerr << "Missing required options." << std::endl;
        print_usage();
        return 1;
    }

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
        if (!fs::exists(res_dir)) {
            if (fs::create_directory(res_dir)) {
                std::cout << "Directory created successfully: " << res_dir << std::endl;
            } else {
                std::cerr << "Failed to create directory: " << res_dir << std::endl;
            }
        } else {
            std::cout << "Directory already exists: " << res_dir << std::endl;
        }
        res_dir += "/";
    }
    
// Setting chains parameters
    
    int Nch, NIter;
    if(mode.compare("association") == 0) {
        Nch = 6;           //number of parallel MCMC chains
        NIter = 1*10000;   //number of step per chain
    } else if (mode.compare("PMTcalibration") == 0){
        Nch = 12;          //number of parallel MCMC chains
        NIter = 1*10000;   //number of step per chain
    }

    // Reading input file
    DataReader data(input_file);
    

    // BEGIN OF THE FIT LOOP
    for (int index = start_ind; index < end_ind+1; index++)
    {
        double x, y;
        // import the L from the array
        double L[4] = {0.};
        L[0] = data.getL1()[index];
        L[1] = data.getL2()[index];
        L[2] = data.getL3()[index];
        L[3] = data.getL4()[index];

        if(mode.compare("association") == 0) {
            x = 0.;
            y = 0.;
        } else if (mode.compare("PMTcalibration") == 0)
        {
            x = data.getXtrue()[index];
            y = data.getYtrue()[index];
        }
        

        // initialize the model
        PMTfit m(mode, Nch, index, L, x, y);

        // Creating Logfile
        if (write_log) {
            BCLog::OpenLog(res_dir + m.GetSafeName() + "_" + index + "_log.txt", BCLog::detail, BCLog::detail);
        }

        // Setting MCMC algorithm and precision
        m.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
        m.SetPrecision(BCEngineMCMC::kMedium);
        if (write_log) {
            BCLog::OutSummary("Model created");
        }
        // Setting prerun iterations to 10^5 (for fast integration, if it does not converge it is saved as not converged)
        if(mode.compare("association") == 0) {
            m.SetNIterationsPreRunMax(100000);
        }
        
        if(mode.compare("PMTcalibration") == 0) {
            m.SetNIterationsPreRunMax(1000000);
            std::vector<double> x0;
            x0.push_back(100.);    // L
            for(int i=0; i<2; i++) { // x and y
                x0.push_back(16.5);
                x0.push_back(16.5);
            }
            for (int i =0; i<4; i++){
                x0.push_back(1.0); // The 4 PMT "calibrations"
            }
            m.SetInitialPositions(x0);
        }
        
        // Setting MC run iterations and number of parallel chains
        m.SetNIterationsRun(NIter);
        m.SetNChains(Nch);



        // Write MCMC on root file (The full chains are not needed for the position reconstruction)
        if (write_chains) {
            m.WriteMarkovChain(res_dir+m.GetSafeName() + "_" + index + "_mcmc.root", "RECREATE");
        }

        // Run MCMC, marginalizing posterior
        m.MarginalizeAll();

        // Run mode finding; by default using Minuit
        m.FindMode(m.GetBestFitParameters());

        if (plot) {
            // Draw all marginalized distributions into a PDF file
            m.PrintAllMarginalized(res_dir+m.GetSafeName() + "_" + index + "_plots.pdf");
        
            // Print summary plots
            m.PrintParameterPlot(res_dir+m.GetSafeName() + "_" + index + "_parameters.pdf");
            m.PrintCorrelationPlot(res_dir+m.GetSafeName() + "_" + index + "_correlation.pdf");
            m.PrintCorrelationMatrix(res_dir+m.GetSafeName() + "_" + index + "_correlationMatrix.pdf");
            m.PrintKnowledgeUpdatePlots(res_dir+m.GetSafeName() + "_" + index + "_update.pdf");
        }
        
        // Print results of the analysis into a text file
        m.PrintSummary();

        // Need to add output

    }
    
        
    
    // std::vector<unsigned> H1Indices = m.GetH1DPrintOrder();
            
    // std::ofstream outfile;
    // outfile.open(outfile_name, std::ios_base::app);

    // // Check if the pre run has converged:
    // int status = m.GetNIterationsConvergenceGlobal();

    // // print results on file
    // if(mode.compare("association") == 0) {
    //     if (status>0){
    //         BCH1D posteriorL = m.GetMarginalized(H1Indices[0]);
    //         BCH1D posteriorx1 = m.GetMarginalized(H1Indices[1]);
    //         BCH1D posteriory1 = m.GetMarginalized(H1Indices[2]);

    //         outfile<<m.getRun(start_ind)<<"\t" <<m.getEvent(start_ind)<<"\t"<<m.getTrigger(start_ind)<<"\t"
    //         <<m.getGE(start_ind)<<"\t" <<m.getIndx(start_ind)<<"\t"
    //         <<posteriorL.GetHistogram()->GetMean()<< "\t"<<posteriorL.GetHistogram()->GetRMS()<< "\t"      // L and L_std
    //         <<(posteriorx1.GetHistogram())->GetMean()<<"\t"<<(posteriorx1.GetHistogram())->GetRMS()<< "\t" // x and x_std
    //         <<(posteriory1.GetHistogram())->GetMean()<<"\t"<<(posteriory1.GetHistogram())->GetRMS()<< "\t" // y and y_std
    //         <<m.getMj2(start_ind) << std::endl;
    //     } else {
    //         outfile<<m.getRun(start_ind)<<"\t" <<m.getEvent(start_ind)<<"\t"<<m.getTrigger(start_ind)<<"\t"
    //         <<m.getGE(start_ind)<<"\t" <<m.getIndx(start_ind)<<"\t"
    //         <<"-100"<< "\t"<<"-100"<< "\t" // L and L_std
    //         <<"-100"<< "\t"<<"-100"<< "\t" // x and x_std
    //         <<"-100"<< "\t"<<"-100"<< "\t" // y and y_std
    //         <<m.getMj2(start_ind)<< std::endl;
            
    //         std::cout<<"==================================================="<<std::endl
    //         << std::endl<< std::endl<< "Did not converge"
    //         <<"==================================================="<<std::endl<<std::endl;
    //     }
        
    // } else if (mode.compare("PMTcalibration") == 0) {
    //     if (status>0){
    //         BCH1D posteriorx1 = m.GetMarginalized(H1Indices[0]);
    //         BCH1D posteriory1 = m.GetMarginalized(H1Indices[1]);
            
    //         BCH1D posc1 = m.GetMarginalized(H1Indices[2]);
    //         BCH1D posc2 = m.GetMarginalized(H1Indices[3]);
    //         BCH1D posc3 = m.GetMarginalized(H1Indices[4]);
    //         BCH1D posc4 = m.GetMarginalized(H1Indices[5]);

    //         outfile<<m.getRun(start_ind)<<"\t"
    //         <<(posteriorx1.GetHistogram())->GetMean()<<"\t"<<(posteriorx1.GetHistogram())->GetRMS()<< "\t" // x1 and x1_std
    //         <<(posteriory1.GetHistogram())->GetMean()<<"\t"<<(posteriory1.GetHistogram())->GetRMS()<< "\t" // y1 and y1_std

    //         <<(posc1.GetHistogram())->GetMean()<<"\t"<<(posc1.GetHistogram())->GetRMS()<< "\t"      // c1 and c1_std
    //         <<(posc2.GetHistogram())->GetMean()<<"\t"<<(posc2.GetHistogram())->GetRMS()<< "\t"      // c2 and c2_std
    //         <<(posc3.GetHistogram())->GetMean()<<"\t"<<(posc3.GetHistogram())->GetRMS()<< "\t"      // c3 and c3_std
    //         <<(posc4.GetHistogram())->GetMean()<<"\t"<<(posc4.GetHistogram())->GetRMS()<<std::endl; // c4 and c4_std
            
    //     } else {
    //         outfile<<m.getRun(start_ind)<<"\t"
    //         <<"-100"<< "\t"<<"-100"<< "\t" <<"-100"<< "\t"<<"-100"<< "\t" // x1, y1 with stds

    //         <<"-100"<< "\t"<<"-100"<< "\t"       // c1 and c1_std
    //         <<"-100"<< "\t"<<"-100"<< "\t"       // c2 and c2_std
    //         <<"-100"<< "\t"<<"-100"<< "\t"       // c3 and c3_std
    //         <<"-100"<< "\t"<<"-100"<< std::endl; // c4 and c4_std
            
    //         std::cout<<"==================================================="<<std::endl
    //         << std::endl<< std::endl<< "Did not converge"
    //         <<"==================================================="<<std::endl<<std::endl;
    //     }
        
    // } else if {
    //     throw std::runtime_error("Unknown reconstruction '"+mode+"'.\n");
    // }
    
    // if (write_log) {
    //     // Close log file
    //     BCLog::OutSummary("Exiting");
    //     BCLog::CloseLog();
    // }
    return 0;
}
