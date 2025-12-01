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
#include <nlohmann/json.hpp>

#include "PMT_association.hpp"
#include "PMT_calibration.hpp"
#include "PMT_FindAlpha.hpp"
#include "PMT_singleEvent.hpp"
#include "helper_lib.hpp"

int main(int argc, char *argv[]) {

    if (argc < 2) {
        std::cerr << "Usage:\n";
        std::cerr << "  " << argv[0] << " <config_file>\n";
        std::cerr << "  " << argv[0] << " --single-point c1 c2 c3 c4\n";
        return 1;
    }

    // std::string configFile = argv[1];
    std::string firstArg = argv[1];
    if (firstArg == "--single-point") {
        // single point mode
        if (argc < 6) { // 4 values needed (the 4 c_i)
            std::cerr << "Error: need 4 values for single point mode\n";
            return 1;
        }

        double values[4];
        for (int i = 2; i < argc; ++i) {
            values[i] = (std::stod(argv[i]));
        }

        // run single point fit
        std::array<double, 9> results;
        results = fitSinglePoint(values);
        // print results as json
        saveSingleEventResults(results);
        
    } else {
        // config file mode
        std::string configFile = argv[1];
        Config config = readConfigFile(configFile);
        
        // check mode validity
        std::string mode = config.mode;
        const std::string known_mode[3] = {"association", "PMTcalibration", "PMTfindalpha"};

        bool compare_r_type = false;
        for(int i=0; i<3;i++){
            if(mode.compare(known_mode[i]) == 0) {
                compare_r_type = true;
            }
        }
        if(compare_r_type) {
            std::cout<<"Initialization of '"<<mode<<" reconstruction'..."<<std::endl;
        } else {
            throw std::invalid_argument("Unknown reconstruction type '"+mode+"'");
        }
        
        // if for mode
        if(mode.compare("association") == 0) {
            std::cout<<"Running PMT association reconstruction..."<<std::endl;
            runAssociationFit(config);

        } else if(mode.compare("PMTcalibration") == 0) {
            std::cout<<"Running PMT calibration..."<<std::endl;
            runCalibrationFit(config);

        } else if(mode.compare("PMTfindalpha") == 0) {
            std::cout<<"Running PMT find alpha..."<<std::endl;
            runFindAlphaFit(config);

        } else { throw std::invalid_argument("Unknown reconstruction type '"+mode+"'");}

    }
    
    return 0;
}
