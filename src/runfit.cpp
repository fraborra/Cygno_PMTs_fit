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
// #include <nlohmann/json.hpp>

#include "PMT_association.hpp"
#include "PMT_calibration.hpp"
#include "PMT_FindAlpha.hpp"
#include "PMT_singleEvent.hpp"
#include "helper_lib.hpp"

int main(int argc, char *argv[]) {

    if (argc < 2) {
        std::cerr << "Usage:\n";
        std::cerr << "  " << argv[0] << " <config_file>\n";
        return 1;
    }

    // std::string configFile = argv[1];
//     std::string firstArg = argv[1];
//     if (firstArg == "--single-point") {
//         // single point mode
//         if (argc < 6) { // 4 values needed (the 4 c_i)
//             std::cerr << "Error: need 4 values for single point mode\n";
//             return 1;
//         }

//         double values[4];
//         for (int i = 2; i < argc; ++i) {
//             values[i] = (std::stod(argv[i]));
//         }
        
//         // create config object for singleEvent:
//         Config config;
//         config.mode = 'association';
//         config.c1 
        
//         // run single point fit
//         SinglePointResult results = runSingleEvent(values);
//         // print results as json
//         saveSingleEventResults(results);
        
//     } else {
        // config file mode
        std::string configFile = argv[1];
        Config config = readConfigFile(configFile);
        
        // check mode validity
        std::string mode = config.mode;
        
        // if for mode
        if(mode == "association") {
            std::cout<<"Running PMT association reconstruction..."<<std::endl;
            runAssociationFit(config);

        } else if(mode == "PMTcalibration") {
            std::cout<<"Running PMT calibration..."<<std::endl;
            runCalibrationFit(config);

        } else if(mode == "PMTfindalpha") {
            std::cout<<"Running PMT find alpha..."<<std::endl;
            runFindAlphaFit(config);

        } else { 
            throw std::invalid_argument("Unknown reconstruction type '"+mode+"'");}

//     }
    
    return 0;
}
