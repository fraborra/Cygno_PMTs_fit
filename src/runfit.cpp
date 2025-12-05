//
//  runfit.cpp
//  LIMEPMTfits
//
//  Created by Stefano Piacentini on 23/09/22.
//  Modified by Francesco Borra on 28/06/23.
//

#include <BAT/BCLog.h>

#include <iostream>
#include <omp.h>
#include <stdexcept>
#include <string>
// #include <nlohmann/json.hpp>

#include "PMT_FindAlpha.hpp"
#include "PMT_association.hpp"
#include "PMT_calibration.hpp"
#include "PMT_singleEvent.hpp"
#include "helper_lib.hpp"

int main(int argc, char *argv[]) {

  if (argc < 2) {
    std::cerr << "Usage:\n";
    std::cerr << "  " << argv[0] << " <config_file>\n";
    return 1;
  }

  std::string configFile = argv[1];
  Config config = readConfigFile(configFile);
  std::string mode = config.mode;

  // if for mode
  if (mode == "association") {
    std::cout << "Running PMT association reconstruction..." << std::endl;
    runAssociationFit(config);

  } else if (mode == "PMTcalibration") {
    std::cout << "Running PMT calibration..." << std::endl;
    runCalibrationFit(config);

  } else if (mode == "PMTfindalpha") {
    std::cout << "Running PMT find alpha..." << std::endl;
    runFindAlphaFit(config);

  } else {
    throw std::invalid_argument("Unknown reconstruction type '" + mode + "'");
  }

  return 0;
}
