//
//  PMT_singleEvent.hpp
//  Francesco Borra, Dic 2025
//

#ifndef PMT_singleEvent_hpp
#define PMT_singleEvent_hpp

#include <BAT/BCModel.h>

#include <fstream>
#include <iostream>
#include <omp.h>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

#include "Math/ProbFunc.h"
#include "PMT_association.hpp"
#include "helper_lib.hpp"

SinglePointResult runSingleEvent(const double q_values[4],
                                 const double calib[4], int NIterPrerun,
                                 int NIter, bool second_round);

#endif /* PMT_singleEvent_hpp */
