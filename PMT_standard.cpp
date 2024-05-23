//
//  PMT_standard.cpp
//  LIMEPMTfits
//
//  Created by Stefano Piacentini on 23/09/22.
//  Modified by Francesco Borra on 28/06/23.
//

#include "PMT_standard.hpp"
#include <TMath.h>
#include <BAT/BCMath.h>
#include <cmath>

PMTfit::PMTfit(const std::string& name,
               std::string r_type, std::string datafile,
               int nth, int Np, int Iprec
//               std::string res_dir
               ) : BCModel(name)
{
    std::cout<<"Starting fit for '"<<r_type<<" reconstruction'"<<std::endl;

    r_type_ = r_type;
    Lmax = 200000;
    
    int cntrl = ReadInputFile(datafile);

    if (cntrl!=0){
      throw std::runtime_error("Couldn't read the input file");
    }

    Npoints  = Np;

    iprec_ = Iprec;
    
    //DEFINING parameters
    if (r_type_.compare("peaks") == 0) {
        AddParameter("L", 0, Lmax, "L", "[a.u.]");
//         GetParameter("L").Fix(1.0);

        for(unsigned int i=(iprec_+1)*Npoints; i<(iprec_+2)*Npoints; i++) {
            AddParameter("x_"+std::to_string(i), 0, 33, "x_"+std::to_string(i), "[cm]");
            AddParameter("y_"+std::to_string(i), 0, 33, "y_"+std::to_string(i), "[cm]");
        }
        
        AddParameter("c1", 0., 2., "c1", "[counts]");
        AddParameter("c2", 0., 2., "c2", "[counts]");
        AddParameter("c3", 0., 2., "c3", "[counts]");
        AddParameter("c4", 0., 2., "c4", "[counts]");
        
        //  FIXING CALIBRATION
        
        // 'Good' calibration
        GetParameter("c1").Fix(1.0);
        GetParameter("c2").Fix(0.965);
        GetParameter("c3").Fix(0.860);
        GetParameter("c4").Fix(0.827);

 	    // 'Bad' calibration
//         GetParameter("c1").Fix(1.0);
//         GetParameter("c2").Fix(1.22);                                                         
//         GetParameter("c3").Fix(0.529);                                                          
//         GetParameter("c4").Fix(0.672);

    } else if (r_type_.compare("slices") == 0){
        AddParameter("L", 0, Lmax, "L", "[a.u.]");

        for(unsigned int i=(iprec_+1)*Npoints; i<(iprec_+2)*Npoints; i++) {
            AddParameter("x_"+std::to_string(i), 0, 33, "x_"+std::to_string(i), "[cm]");
            AddParameter("y_"+std::to_string(i), 0, 33, "y_"+std::to_string(i), "[cm]");
        }
        
        AddParameter("c1", 0., 2., "c1", "[counts]");
        AddParameter("c2", 0., 2., "c2", "[counts]");
        AddParameter("c3", 0., 2., "c3", "[counts]");
        AddParameter("c4", 0., 2., "c4", "[counts]");

        // 'Good' calibration
        GetParameter("c1").Fix(1.0);
        GetParameter("c2").Fix(0.965);
        GetParameter("c3").Fix(0.860);
        GetParameter("c4").Fix(0.827);
        
    } else if (r_type_.compare("PMTcalibration") == 0){
        AddParameter("L", 0, Lmax, "L", "[a.u.]");
        GetParameter("L").Fix(40000.0);
        
        // Alternativa
//         GetParameter("L").Fix(1.0);

        for(unsigned int i=(iprec_+1)*Npoints; i<(iprec_+2)*Npoints; i++) {
            AddParameter("x_"+std::to_string(i), 0, 33, "x_"+std::to_string(i), "[cm]");
            AddParameter("y_"+std::to_string(i), 0, 33, "y_"+std::to_string(i), "[cm]");
        }
        
        AddParameter("c1", 0., 2., "c1", "[counts]");
        AddParameter("c2", 0., 2., "c2", "[counts]");
        AddParameter("c3", 0., 2., "c3", "[counts]");
        AddParameter("c4", 0., 2., "c4", "[counts]");
        
        // Alternativa
       
//        AddParameter("c1", 1000., 80000., "c1", "[counts]");
//        AddParameter("c2", 1000., 80000., "c2", "[counts]");
//        AddParameter("c3", 1000., 80000., "c3", "[counts]");
//        AddParameter("c4", 1000., 80000., "c4", "[counts]");

    } else if (r_type_.compare("POScomparison") == 0){
        throw std::runtime_error("Reconstruction '"+r_type_+"' not implemented yet.\n");
    } else {
        throw std::runtime_error("Unknown model '"+r_type_+"'.\n");
    }
    
    omp_set_dynamic(0);
    omp_set_num_threads(nth);
}



double PMTfit::LogLikelihood(const std::vector<double>& pars) {

    double LL = 0.;
    
//    if(r_type_.compare("standard")==0) {

        for(unsigned int i=0; i<Npoints; i++) {
            for(unsigned int j=0; j<4; j++) {
                double Lij = data[j][i+(iprec_+1)*Npoints]; // here the data
                double sLij = 0.1*Lij;                      // for now set to 10% of the integral
                
                int k = 2*Npoints +1 +j;
                
                double tmp = sqrt(D2(pars[2*i+1], pars[2*i+2], j));
                LL += BCMath::LogGaus(Lij,                             // x, namely Lij
                                      (pars[0]*pars[k])/(pow(tmp, 4)), // mu, namely the light computed in the step
                                      sLij,                            // sigma
                                      true                             // norm factor
                                      );
            }
        }
//     } else if (r_type_.compare("differentL") == 0){
//         throw std::runtime_error("Model '"+r_type_+"' not implemented yet.\n");
    // } else {
    //     throw std::runtime_error("Unknown model '"+r_type_+"'.\n");
    // }

    return LL;
}




double PMTfit::LogAPrioriProbability(const std::vector<double>& pars) {
    double LL = 0.;
        //flat priors everywhere
        if(pars[0]<0 || pars[0]>Lmax) {
            LL += log(0.0);
        } else {
            LL += log(1.0/Lmax);
        }
    
    for(unsigned int i=1; i<2*Npoints+1; i++) {
        if(pars[i]<0 || pars[i]>33.) {
            LL += log(0.0);
        } else {
            LL += log(1.0/33.);
        }
    }
    
    for(unsigned int j=(pars.size()-4); j<pars.size(); j++){
            if(pars[j]<0. || pars[j]>2.) {
                LL += log(0.0);
            } else {
                LL += log(1.0/2.);
            }
        }
    // Alternativa
//     for(unsigned int j=(pars.size()-4); j<pars.size(); j++){
//         if(pars[j]<1000. || pars[j]>80000.) {
//             LL += log(0.0);
//         } else {
//             LL += log(1.0/79000.);
//         }
//     }

    return LL;
}


double PMTfit::D2(double x, double y, int i) {
    if (i == 0) {
        return (x - x1)*(x - x1) + (y - y1)*(y - y1) + zGEM*zGEM;
    } else if (i == 1) {
        return (x - x2)*(x - x2) + (y - y2)*(y - y2) + zGEM*zGEM;
    } else if (i == 2) {
        return (x - x3)*(x - x3) + (y - y3)*(y - y3) + zGEM*zGEM;
    } else if (i == 3) {
        return (x - x4)*(x - x4) + (y - y4)*(y - y4) + zGEM*zGEM;
    } else {
        throw std::runtime_error("Uknown value of PMT index.\n");
    }
    return 0.;
}




int PMTfit::ReadInputFile(std::string filename) {

    std::ifstream myfile;
    myfile.open(filename);

    std::string runstr, evstr, trgstr, L1str, L2str, L3str, L4str;
    std::string GEstr, indxstr, mj2str;                   //Analisi mia
    std::string xx, yy, sx, integ, sL1, sL2, sL3, sL4;    //Analisi di Matteo
    std::string p1str, p2str, p3str, p4str, w1str, w2str, w3str, w4str;
    
    if(myfile.is_open()) {
        if(r_type_.compare("peaks") == 0) {
            while(myfile >>                              //Analisi mia
                  runstr >> evstr >> trgstr >> GEstr >> indxstr >>
                  L1str >> L2str >> L3str >> L4str >> mj2str)
            {
                //       while(myfile >>                                 //Analisi di Matteo
                //             runstr >> evstr >> trgstr >>
                //             xx >> yy >> sx >> integ >>
                //             L1str >> L2str >> L3str >> L4str>>
                //             sL1 >> sL2 >> sL3 >> sL4
                //             ) {
                
                run.push_back(stoi(runstr));
                event.push_back(stoi(evstr));
                trigger.push_back(stoi(trgstr));
                
                GE.push_back(stoi(GEstr));                //Analisi mia
                indx.push_back(stoi(indxstr));
                mj2.push_back(stoi(mj2str));
                //           time.push_back(stod(timestr));
                
                //          xtrue.push_back(stod(xx));                     //Analisi di Matteo
                //          ytrue.push_back(stod(yy));
                //          sX.push_back(stod(sx));
                //          integral.push_back(stod(integ));
                
                double L1 = stod(L1str);
                data[0].push_back(L1);
                
                double L2 = stod(L2str);
                data[1].push_back(L2);
                
                double L3 = stod(L3str);
                data[2].push_back(L3);
                
                double L4 = stod(L4str);
                data[3].push_back(L4);
                
            }
        } else if (r_type_.compare("slices") == 0){
            while(myfile >>
                  runstr >> evstr >> trgstr >> GEstr >> indxstr >>
                  L1str >> L2str >> L3str >> L4str >> mj2str
                  ) {
                
                run.push_back(stoi(runstr));
                event.push_back(stoi(evstr));
                trigger.push_back(stoi(trgstr));
                GE.push_back(stoi(GEstr));
                indx.push_back(stoi(indxstr));
                
                double L1 = stod(L1str);
                data[0].push_back(L1);
                
                double L2 = stod(L2str);
                data[1].push_back(L2);
                
                double L3 = stod(L3str);
                data[2].push_back(L3);
                
                double L4 = stod(L4str);
                data[3].push_back(L4);
            }
            
        } else if (r_type_.compare("PMTcalibration") == 0){
//             int count = 0;
            while(myfile >>
                  runstr >> evstr >> trgstr >> GEstr >> indxstr >>
                  L1str >> L2str >> L3str >> L4str >> mj2str>>
                  p1str >> p2str >> p3str >> p4str >>
                  w1str >> w2str >> w3str >> w4str
                  ) {
//                 if (count<10){
//                 std::cout<<runstr<< "\t" <<evstr<< "\t" <<trgstr<< "\t" <<GEstr<< "\t" <<indxstr
//                     << "\t" <<L1str<< "\t" <<L2str<< "\t" <<L3str<< "\t" <<L4str<<std::endl;
//                 count +=1;
//                 }
                
                    
                run.push_back(stoi(runstr));
                event.push_back(stoi(evstr));
                trigger.push_back(stoi(trgstr));
                GE.push_back(stoi(GEstr));
                indx.push_back(stoi(indxstr));
                
                double L1 = stod(L1str);
                data[0].push_back(L1);
                
                double L2 = stod(L2str);
                data[1].push_back(L2);
                
                double L3 = stod(L3str);
                data[2].push_back(L3);
                
                double L4 = stod(L4str);
                data[3].push_back(L4);
            }
        } else if (r_type_.compare("POScomparison") == 0){
            throw std::runtime_error("Unknown model '"+r_type_+"'.\n");
        }
    } else {
        throw std::runtime_error("Could not open the file\n");
    }
    myfile.close();

    return 0;
}

// "returners"
int PMTfit::getRun(int i){
    return run[i];
}
int PMTfit::getEvent(int i){
    return event[i];
}
int PMTfit::getTrigger(int i){
    return trigger[i];
}
int PMTfit::getMj2(int i){
    return mj2[i];
}
int PMTfit::getGE(int i){
    return GE[i];
}
int PMTfit::getIndx(int i){
    return indx[i];
}
double PMTfit::getXtrue(int i){
    return xtrue[i];
}
double PMTfit::getYtrue(int i){
    return ytrue[i];
}
double PMTfit::getSigX(int i){
    return sX[i];
}
double PMTfit::getIntegral(int i){
    return integral[i];
}
