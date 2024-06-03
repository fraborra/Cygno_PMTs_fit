#include "helper_lib.hpp"
void how_to_use() {
    std::cout << "Usage: program -m mode -i input_file -s start_ind -e end_ind -o output_file [-p] [-c] [-l] [-h]" << std::endl;
}

// DataReader class

// "returners"
const std::vector<int>& DataReader::getRun(){
    return run;
}
const std::vector<int>& DataReader::getEvent(){
    return event;
}
const std::vector<int>& DataReader::getTrigger(){
    return trigger;
}
const std::vector<int>& DataReader::getIndx(){
    return indx;
}
const std::vector<double>& DataReader::getXtrue(){
    return xtrue;
}
const std::vector<double>& DataReader::getYtrue(){
    return ytrue;
}
const std::vector<double>& DataReader::getL1(){
    return L1;
}
const std::vector<double>& DataReader::getL2(){
    return L2;
}
const std::vector<double>& DataReader::getL3(){
    return L3;
}
const std::vector<double>& DataReader::getL4(){
    return L4;
}

//File reader
void DataReader::readFile(const std::string& input_file, const std::string& mode){
    std::ifstream file(input_file);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    std::string runstr, evstr, trgstr, indxstr, L1str, L2str, L3str, L4str;
    std::string xtruestr, ytruestr;                   

    if(file.is_open()) {
        if(mode.compare("association") == 0) {
            while(file >>                              //Analisi mia
                  runstr >> evstr >> trgstr >> indxstr >>
                  L1str >> L2str >> L3str >> L4str)
            {                
                run.push_back(stoi(runstr));
                event.push_back(stoi(evstr));
                trigger.push_back(stoi(trgstr));
                indx.push_back(stoi(indxstr));
                
                L1.push_back(stod(L1str));
                L2.push_back(stod(L2str));
                L3.push_back(stod(L3str));
                L4.push_back(stod(L4str));
            }
        } else if (mode.compare("PMTcalibration") == 0){
            while(file >>
                  runstr >> evstr >> trgstr >> indxstr >>
                  L1str >> L2str >> L3str >> L4str >> xtruestr >> ytruestr
                  ) {                
                run.push_back(stoi(runstr));
                event.push_back(stoi(evstr));
                trigger.push_back(stoi(trgstr));
                indx.push_back(stoi(indxstr));
                
                L1.push_back(stod(L1str));
                L2.push_back(stod(L2str));
                L3.push_back(stod(L3str));
                L4.push_back(stod(L4str));
                xtrue.push_back(stod(xtruestr));
                ytrue.push_back(stod(ytruestr));
            }
        } else {
            throw std::runtime_error("No matched mode for file readout\n");
        }
    } else {
        throw std::runtime_error("Could not open the file\n");
    }
    file.close();

    return;
}