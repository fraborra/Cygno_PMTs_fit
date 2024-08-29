#include "helper_lib.hpp"
#include <sstream>

void how_to_use() {
    std::cout << "Usage: program -m mode -i input_file -s start_ind -e end_ind -o output_file [-p] [-c] [-l] [-h]" << std::endl;
}

// Config file reading
bool stringToBool(const std::string& str) {
    return str == "true" || str == "1";
}

// function to read config file
Config readConfigFile(const std::string& filename) {
    Config config;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Impossible to open config file: " << filename << std::endl;
        return config;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string key, value;

        if (std::getline(ss, key, '=') && std::getline(ss, value)) {
            if (key == "mode") config.mode = value;
            else if (key == "input_file") config.input_file = value;
            else if (key == "start_ind") config.start_ind = std::stoi(value);
            else if (key == "end_ind") config.end_ind = std::stoi(value);
            else if (key == "nPoints") config.nPoints = std::stoi(value);
            else if (key == "output_file") config.output_file = value;
            else if (key == "plot") config.plot = stringToBool(value);
            else if (key == "write_chains") config.write_chains = stringToBool(value);
            else if (key == "write_log") config.write_log = stringToBool(value);
            else if (key == "print_summary") config.print_summary = stringToBool(value);
        }
    }

    file.close();
    return config;
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
        } else if (mode.compare("PMTfindalpha") == 0){
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
        }
	else {
            throw std::runtime_error("No matched mode for file readout\n");
        }
    } else {
        throw std::runtime_error("Could not open the file\n");
    }
    file.close();

    return;
}
