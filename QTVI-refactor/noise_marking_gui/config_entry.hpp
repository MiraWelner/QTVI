
#pragma once
/**
 * @file   config_entry.hpp
 * @brief  a struct containing all the 
 */

#include <string>

struct config_entry {
    std::string dataType;             // "MESA" / "Bittium" / "CHAOS"
    std::string mainExt;              // ".edf" / ".dat"
    std::string sleepExt;             // ".XML" or empty

    double      ecgRate = 0.0;
    double      ppgRate = 0.0;
    double      finalSamplingRate = 1000.0;

    std::string original_file_path;
    std::string bin_file_path;       
    std::string noise_data_path;     
    std::string annealed_data_path;  
    std::string r_peak_data_path;
    std::string template_path;



    //this stuff is not in the config but is rather in input_file_handler.cpp. 
    // it is put here so it can easily be added into the config and because it is dependent on the type of file
    std::string ecg1Label;
    std::string ecg2Label;
    std::string ecg3Label;
    std::string ppgLabel;

    std::string accelXLabel;
    std::string accelYLabel;
    std::string accelZLabel;

    double bin_length_minutes;
};