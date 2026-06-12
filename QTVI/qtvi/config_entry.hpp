#pragma once
/**
 * @file   config_entry.hpp
 * @brief  a struct containing all the rates, paths, and labels
 */

#include <string>

struct config_entry {
    //these params are always set by the config.csv, if they are not set by the config.csv, then the file will fail to load
    std::string dataset_type;
    std::string original_file_extention;
    std::string sleep_file_extention;
    double      ecg_rate;
    double      ppg_rate;
    double      central_venous_pressure_rate;
    double      arterial_blood_pressure_rate;
    double      resp_rate;
    double      target_sampling_rate; //this rate is not defined by the input data, rather it is the rate that is upsampled to
    double      bin_length_minutes;
    double      blanking_period; // this is the percent of time after the r peak in which a new r peak cannot occur
    double      height_threshold_percent; //threshold for r peaks - but not raw threshold, the percent of the height distribution. So 1.0 means no R peaks can ever be found, 0.5 means the threshold is the median



    //the output_path is never directly used, but is rather the
    //parent folder of all the below paths. It is its own entry in the config
    // becuase it is can be set by the user, the below folders are derived from it
    std::string input_path;          // raw source folder (.dat / .edf input);
    std::string output_path;
    std::string bin_file_path;
    std::string noise_data_path;
    std::string annealed_data_path;
    std::string r_peak_data_path;
    std::string template_path;
    std::string qtvi_marker_path;
    std::string snapshot_path;
    std::string log_path;


    // these parameters are not in the config but rather in input_file_handler.cpp. 
    // different filetypes have different labels for the same type of signal - for example, ECG_1 vs EKG is talking about the same thing
    std::string ecg1Label;
    std::string ecg2Label;
    std::string ecg3Label;
    std::string ppgLabel;
    std::string accelXLabel;
    std::string accelYLabel;
    std::string accelZLabel;
    std::string cvpLabel;
    std::string respLabel;
    std::string eeg1Label;
    std::string eeg2Label;
    std::string eeg3Label;
    std::string eeg4Label;
    std::string abpLabel;
    std::string artLabel;
    std::string artPulmLabel;


};