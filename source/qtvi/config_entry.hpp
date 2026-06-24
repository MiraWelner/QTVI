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
    std::string ecg_1_label;
    std::string ecg_2_label;
    std::string ecg_3_label;
    std::string ppg_label;
    std::string accel_x_label = "Accelerometer_X";
    std::string accel_y_label = "Accelerometer_Y";
    std::string accel_z_label = "Accelerometer_Z";
    std::string cvp_label = "NLS_NOM_PRESS_BLD_VEN_CENT";
    std::string resp_label = "NLS_NOM_RESP";
    std::string temp_label = "DEV_Temperature";
    std::string marker_label = "Marker";
    std::string pacemaker_label = "Pacemaker_events";
    std::string eeg_1_label = "NLS_EEG_NAMES_EEG_CHAN1";
    std::string eeg_2_label = "NLS_EEG_NAMES_EEG_CHAN2";
    std::string eeg_3_label = "NLS_EEG_NAMES_EEG_CHAN3";
    std::string eeg_4_label = "NLS_EEG_NAMES_EEG_CHAN4";
    std::string abp_label = "NLS_NOM_PRESS_BLD_ART_ABP";
    std::string art_label = "NLS_NOM_PRESS_BLD_ART";
    std::string art_pulm_label = "NLS_NOM_PRESS_BLD_ART_PULM";


};