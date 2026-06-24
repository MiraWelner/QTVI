#pragma once
/**
 * @file   config_entry.hpp
 * @brief  a struct containing all the rates, paths, and labels
 */

#include <string>

struct config_entry {
    //params set by the config.csv
    std::string dataset_type;
    std::string original_file_extention;
    std::string sleep_file_extention;
    std::string input_path;
    std::string output_path;

    double      ecg_rate;
    double      ppg_rate;
    double      central_venous_pressure_rate;
    double      resp_rate;
    double      pacemaker_event_rate;
    double      arterial_blood_pressure_rate;
    double      accel_rate;
    double      temp_rate;
    double      marker_rate;
    double      sleep_state_length;
    double      low_sample_rate;//the upsample rate of marker, temperature, and other things that are originally 1hz
    double      high_upsample_rate;//the upsample rate of ecg,ppg,etc

    /*
    Different filetypes have different terms for the same type of signal,
    If only one filetype has a given type of data (ie only bittium has accelration) then the label
    name is set here. Otherwise, it is set in the apply_dataset_specific_channel_labels function in config_loader.cpp
    */
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