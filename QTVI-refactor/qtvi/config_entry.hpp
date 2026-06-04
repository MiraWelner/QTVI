#pragma once
/**
 * @file   config_entry.hpp
 * @brief  a struct containing all the rates, paths, and labels
 */

#include <string>

struct config_entry {
    std::string dataType;             // "MESA" / "Bittium" / "CHAOS"
    std::string mainExt;              // ".edf" / ".dat"
    std::string sleepExt;             // ".XML" or empty

    double      ecgRate = 0.0;          // ECG1/2/3 native rate, also = .dat row rate for CHAOS
    double      ppgRate = 0.0;
    double      finalSamplingRate = 1000.0;

    // Native sampling rates (Hz) for the remaining channels written by
    // file_to_bin. Zero means "absent" (channel emits a placeholder).
    // ART and ART_PULM share abpRate since they're variants on the same
    // arterial-pressure sensor.
    double      cvpRate = 0.0;          // central venous pressure (CH_PRES)
    double      abpRate = 0.0;          // also used for ART and ART_PULM
    double      respRate = 0.0;


    //the output_path is never directly used, but is rather the
    //parent folder of all the below paths. It is its own entry in the config
    // becuase it is can be set by the user, the below folders are derived from it
    std::string input_path;          // raw source folder (.dat / .edf input);
    // read by file_to_bin, unused by the GUI
    std::string output_path;
    std::string bin_file_path;
    std::string noise_data_path;
    std::string annealed_data_path;
    std::string r_peak_data_path;
    std::string template_path;
    std::string qtvi_marker_path;
    std::string snapshot_path;



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


    double bin_length_minutes;
};