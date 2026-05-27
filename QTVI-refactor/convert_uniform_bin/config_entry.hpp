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

    double      ecgRate = 0.0;
    double      ppgRate = 0.0;
    double      finalSamplingRate = 0.0;

    // Native sampling rates (Hz) for the remaining channels. Zero means
    // "absent" -- the corresponding channel slot in the .bin gets a
    // missing-channel placeholder.
    // ART and ART_PULM share abpRate since they're variants on the
    // same arterial-pressure sensor.
    double      cvpRate = 0.0;
    double      abpRate = 0.0;
    double      respRate = 0.0;

    std::string input_path;
    std::string output_path;


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