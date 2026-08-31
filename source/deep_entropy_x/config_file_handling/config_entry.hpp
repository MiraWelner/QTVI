#pragma once
/**
 * @file   config_entry.hpp
 * @brief  a struct containing all the rates, paths, and labels
 *
 *         Every signal channel carries an explicit (raw_rate,
 *         upsample_rate) pair straight from config.csv. There is no
 *         single global upsample rate -- each channel is upsampled to its
 *         own configured target rate. Channels with no configured rate
 *         (or no label) are written as missing-channel placeholders.
 */

#include <string>

struct config_entry {
    //params set by the config.csv
    std::string dataset_type;
    std::string main_file_extention;
    std::string sleep_file_extention;
    std::string input_path;
    std::string output_path;

    // Per-channel native (raw) rate + target (upsample) rate, in Hz.
    double ecg_raw_rate = 0.0, ecg_upsample_rate = 0.0;
    double ppg_raw_rate = 0.0, ppg_upsample_rate = 0.0;
    double cvp_raw_rate = 0.0, cvp_upsample_rate = 0.0;
    double pres_raw_rate = 0.0, pres_upsample_rate = 0.0;
    double abp_raw_rate = 0.0, abp_upsample_rate = 0.0;
    double art_raw_rate = 0.0, art_upsample_rate = 0.0;
    double art_pulm_raw_rate = 0.0, art_pulm_upsample_rate = 0.0;
    double accel_raw_rate = 0.0, accel_upsample_rate = 0.0;
    double temp_raw_rate = 0.0, temp_upsample_rate = 0.0;
    double marker_raw_rate = 0.0, marker_upsample_rate = 0.0;
    double resp_raw_rate = 0.0, resp_upsample_rate = 0.0;
    double pacemaker_raw_rate = 0.0, pacemaker_upsample_rate = 0.0;
    double eeg_raw_rate = 0.0, eeg_upsample_rate = 0.0;
    double eeg_4_raw_rate = 0.0, eeg_4_upsample_rate = 0.0;
    double flow_raw_rate = 0.0, flow_upsample_rate = 0.0;
    double snore_raw_rate = 0.0, snore_upsample_rate = 0.0;
    double thor_raw_rate = 0.0, thor_upsample_rate = 0.0;
    double abdo_raw_rate = 0.0, abdo_upsample_rate = 0.0;
    double leg_raw_rate = 0.0, leg_upsample_rate = 0.0;
    double auxac_raw_rate = 0.0, auxac_upsample_rate = 0.0;
    double therm_raw_rate = 0.0, therm_upsample_rate = 0.0;
    double pos_raw_rate = 0.0, pos_upsample_rate = 0.0;
    double oxstatus_raw_rate = 0.0, oxstatus_upsample_rate = 0.0;
    double spo2_raw_rate = 0.0, spo2_upsample_rate = 0.0;
    double hr_raw_rate = 0.0, hr_upsample_rate = 0.0;
    double dhr_raw_rate = 0.0, dhr_upsample_rate = 0.0;
    double eog_l_raw_rate = 0.0, eog_l_upsample_rate = 0.0;
    double eog_r_raw_rate = 0.0, eog_r_upsample_rate = 0.0;
    double emg_raw_rate = 0.0, emg_upsample_rate = 0.0;

    double sleepstate_length = 0.0;
    double blanking_period = 0.0;
    double threshold = 0.0;
    double bin_size_minutes = 0.0;

    // Output subpaths used by the marking / viewer pipeline. output_path is
    // the user-set parent; the rest are derived from it by deriveSubpaths()
    // in config_loader. Ignored by the bin maker.
    std::string bin_file_path;
    std::string noise_data_path;
    std::string annealed_data_path;
    std::string r_peak_data_path;
    std::string template_path;
    std::string qtvi_marker_path;
	std::string quality_metric;
    std::string snapshot_path;
    std::string log_path;
    std::string training_log;
    std::string vcg_output;
    std::string bin_archive_path;

    /*
    Different filetypes have different terms for the same type of signal,
    If only one filetype has a given type of data (ie only bittium has accelration) then the label
    name is set here. Otherwise, it is set in the apply_dataset_specific_channel_labels function in config_loader.cpp
    */
    std::string ecg_1_label;
    std::string ecg_2_label;
    std::string ecg_3_label;
    std::string ppg_label;
    std::string eeg_1_label;
    std::string eeg_2_label;
    std::string eeg_3_label;
    std::string accel_x_label = "Accelerometer_X";
    std::string accel_y_label = "Accelerometer_Y";
    std::string accel_z_label = "Accelerometer_Z";
    std::string cvp_label = "NLS_NOM_PRESS_BLD_VEN_CENT";
    std::string resp_label = "NLS_NOM_RESP";
    std::string temp_label = "DEV_Temperature";
    std::string marker_label = "Marker";
    std::string pacemaker_label = "Pacemaker_events";
    std::string eeg_4_label = "NLS_EEG_NAMES_EEG_CHAN4";
    std::string abp_label = "NLS_NOM_PRESS_BLD_ART_ABP";
    std::string art_label = "NLS_NOM_PRESS_BLD_ART";
    std::string art_pulm_label = "NLS_NOM_PRESS_BLD_ART_PULM";
    std::string eog_l_label = "EOG-L";
    std::string eog_r_label = "EOG-R";
    std::string emg_label = "EMG";
    std::string pres_label = "Pres";
    std::string flow_label = "Flow";
    std::string snore_label = "Snore";
    std::string thor_label = "Thor";
    std::string abdo_label = "Abdo";
    std::string leg_label = "Leg";
    std::string auxac_label = "Aux_AC";
    std::string therm_label = "Therm";
    std::string pos_label = "Pos";
    std::string oxstatus_label = "OxStatus";
    std::string spo2_label = "SpO2";
    std::string hr_label = "HR";
    std::string dhr_label = "DHR";


    //determine what r peak finding method to utilize
    bool use_consensus_rpeak = true;
    // --- Filtering options ---
    // Powerline notch filter. 0 = disabled; valid enabled values are 50 or 60 Hz.
    int notch_filter_hz = 0;
    // Waveform high-pass cutoff in Hz. 0 = disabled; default 0.5 when enabled.
    double waveform_highpass_hz = 0.5;

    // --- Subject demographics (stored only; no downstream use yet) ---
    int    age = 0;
    std::string sex;               // stored verbatim, not interpreted
    double weight_kg = 0.0;
    double height_cm = 0.0;
    int    hr_rest = 0;
    int    hr_max = 0;
};