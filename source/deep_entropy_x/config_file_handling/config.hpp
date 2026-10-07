#pragma once
/**
 * @file   config.hpp
 * @brief  config_entry is a struct containing all the rates, paths, lables, etc in the config.csv file.
 *
 *         SignalRates lives here too, because it is nothing but a projection of
 *         config_entry: five upsample rates and the two morphology half-window
 *         fractions, and no more than that. It was in template_structs.hpp,
 *         which is a header of template-generation intermediates and has
 *         nothing to do with either -- the generation stage was simply the only
 *         consumer, so the type ended up filed under its consumer rather than
 *         its source. config_entry::signalRates() is the one place it is built.
 */

#include <string>

 // Per-channel sample rates (Hz), threaded through the template-generation
 // pipeline. Each channel is at its own rate; the slicer converts between them
 // via the ratio channelRate / ecgRate. A rate of 0 means the channel is absent
 // from this dataset (skip it).
 //
 // A PROJECTION, NOT A SUBSET OF THE CONFIG. The generation stage takes this
 // rather than a config_entry on purpose: seven doubles it needs instead of a
 // hundred fields it does not, so nothing down there can reach a path, a label
 // or a filter cutoff.
struct SignalRates {
    double ecg = 0.0;
    double ppg = 0.0;
    double abp = 0.0;
    double art = 0.0;
    double artPulm = 0.0;
    double morph_halfwin_ecg_pct_rr = 0.0;
    double morph_halfwin_ppg_pct_rr = 0.0;
    // PULSE WINDOW OFFSETS, ms: the measured, confident ECG-to-channel
    // hardware lag (channel_offset), or 0. Where each beat's pulse landmarks
    // are SEARCHED for -- [R + lag, next R + lag] -- and nothing else; the
    // signals are never shifted. Set by analysis_job::prepare after the lag is
    // measured; 0 for every non-CHAOS record.
    double ppg_lag_ms = 0.0;
    double arterial_lag_ms = 0.0;
};

struct config_entry {
    //params set by the config.csv
    std::string dataset_type;
    std::string main_file_extension;
    std::string sleep_file_extension;
    std::string input_path;
    std::string output_path;
    std::string cgm_folder;

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

	//thresholds and other parameters
    double sleepstate_length = 0.0;
    double blanking_period = 0.0;
    double threshold = 0.0;
    double bin_size_minutes = 0.0;
    double cgm_bin_minutes = 60.0;
    double morph_threshold_ecg = 0.0;
    double morph_threshold_ppg = 0.0;
    double ppg_pearson_threshold = 0.0;
    double ecg_tukey_fence = 1.5;
    double ppg_tukey_fence = 1.5;
    int min_beats_template_ecg = 0;
    int min_beats_template_ppg = 0;
    double region_around_Rpeak_for_morphology_split = 0.0;
    double region_around_PPGPeak_for_morphology_split = 0.0;

    // Output subpaths used by the marking / viewer pipeline. output_path is
    // the user-set parent; the rest are derived from it by deriveSubpaths()
    // in config_loader. Ignored by the bin maker.
    std::string noise_data_path;
    std::string annealed_data_path;
    std::string r_peak_data_path;
    std::string template_path;
    std::string fiducial_marker_locations;
    std::string snapshot_path;
    std::string log_path;
    std::string training_log;
    std::string vcg_output;
    std::string cgm_output_path;

    // Different filetypes have different terms for the same type of signal,
    // If only one filetype has a given type of data (ie only bittium has accelration) then the label
    // name is set here. Otherwise, it is set in the apply_dataset_specific_channel_labels function in config_loader.cpp
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


    bool use_consensus_rpeak = true;
    bool exclude_beat_after_ectopic = true;
    //if you are trying to reload somebody's old markings and they don't work bc the splits are different - set this to true (in config.csv not here)
    bool override_morphology = false;
    double notch_filter_hz = 0.0;
    double waveform_highpass_hz = 0.0;

    // --- Subject demographics (stored only; no downstream use yet) ---
    int    age = 0;
    std::string sex;
    double weight_kg = 0.0;
    double height_cm = 0.0;
    double hr_rest = 0.0;
    double hr_max = 0.0;


    SignalRates signalRates() const {
        //the struct that carries all signal rates, as well as the regions around the peaks for which
        //morphologies are split based on
        SignalRates r;
        r.ecg = ecg_upsample_rate;
        r.ppg = ppg_upsample_rate;
        r.abp = abp_upsample_rate;
        r.art = art_upsample_rate;
        r.artPulm = art_pulm_upsample_rate;
        r.morph_halfwin_ecg_pct_rr = region_around_Rpeak_for_morphology_split;
        r.morph_halfwin_ppg_pct_rr = region_around_PPGPeak_for_morphology_split;
        return r;
    }
};