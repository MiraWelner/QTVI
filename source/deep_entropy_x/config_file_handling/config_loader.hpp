#pragma once
/**
 * @file   config_loader.hpp
 * @brief  Loads the config.csv file, parses it based on the datset type selected by the user, and fills up a config_entry struct with the relevant paths, rates, and channel labels.
 *         The channel labels (eg. "ECG_1" vs "EKG") are dataset-specific but not in the config file, so they are assigned in apply_dataset_specific_channel_labels() based on the dataset type.
 *         The output paths are either found in the config file, or prompted for manually if the config file cells are blank. The output_path is used to create the subfolders where the
 *         specific types of output are found.

 */

#include "config_file_handling/config.hpp"

#include <QFileDialog>

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>

namespace config_loader_detail {

    inline double stod_or_default(const std::string& s, double default_value) {
        //string to double, or return default if the string is not a number
        try { return std::stod(s); }
        catch (...) { return default_value; }
    }

    inline std::string normalize_key(std::string s) {
        for (char& ch : s)
            ch = (char)std::tolower((unsigned char)ch);
        return s;
    }

    inline std::vector<std::string> parse_config_row(const std::string& line) {
        // The config file is a csv - this just is a util for loading a .csv row
        std::vector<std::string> fields;
        std::string cur;
        for (char c : line) {
            if (c == ',') { fields.push_back(cur); cur.clear(); }
            else          cur += c;
        }
        fields.push_back(cur);
        for (auto& f : fields) {
            size_t first = f.find_first_not_of(" \t\r\n");
            if (first == std::string::npos) { f.clear(); continue; }
            size_t last = f.find_last_not_of(" \t\r\n");
            f = f.substr(first, last - first + 1);
        }
        return fields;
    }

    inline bool parseBool(const std::string& s, bool dflt) {
        //this is for the the r consensus variable in the config.csv file
        if (s.empty()) return dflt;
        std::string t; for (char c : s) t += std::tolower((unsigned char)c);
        if (t == "1" || t == "true" || t == "yes") return true;
        if (t == "0" || t == "false" || t == "no") return false;
        return dflt;
    }


    inline void apply_dataset_specific_channel_labels(config_entry& cfg) {
        /* The ECG and PPG channels exist in different datasets, but their names are
            different so they are set here. */
        if (cfg.dataset_type == "MESA") {
            cfg.ecg_1_label = "EKG";
            cfg.ppg_label = "Pleth";
            cfg.eeg_1_label = "EEG1";
            cfg.eeg_2_label = "EEG2";
            cfg.eeg_3_label = "EEG3";
        }
        else if (cfg.dataset_type == "BITTIUM") {
            cfg.ecg_1_label = "ECG_1";
            cfg.ecg_2_label = "ECG_2";
            cfg.ecg_3_label = "ECG_3";
        }

        else if (cfg.dataset_type == "CHAOS") {
            cfg.ecg_1_label = "NLS_NOM_ECG_ELEC_POTL_I";
            cfg.ecg_2_label = "NLS_NOM_ECG_ELEC_POTL_II";
            cfg.ecg_3_label = "NLS_NOM_ECG_ELEC_POTL_III";
            cfg.ppg_label = "NLS_NOM_PULS_OXIM_PLETH";
            cfg.eeg_1_label = "NLS_EEG_NAMES_EEG_CHAN1";
            cfg.eeg_2_label = "NLS_EEG_NAMES_EEG_CHAN2";
            cfg.eeg_3_label = "NLS_EEG_NAMES_EEG_CHAN3";
        }
        else if (cfg.dataset_type == "SHHS") {
            cfg.ecg_1_label = "ECG";
            cfg.eeg_1_label = "EEG";
            cfg.eeg_2_label = "EEG(sec)";
            cfg.eog_l_label = "EOG(L)";
            cfg.eog_r_label = "EOG(R)";
            cfg.emg_label = "EMG";
            cfg.flow_label = "AIRFLOW";
            cfg.thor_label = "THOR RES";
            cfg.abdo_label = "ABDO RES";
            cfg.pos_label = "POSITION";
            cfg.oxstatus_label = "OX STAT";
            cfg.spo2_label = "SaO2";
            cfg.hr_label = "PR";
            cfg.ecg_2_label.clear();
            cfg.ecg_3_label.clear();
            cfg.ppg_label.clear();
        }
    }


    inline void create_output_folders(config_entry& cfg) {
        //make the subfolders to organize the output data
        auto create_subfolder = [&cfg](const char* name) {
            const std::string p = cfg.output_path + "/" + name + "/";
            std::error_code ec;
            std::filesystem::create_directories(p, ec);
            return p;
            };

        cfg.annealed_data_path = create_subfolder("annealed_output");
        cfg.noise_data_path = create_subfolder("noise_marking_output");
        cfg.r_peak_data_path = create_subfolder("r_peak_finding_output");
        cfg.template_path = create_subfolder("template_outputs");
        cfg.fiducial_marker_locations = create_subfolder("fiducial_marker_locations");
        cfg.logs = create_subfolder("logs");
        cfg.snapshot_path = create_subfolder("snapshot_path");
        cfg.vcg_output = create_subfolder("vcg_output");
        // Only when there is CGM to process, so other datasets get no empty folder.
        if (cfg.dataset_type == "BITTIUM" && !cfg.cgm_folder.empty())
            cfg.cgm_output_path = create_subfolder("cgm_output");
    }
    inline bool prompt_for_missing_folders(config_entry& cfg) {
        // If the input or output folder is not in the config.csv (i.e. its field is empty), prompt the user to select it.
        const std::vector<std::pair<const char*, std::string*>> fields = {
            { "Bin Files:", &cfg.input_path },
            { "Output", &cfg.output_path },
        };

        for (const auto& [label, fieldPtr] : fields) {
            if (!fieldPtr->empty()) continue;

            QString title = QString("%1 (%2)").arg(label, QString::fromStdString(cfg.dataset_type));
            QString chosen = QFileDialog::getExistingDirectory(
                nullptr, title, QString(),
                QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
            if (chosen.isEmpty()) return false;
            *fieldPtr = chosen.toStdString();
        }
        return true;
    }
}

// Creates the output subfolders under cfg.output_path and fills in their paths.
// Call it once the input / output folders are final.
using config_loader_detail::create_output_folders;

inline bool load_config(int dataType, config_entry& out) {
    constexpr const char* CONFIG_PATH = "config.csv";
    using namespace config_loader_detail;
    std::ifstream file(CONFIG_PATH);
    if (!file.is_open()) {
        std::cerr << "ERROR: cannot open " << CONFIG_PATH << "\n";
        return false;
    }
    std::string user_selected_dataset = (dataType == 1) ? "MESA"
        : (dataType == 2) ? "BITTIUM"
        : (dataType == 3) ? "CHAOS"
        : (dataType == 4) ? "SHHS" : "";
    if (user_selected_dataset.empty()) {
        std::cerr << "ERROR: unknown dataset selection " << dataType << "\n";
        return false;
    }

    // ---- THE LAYOUT: ONE ROW PER SETTING, ONE COLUMN PER DATASET ------------
    //
    //       parameter,mesa,bittium,chaos,shhs
    //       ecg_raw_rate,256,500,500,125
    //       ...
    //
    // The header's first cell is "parameter"; each later header cell names a
    // dataset. The chosen dataset's column becomes one map, setting name ->
    // value, which get_value_from_config reads. Setting and dataset names are
    // matched case-insensitively, surrounding spaces ignored. A row whose
    // first cell starts with '#' is a comment; a blank row is skipped.
    std::string header;
    if (!std::getline(file, header)) {
        std::cerr << "ERROR: " << CONFIG_PATH << " is empty\n";
        return false;
    }
    std::vector<std::string> headerFields = parse_config_row(header);
    // Excel's "CSV UTF-8" save puts a byte-order mark in front of the first
    // cell, which would make "parameter" not match.
    if (!headerFields.empty() && headerFields[0].rfind("\xEF\xBB\xBF", 0) == 0)
        headerFields[0].erase(0, 3);
    const std::string firstCell = headerFields.empty() ? std::string() : normalize_key(headerFields[0]);
    auto upper = [](std::string v) {
        std::transform(v.begin(), v.end(), v.begin(), ::toupper);
        return v;
        };

    if (firstCell != "parameter") {
        std::cerr << "ERROR: " << CONFIG_PATH << " must start with a 'parameter' column "
            "(one row per setting, one column per dataset); its first cell is '"
            << (headerFields.empty() ? std::string() : headerFields[0]) << "'\n";
        return false;
    }
    int dcol = -1;
    for (int i = 1; i < (int)headerFields.size(); ++i)
        if (upper(headerFields[i]) == user_selected_dataset) { dcol = i; break; }
    if (dcol < 0) {
        std::cerr << "ERROR: no " << user_selected_dataset << " column in " << CONFIG_PATH << "\n";
        return false;
    }

    std::unordered_map<std::string, std::string> values;
    std::string line;
    while (std::getline(file, line)) {
        const std::vector<std::string> row = parse_config_row(line);
        if (row.empty() || row[0].empty() || row[0][0] == '#') continue;
        const std::string key = normalize_key(row[0]);
        if (values.count(key))
            std::cerr << "WARNING: " << CONFIG_PATH << " lists '" << row[0]
            << "' more than once; the last one is used\n";
        values[key] = (dcol < (int)row.size()) ? row[dcol] : std::string();
    }
    auto get_value_from_config = [&](const std::string& name) -> std::string {
        const auto it = values.find(normalize_key(name));
        return it == values.end() ? std::string() : it->second;
        };

    {

        out.dataset_type = user_selected_dataset;
        out.main_file_extension = get_value_from_config("main_file_extension");
        out.sleep_file_extension = get_value_from_config("sleep_file_extension");

        // NUMERIC SETTINGS: config.csv key, config_entry field, value when blank.
        // One row per setting, so a key and the field it fills sit side by side.
        struct NumericKey { const char* key; double config_entry::* field; double whenBlank; };
        static const NumericKey kNumeric[] = {
            { "ecg_raw_rate",                               &config_entry::ecg_raw_rate, 0.0 },
            { "ecg_upsample_rate",                          &config_entry::ecg_upsample_rate, 0.0 },
            { "ppg_raw_rate",                               &config_entry::ppg_raw_rate, 0.0 },
            { "ppg_upsample_rate",                          &config_entry::ppg_upsample_rate, 0.0 },
            { "cvp_raw_rate",                               &config_entry::cvp_raw_rate, 0.0 },
            { "cvp_upsample_rate",                          &config_entry::cvp_upsample_rate, 0.0 },
            { "pres_raw_rate",                              &config_entry::pres_raw_rate, 0.0 },
            { "pres_upsample_rate",                         &config_entry::pres_upsample_rate, 0.0 },
            { "abp_raw_rate",                               &config_entry::abp_raw_rate, 0.0 },
            { "abp_upsample_rate",                          &config_entry::abp_upsample_rate, 0.0 },
            { "art_raw_rate",                               &config_entry::art_raw_rate, 0.0 },
            { "art_upsample_rate",                          &config_entry::art_upsample_rate, 0.0 },
            { "art_pulm_raw_rate",                          &config_entry::art_pulm_raw_rate, 0.0 },
            { "art_pulm_upsample_rate",                     &config_entry::art_pulm_upsample_rate, 0.0 },
            { "accel_raw_rate",                             &config_entry::accel_raw_rate, 0.0 },
            { "accel_upsample_rate",                        &config_entry::accel_upsample_rate, 0.0 },
            { "accel_epoch_sec",                            &config_entry::accel_epoch_sec, 30.0 },
            { "accel_valid_epoch_pct",                      &config_entry::accel_valid_epoch_pct, 80.0 },
            { "temp_raw_rate",                              &config_entry::temp_raw_rate, 0.0 },
            { "temp_upsample_rate",                         &config_entry::temp_upsample_rate, 0.0 },
            { "marker_raw_rate",                            &config_entry::marker_raw_rate, 0.0 },
            { "marker_upsample_rate",                       &config_entry::marker_upsample_rate, 0.0 },
            { "resp_raw_rate",                              &config_entry::resp_raw_rate, 0.0 },
            { "resp_upsample_rate",                         &config_entry::resp_upsample_rate, 0.0 },
            { "pacemaker_event_raw_rate",                   &config_entry::pacemaker_raw_rate, 0.0 },
            { "pacemaker_event_upsample_rate",              &config_entry::pacemaker_upsample_rate, 0.0 },
            { "eeg_raw_rate",                               &config_entry::eeg_raw_rate, 0.0 },
            { "eeg_upsample_rate",                          &config_entry::eeg_upsample_rate, 0.0 },
            { "eogl_raw_rate",                              &config_entry::eog_l_raw_rate, 0.0 },
            { "eogl_upsample_rate",                         &config_entry::eog_l_upsample_rate, 0.0 },
            { "eogr_raw_rate",                              &config_entry::eog_r_raw_rate, 0.0 },
            { "eogr_upsample_rate",                         &config_entry::eog_r_upsample_rate, 0.0 },
            { "emg_raw_rate",                               &config_entry::emg_raw_rate, 0.0 },
            { "emg_upsample_rate",                          &config_entry::emg_upsample_rate, 0.0 },
            { "flow_raw_rate",                              &config_entry::flow_raw_rate, 0.0 },
            { "flow_upsample_rate",                         &config_entry::flow_upsample_rate, 0.0 },
            { "snore_raw_rate",                             &config_entry::snore_raw_rate, 0.0 },
            { "snore_upsample_rate",                        &config_entry::snore_upsample_rate, 0.0 },
            { "thor_raw_rate",                              &config_entry::thor_raw_rate, 0.0 },
            { "thor_upsample_rate",                         &config_entry::thor_upsample_rate, 0.0 },
            { "abdo_raw_rate",                              &config_entry::abdo_raw_rate, 0.0 },
            { "abdo_upsample_rate",                         &config_entry::abdo_upsample_rate, 0.0 },
            { "leg_raw_rate",                               &config_entry::leg_raw_rate, 0.0 },
            { "leg_upsample_rate",                          &config_entry::leg_upsample_rate, 0.0 },
            { "auxac_raw_rate",                             &config_entry::auxac_raw_rate, 0.0 },
            { "auxac_upsample_rate",                        &config_entry::auxac_upsample_rate, 0.0 },
            { "therm_raw_rate",                             &config_entry::therm_raw_rate, 0.0 },
            { "therm_upsample_rate",                        &config_entry::therm_upsample_rate, 0.0 },
            { "pos_raw_rate",                               &config_entry::pos_raw_rate, 0.0 },
            { "pos_upsample_rate",                          &config_entry::pos_upsample_rate, 0.0 },
            { "oxstatus_raw_rate",                          &config_entry::oxstatus_raw_rate, 0.0 },
            { "oxstatus_upsample_rate",                     &config_entry::oxstatus_upsample_rate, 0.0 },
            { "spo2_raw_rate",                              &config_entry::spo2_raw_rate, 0.0 },
            { "spo2_upsample_rate",                         &config_entry::spo2_upsample_rate, 0.0 },
            { "hr_raw_rate",                                &config_entry::hr_raw_rate, 0.0 },
            { "hr_upsample_rate",                           &config_entry::hr_upsample_rate, 0.0 },
            { "dhr_raw_rate",                               &config_entry::dhr_raw_rate, 0.0 },
            { "dhr_upsample_rate",                          &config_entry::dhr_upsample_rate, 0.0 },
            { "sleepstate_length",                          &config_entry::sleepstate_length, 0.0 },
            { "blanking_period",                            &config_entry::blanking_period, 0.0 },
            { "threshold",                                  &config_entry::threshold, 0.0 },
            { "bin_size_minutes",                           &config_entry::bin_size_minutes, 0.0 },
            { "morph_threshold_ecg",                        &config_entry::morph_threshold_ecg, 0.0 },
            { "morph_threshold_ppg",                        &config_entry::morph_threshold_ppg, 0.0 },
            { "ppg_pearson_threshold",                        &config_entry::ppg_pearson_threshold, 0.0 },
            { "ecg_tukey_fence",                            &config_entry::ecg_tukey_fence, 1.5 },
            { "ppg_tukey_fence",                            &config_entry::ppg_tukey_fence, 1.5 },
            { "region_around_Rpeak_for_morphology_split",   &config_entry::region_around_Rpeak_for_morphology_split, 0.0 },
            { "region_around_PPGPeak_for_morphology_split", &config_entry::region_around_PPGPeak_for_morphology_split, 0.0 },
        };
        for (const NumericKey& n : kNumeric)
            out.*n.field = stod_or_default(get_value_from_config(n.key), n.whenBlank);
        out.min_beats_template_ecg = static_cast<int>(stod_or_default(get_value_from_config("min_beats_template_ecg"), 0));
        out.min_beats_template_ppg = static_cast<int>(stod_or_default(get_value_from_config("min_beats_template_ppg"), 0));

        out.input_path = get_value_from_config("original_file_path");
        out.output_path = get_value_from_config("output_folder");
        out.cgm_folder = get_value_from_config("cgm_folder");
        if (!out.cgm_folder.empty() && out.dataset_type != "BITTIUM")
            std::cerr << "NOTE: cgm_folder is set for " << out.dataset_type
            << "; CGM is processed for BITTIUM only, ignoring it\n";
        out.cgm_bin_minutes = stod_or_default(get_value_from_config("cgm_bin_minutes"), 60.0);
        out.use_consensus_rpeak = parseBool(get_value_from_config("use_consensus_rpeak"), true);
        out.exclude_beat_after_ectopic = parseBool(get_value_from_config("exclude_beat_after_ectopic"), true);
        out.override_morphology = parseBool(get_value_from_config("override_morphology"), false);
        out.notch_filter_hz = stod_or_default(get_value_from_config("notch_filter_hz"), 0.0); //the spec limits notch filter to 0 (none) 50, or 60
        if (out.notch_filter_hz != 0.0 &&
            out.notch_filter_hz != 50.0 &&
            out.notch_filter_hz != 60.0) {
            std::cerr << "WARNING: notch_filter_hz=" << out.notch_filter_hz
                << " is not 0, 50 or 60; disabling notch filter\n";
            out.notch_filter_hz = 0.0;
        }
        out.waveform_highpass_hz = stod_or_default(get_value_from_config("waveform_highpass_hz"), 0.0);

        // --- Subject demographics (stored only, ignored downstream for now) ---
        out.age = stod_or_default(get_value_from_config("age"), 0);
        out.sex = get_value_from_config("sex");
        out.weight_kg = stod_or_default(get_value_from_config("weight_kg"), 0.0);
        out.height_cm = stod_or_default(get_value_from_config("height_cm"), 0.0);
        out.hr_rest = stod_or_default(get_value_from_config("hr_rest"), 0.0);
        out.hr_max = stod_or_default(get_value_from_config("hr_max"), 0.0);

        apply_dataset_specific_channel_labels(out);

        bool ok = prompt_for_missing_folders(out);
        if (ok) create_output_folders(out);
        return ok;
    }
}