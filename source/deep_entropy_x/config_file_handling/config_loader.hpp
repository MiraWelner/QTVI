#pragma once
/**
 * @file   config_loader.hpp
 * @brief  Loads the config.csv file, parses it based on the datset type selected by the user, and fills up a config_entry struct with the relevant paths, rates, and channel labels.
 *         The channel labels (eg. "ECG_1" vs "EKG") are dataset-specific but not in the config file, so they are assigned in apply_dataset_specific_channel_labels() based on the dataset type.
 *         The output paths are either found in the config file, or prompted for manually if the config file cells are blank. The output_path is used to create the subfolders where the
 *         specific types of output are found.
 
 */

#include "config.hpp"

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
            cfg.eeg_3_label = "EEG3 ";
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
        cfg.quality_metric = create_subfolder("quality_metric");
        cfg.training_log = create_subfolder("training_log");
        cfg.snapshot_path = create_subfolder("snapshot_path");
        cfg.vcg_output = create_subfolder("vcg_output");
    }

    inline bool manually_select_folder(config_entry& cfg) {
        // If the input or output folder is not in the config.csv (i.e. its field is empty), prompt the user to select it.
        if (cfg.bin_file_path.empty() && !cfg.input_path.empty())
            cfg.bin_file_path = cfg.input_path;
        const std::vector<std::pair<const char*, std::string*>> fields = {
            { "Bin Files:", &cfg.bin_file_path },
            { "Output", &cfg.output_path },
        };

        bool outputChanged = false;
        for (const auto& [label, fieldPtr] : fields) {
            if (!fieldPtr->empty()) continue;

            QString title = QString("%1 (%2)").arg(label, QString::fromStdString(cfg.dataset_type));
            QString chosen = QFileDialog::getExistingDirectory(
                nullptr, title, QString(),
                QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
            if (chosen.isEmpty()) return false;
            *fieldPtr = chosen.toStdString();
            if (fieldPtr == &cfg.output_path) outputChanged = true;
        }

        if (outputChanged) create_output_folders(cfg);
        return true;
    }
}

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

    // Read the header row so column names can be mapped to indices.
    std::string header;
    if (!std::getline(file, header)) {
        std::cerr << "ERROR: " << CONFIG_PATH << " is empty\n";
        return false;
    }
    std::vector<std::string> headerFields = parse_config_row(header);
    std::unordered_map<std::string, int> col;
    for (int i = 0; i < (int)headerFields.size(); ++i) {
        col[normalize_key(headerFields[i])] = i;
    }
    std::string line;
    while (std::getline(file, line)) {
        std::vector<std::string> row = parse_config_row(line);
        auto cell = [&](const std::string& name) -> std::string {
            // normalize_key on the LOOKUP too, not just the header. Without
            // it a mixed-case name can never match a lowercased key, and the
            // miss is indistinguishable from a blank cell -- which is how
            // region_around_Rpeak_for_morphology_split silently defaulted to
            // 0 ("split on the whole beat") whatever the config said.
            auto it = col.find(normalize_key(name));
            if (it == col.end() || it->second >= (int)row.size()) return {};
            return row[it->second];
            };
        std::string rowName = cell("data_type");
        std::transform(rowName.begin(), rowName.end(), rowName.begin(), ::toupper);
        if (rowName != user_selected_dataset) continue;


        out.dataset_type = user_selected_dataset;
        out.main_file_extension = cell("main_file_extension");
        out.sleep_file_extension = cell("sleep_file_extension");

        out.ecg_raw_rate = stod_or_default(cell("ecg_raw_rate"), 0.0);
        out.ecg_upsample_rate = stod_or_default(cell("ecg_upsampled_rate"), 0.0);
        out.ppg_raw_rate = stod_or_default(cell("ppg_raw_rate"), 0.0);
        out.ppg_upsample_rate = stod_or_default(cell("ppg_upsampled_rate"), 0.0);
        out.cvp_raw_rate = stod_or_default(cell("cvp_raw_rate"), 0.0);
        out.cvp_upsample_rate = stod_or_default(cell("cvp_upsampled_rate"), 0.0);
        out.pres_raw_rate = stod_or_default(cell("pres_raw_rate"), 0.0);
        out.pres_upsample_rate = stod_or_default(cell("pres_upsampled_rate"), 0.0);
        out.abp_raw_rate = stod_or_default(cell("abp_raw_rate"), 0.0);
        out.abp_upsample_rate = stod_or_default(cell("abp_upsampled_rate"), 0.0);
        out.art_raw_rate = stod_or_default(cell("art_raw_rate"), 0.0);
        out.art_upsample_rate = stod_or_default(cell("art_upsampled_rate"), 0.0);
        out.art_pulm_raw_rate = stod_or_default(cell("art_pulm_raw_rate"), 0.0);
        out.art_pulm_upsample_rate = stod_or_default(cell("art_pulm_upsampled_rate"), 0.0);
        out.accel_raw_rate = stod_or_default(cell("accel_raw_rate"), 0.0);
        out.accel_upsample_rate = stod_or_default(cell("accel_upsampled_rate"), 0.0);
        out.temp_raw_rate = stod_or_default(cell("temp_raw_rate"), 0.0);
        out.temp_upsample_rate = stod_or_default(cell("temp_upsampled_rate"), 0.0);
        out.marker_raw_rate = stod_or_default(cell("marker_raw_rate"), 0.0);
        out.marker_upsample_rate = stod_or_default(cell("marker_upsampled_rate"), 0.0);
        out.resp_raw_rate = stod_or_default(cell("resp_raw_rate"), 0.0);
        out.resp_upsample_rate = stod_or_default(cell("resp_upsampled_rate"), 0.0);
        out.pacemaker_raw_rate = stod_or_default(cell("pacemaker_event_raw_rate"), 0.0);
        out.pacemaker_upsample_rate = stod_or_default(cell("pacemaker_event_upsampled_rate"), 0.0);
        out.eeg_raw_rate = stod_or_default(cell("eeg_raw_rate"), 0.0);
        out.eeg_upsample_rate = stod_or_default(cell("eeg_upsampled_rate"), 0.0);
        out.eog_l_raw_rate = stod_or_default(cell("eogl_raw_rate"), 0.0);
        out.eog_l_upsample_rate = stod_or_default(cell("eogl_upsampled_rate"), 0.0);
        out.eog_r_raw_rate = stod_or_default(cell("eogr_raw_rate"), 0.0);
        out.eog_r_upsample_rate = stod_or_default(cell("eogr_upsampled_rate"), 0.0);
        out.emg_raw_rate = stod_or_default(cell("emg_raw_rate"), 0.0);
        out.emg_upsample_rate = stod_or_default(cell("emg_upsampled_rate"), 0.0);
        out.flow_raw_rate = stod_or_default(cell("flow_raw_rate"), 0.0);
        out.flow_upsample_rate = stod_or_default(cell("flow_upsampled_rate"), 0.0);
        out.snore_raw_rate = stod_or_default(cell("snore_raw_rate"), 0.0);
        out.snore_upsample_rate = stod_or_default(cell("snore_upsampled_rate"), 0.0);
        out.thor_raw_rate = stod_or_default(cell("thor_raw_rate"), 0.0);
        out.thor_upsample_rate = stod_or_default(cell("thor_upsampled_rate"), 0.0);
        out.abdo_raw_rate = stod_or_default(cell("abdo_raw_rate"), 0.0);
        out.abdo_upsample_rate = stod_or_default(cell("abdo_upsampled_rate"), 0.0);
        out.leg_raw_rate = stod_or_default(cell("leg_raw_rate"), 0.0);
        out.leg_upsample_rate = stod_or_default(cell("leg_upsampled_rate"), 0.0);
        out.auxac_raw_rate = stod_or_default(cell("auxac_raw_rate"), 0.0);
        out.auxac_upsample_rate = stod_or_default(cell("auxac_upsampled_rate"), 0.0);
        out.therm_raw_rate = stod_or_default(cell("therm_raw_rate"), 0.0);
        out.therm_upsample_rate = stod_or_default(cell("therm_upsampled_rate"), 0.0);
        out.pos_raw_rate = stod_or_default(cell("pos_raw_rate"), 0.0);
        out.pos_upsample_rate = stod_or_default(cell("pos_upsampled_rate"), 0.0);
        out.oxstatus_raw_rate = stod_or_default(cell("oxstatus_raw_rate"), 0.0);
        out.oxstatus_upsample_rate = stod_or_default(cell("oxstatus_upsampled_rate"), 0.0);
        out.spo2_raw_rate = stod_or_default(cell("spo2_raw_rate"), 0.0);
        out.spo2_upsample_rate = stod_or_default(cell("spo2_upsampled_rate"), 0.0);
        out.hr_raw_rate = stod_or_default(cell("hr_raw_rate"), 0.0);
        out.hr_upsample_rate = stod_or_default(cell("hr_upsampled_rate"), 0.0);
        out.dhr_raw_rate = stod_or_default(cell("dhr_raw_rate"), 0.0);
        out.dhr_upsample_rate = stod_or_default(cell("dhr_upsampled_rate"), 0.0);
        out.sleepstate_length = stod_or_default(cell("sleepstate_length"), 0.0);
        out.blanking_period = stod_or_default(cell("blanking_period"), 0.0);
        out.threshold = stod_or_default(cell("threshold"), 0.0);
        out.bin_size_minutes = stod_or_default(cell("bin_size_minutes"), 0.0);
        out.ecg_match_floor = stod_or_default(cell("ecg_match_floor"), 0.0);
        out.ppg_match_floor = stod_or_default(cell("ppg_match_floor"), 0.0);
        out.ppg_fit_error_pct = stod_or_default(cell("ppg_fit_error_pct"), 0.0);
        out.min_beats_template_ecg = stod_or_default(cell("min_beats_template_ecg"), 0);
        out.min_beats_template_ppg = stod_or_default(cell("min_beats_template_ppg"), 0);
        out.region_around_Rpeak_for_morphology_split = stod_or_default(cell("region_around_Rpeak_for_morphology_split"), 0.0);
        out.region_around_PPGPeak_for_morphology_split = stod_or_default(cell("region_around_PPGPeak_for_morphology_split"), 0.0);
        out.input_path = cell("original_file_path");
        out.output_path = cell("output_folder");
        out.use_consensus_rpeak = parseBool(cell("use_consensus_rpeak"), true);
        out.notch_filter_hz = stod_or_default(cell("notch_filter_hz"), 0); //the spec limits notch filter to 0 (none) 50, or 60
        if (out.notch_filter_hz != 0 &&
            out.notch_filter_hz != 50 &&
            out.notch_filter_hz != 60) {
            std::cerr << "WARNING: notch_filter_hz=" << out.notch_filter_hz
                << " is not 50 or 60; disabling notch filter\n";
            out.notch_filter_hz = 0;
        }
        out.waveform_highpass_hz = stod_or_default(cell("waveform_highpass_hz"), 0.0);

        // --- Subject demographics (stored only, ignored downstream for now) ---
        out.age = stod_or_default(cell("age"), 0);
        out.sex = cell("sex");
        out.weight_kg = stod_or_default(cell("weight_kg"), 0.0);
        out.height_cm = stod_or_default(cell("height_cm"), 0.0);
        out.hr_rest = stod_or_default(cell("hr_rest"), 0);
        out.hr_max = stod_or_default(cell("hr_max"), 0);

        apply_dataset_specific_channel_labels(out);

        bool ok = manually_select_folder(out);
        if (ok) create_output_folders(out);
        return ok;
    }
    std::cerr << "ERROR: no row with data_type=" << user_selected_dataset
        << " in " << CONFIG_PATH << "\n";
    return false;
}