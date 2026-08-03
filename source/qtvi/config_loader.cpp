/**
 * @file   config_loader.cpp
 * @brief  Loads the config.csv file, parses it based on the datset type selected by the user, and fills up a config_entry struct with the relevant paths, rates, and channel labels.
 *         The channel labels (eg. "ECG_1" vs "EKG") are dataset-specific but not in the config file, so they are assigned in apply_dataset_specific_channel_labels() based on the dataset type.
 *         The output paths are either found in the config file, or prompted for manually if the config file cells are blank. The output_path is used to create the subfolders where the
 *         specific types of cfg put are found.
 */

#include "config_loader.hpp"

#include <QDir>
#include <QDirIterator>
#include <QFileInfo>
#include <QFileDialog>
#include <unordered_map>


#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <format>
#include <vector>
#include <optional>


static const std::string CONFIG_PATH = "config.csv";

namespace {
    static double stod_or_zero(const std::string& s) {
        if (s.empty()) return 0.0;
        try { return std::stod(s); }
        catch (...) { return 0.0; }
    }

    static int stoi_or_default(const std::string& s, int dflt) {
        if (s.empty()) return dflt;
        try { return std::stoi(s); }
        catch (...) { return dflt; }
    }

    static double stod_or_default(const std::string& s, double dflt) {
        if (s.empty()) return dflt;
        try { return std::stod(s); }
        catch (...) { return dflt; }
    }


    std::vector<std::string> parse_config_row(const std::string& line) {
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

    static bool parseBool(const std::string& s, bool dflt) {
        //this is for the the r consensus variable in the config.csv file
        if (s.empty()) return dflt;
        std::string t; for (char c : s) t += std::tolower((unsigned char)c);
        if (t == "1" || t == "true" || t == "yes") return true;
        if (t == "0" || t == "false" || t == "no") return false;
        return dflt;
    }


    void apply_dataset_specific_channel_labels(config_entry& cfg) {
        /*
            The ECG and PPG channels exist in different datasets, but their names are
            different so they are set here.
        */
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
    }


    void deriveSubpaths(config_entry & cfg) {
        /*
        This script yeilds different types of outputs, each of which are saved in their own
        subfolder of output_path          
        */
        cfg.snapshot_path = cfg.output_path + "/saved_plot_snapshots/";
        cfg.annealed_data_path = cfg.output_path + "/annealed_output/";
        cfg.noise_data_path = cfg.output_path + "/noise_marking_output/";
        cfg.r_peak_data_path = cfg.output_path + "/r_peak_finding_output/";
        cfg.template_path = cfg.output_path + "/template_outputs/";
        cfg.qtvi_marker_path = cfg.output_path + "/qtvi_marker_path/";
        cfg.quality_metric = cfg.output_path + "/quality_metric/";

        cfg.snapshot_path = cfg.output_path + "/snapshot_path/";

    }

    bool manually_select_folder(config_entry& cfg) {
        /*
            If a folder is not in the config.csv (i.e. its field is empty),
            prompt the user to select it.
        */
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

        if (outputChanged) deriveSubpaths(cfg);
        return true;
    }
}

bool load_config(int dataType, config_entry& out) {
    /*
        Take a pointer to a config_entry and fill it with the values from the config.csv file. The return
        bool indicates if it happened correctly, the dataType is an int that indicates which dataset to load (1 = MESA, 2 = BITTIUM, 3 = CHAOS)
    */
    std::ifstream file(CONFIG_PATH);
    if (!file.is_open()) {
        std::cerr << "ERROR: cannot open " << CONFIG_PATH << "\n";
        return false;
    }
    std::string user_selected_dataset = (dataType == 1) ? "MESA"
        : (dataType == 2) ? "BITTIUM"
        : (dataType == 3) ? "CHAOS" : "";

    // Read the header row so column names can be mapped to indices.
    std::string header;
    if (!std::getline(file, header)) {
        std::cerr << "ERROR: " << CONFIG_PATH << " is empty\n";
        return false;
    }
    std::vector<std::string> headerFields = parse_config_row(header);
    std::unordered_map<std::string, int> col;
    for (int i = 0; i < (int)headerFields.size(); ++i) {
        std::string h = headerFields[i];
        std::transform(h.begin(), h.end(), h.begin(), ::tolower);
        col[h] = i;
    }
    std::string line;
    while (std::getline(file, line)) {
        std::vector<std::string> row = parse_config_row(line);
        auto cell = [&](const std::string& name) -> std::string {
            auto it = col.find(name);
            if (it == col.end() || it->second >= (int)row.size()) return {};
            return row[it->second];
            };
        std::string rowName = cell("data_type");
        std::transform(rowName.begin(), rowName.end(), rowName.begin(), ::toupper);
        if (rowName != user_selected_dataset) continue;

        
        out.dataset_type = user_selected_dataset;
        out.main_file_extention = cell("main_file_extention");
        out.sleep_file_extention = cell("sleep_file_extention");

        out.ecg_raw_rate = stod_or_zero(cell("ecg_raw_rate"));
        out.ecg_upsample_rate = stod_or_zero(cell("ecg_upsampled_rate"));
        out.ppg_raw_rate = stod_or_zero(cell("ppg_raw_rate"));
        out.ppg_upsample_rate = stod_or_zero(cell("ppg_upsampled_rate"));
        out.cvp_raw_rate = stod_or_zero(cell("cvp_raw_rate"));
        out.cvp_upsample_rate = stod_or_zero(cell("cvp_upsampled_rate"));
        out.pres_raw_rate = stod_or_zero(cell("pres_raw_rate"));
        out.pres_upsample_rate = stod_or_zero(cell("pres_upsampled_rate"));
        out.abp_raw_rate = stod_or_zero(cell("abp_raw_rate"));
        out.abp_upsample_rate = stod_or_zero(cell("abp_upsampled_rate"));
        out.art_raw_rate = stod_or_zero(cell("art_raw_rate"));
        out.art_upsample_rate = stod_or_zero(cell("art_upsampled_rate"));
        out.art_pulm_raw_rate = stod_or_zero(cell("art_pulm_raw_rate"));
        out.art_pulm_upsample_rate = stod_or_zero(cell("art_pulm_upsampled_rate"));
        out.accel_raw_rate = stod_or_zero(cell("accel_raw_rate"));
        out.accel_upsample_rate = stod_or_zero(cell("accel_upsampled_rate"));
        out.temp_raw_rate = stod_or_zero(cell("temp_raw_rate"));
        out.temp_upsample_rate = stod_or_zero(cell("temp_upsampled_rate"));
        out.marker_raw_rate = stod_or_zero(cell("marker_raw_rate"));
        out.marker_upsample_rate = stod_or_zero(cell("marker_upsampled_rate"));
        out.resp_raw_rate = stod_or_zero(cell("resp_raw_rate"));
        out.resp_upsample_rate = stod_or_zero(cell("resp_upsampled_rate"));
        out.pacemaker_raw_rate = stod_or_zero(cell("pacemaker_event_raw_rate"));
        out.pacemaker_upsample_rate = stod_or_zero(cell("pacemaker_event_upsampled_rate"));
        out.eeg_raw_rate = stod_or_zero(cell("eeg_raw_rate"));
        out.eeg_upsample_rate = stod_or_zero(cell("eeg_upsampled_rate"));
        out.eog_l_raw_rate = stod_or_zero(cell("eogl_raw_rate"));
        out.eog_l_upsample_rate = stod_or_zero(cell("eogl_upsampled_rate"));
        out.eog_r_raw_rate = stod_or_zero(cell("eogr_raw_rate"));
        out.eog_r_upsample_rate = stod_or_zero(cell("eogr_upsampled_rate"));
        out.emg_raw_rate = stod_or_zero(cell("emg_raw_rate"));
        out.emg_upsample_rate = stod_or_zero(cell("emg_upsampled_rate"));
        out.flow_raw_rate = stod_or_zero(cell("flow_raw_rate"));
        out.flow_upsample_rate = stod_or_zero(cell("flow_upsampled_rate"));
        out.snore_raw_rate = stod_or_zero(cell("snore_raw_rate"));
        out.snore_upsample_rate = stod_or_zero(cell("snore_upsampled_rate"));
        out.thor_raw_rate = stod_or_zero(cell("thor_raw_rate"));
        out.thor_upsample_rate = stod_or_zero(cell("thor_upsampled_rate"));
        out.abdo_raw_rate = stod_or_zero(cell("abdo_raw_rate"));
        out.abdo_upsample_rate = stod_or_zero(cell("abdo_upsampled_rate"));
        out.leg_raw_rate = stod_or_zero(cell("leg_raw_rate"));
        out.leg_upsample_rate = stod_or_zero(cell("leg_upsampled_rate"));
        out.auxac_raw_rate = stod_or_zero(cell("auxac_raw_rate"));
        out.auxac_upsample_rate = stod_or_zero(cell("auxac_upsampled_rate"));
        out.therm_raw_rate = stod_or_zero(cell("therm_raw_rate"));
        out.therm_upsample_rate = stod_or_zero(cell("therm_upsampled_rate"));
        out.pos_raw_rate = stod_or_zero(cell("pos_raw_rate"));
        out.pos_upsample_rate = stod_or_zero(cell("pos_upsampled_rate"));
        out.oxstatus_raw_rate = stod_or_zero(cell("oxstatus_raw_rate"));
        out.oxstatus_upsample_rate = stod_or_zero(cell("oxstatus_upsampled_rate"));
        out.spo2_raw_rate = stod_or_zero(cell("spo2_raw_rate"));
        out.spo2_upsample_rate = stod_or_zero(cell("spo2_upsampled_rate"));
        out.hr_raw_rate = stod_or_zero(cell("hr_raw_rate"));
        out.hr_upsample_rate = stod_or_zero(cell("hr_upsampled_rate"));
        out.dhr_raw_rate = stod_or_zero(cell("dhr_raw_rate"));
        out.dhr_upsample_rate = stod_or_zero(cell("dhr_upsampled_rate"));
        out.sleepstate_length = stod_or_zero(cell("sleepstate_length"));
        out.blanking_period = stod_or_zero(cell("blanking_period"));
        out.threshold = stod_or_zero(cell("threshold"));
        out.bin_size_minutes = stod_or_zero(cell("bin_size_minutes"));
        out.input_path = cell("original_file_path");
        out.output_path = cell("output_folder");
        out.use_consensus_rpeak = parseBool(cell("use_consensus_rpeak"), true);
        // --- Filtering options ---
        // Notch: blank -> 0 (disabled). Only 0/50/60 are valid; anything else
        // warns and falls back to disabled.
        out.notch_filter_hz = stoi_or_default(cell("notch_filter_hz"), 0);
        if (out.notch_filter_hz != 0 &&
            out.notch_filter_hz != 50 &&
            out.notch_filter_hz != 60) {
            std::cerr << "WARNING: notch_filter_hz=" << out.notch_filter_hz
                << " is not 50 or 60; disabling notch filter\n";
            out.notch_filter_hz = 0;
        }

        // High-pass: blank -> documented default 0.5 (i.e. enabled at 0.5 Hz).
        // An explicit 0 in the cell means disabled.
        // NOTE: if a blank cell should instead mean *disabled*, change the
        // default below from 0.5 to 0.0.
        out.waveform_highpass_hz = stod_or_default(cell("waveform_highpass_hz"), 0.5);

        // --- Subject demographics (stored only, ignored downstream for now) ---
        out.age = stoi_or_default(cell("age"), 0);
        out.sex = cell("sex");                       // stored verbatim
        out.weight_kg = stod_or_default(cell("weight_kg"), 0.0);
        out.height_cm = stod_or_default(cell("height_cm"), 0.0);
        out.hr_rest = stoi_or_default(cell("hr_rest"), 0);
        out.hr_max = stoi_or_default(cell("hr_max"), 0);

        apply_dataset_specific_channel_labels(out);

        bool ok = manually_select_folder(out);
        if (ok) deriveSubpaths(out);   // always derive, whether paths came from
        // the config or the manual prompt
        return ok;
    }
    return false;
}