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

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <format>
#include <vector>
#include <optional>


static const std::string CONFIG_PATH = "config.csv";

namespace {
    std::vector<std::string> parse_csv_row(const std::string& line) {
        std::vector<std::string> fields;
        std::string cur;
        for (char c : line) {
            if (c == ',') { fields.push_back(cur); cur.clear(); }
            else          cur += c;
        }
        fields.push_back(cur);
        for (auto& f : fields) {
            int first = f.find_first_not_of(" \t\r\n");
            if (first == std::string::npos) { f.clear(); continue; }
            int last = f.find_last_not_of(" \t\r\n");
            f = f.substr(first, last - first + 1);
        }
        return fields;
    }

    void apply_dataset_specific_channel_labels(config_entry& cfg) {
        /*
            Each dataset uses different names for the same channel. For example,
            the EKG channel in MESA is equivalent to the NLS_NOM_ECG_ELEC_POTL_I
            label in CHAOS and the ECG_1 label in Bittium. Each channel is a
            feature of the config_entry cfg, and they are assigned here.
        */
        if (cfg.dataset_type == "MESA") {
            cfg.ecg1Label = "EKG";
            cfg.ppgLabel = "Pleth";
        }
        else if (cfg.dataset_type == "BITTIUM") {
            cfg.ecg1Label = "ECG_1";
            cfg.ecg2Label = "ECG_2";
            cfg.ecg3Label = "ECG_3";
            cfg.accelXLabel = "Accelerometer_X";
            cfg.accelYLabel = "Accelerometer_Y";
            cfg.accelZLabel = "Accelerometer_Z";
        }
        else if (cfg.dataset_type == "CHAOS") {
            cfg.ecg1Label = "NLS_NOM_ECG_ELEC_pOTL_I";
            cfg.ecg2Label = "NLS_NOM_ECG_ELEC_pOTL_II";
            cfg.ecg3Label = "NLS_NOM_ECG_ELEC_pOTL_III";
            cfg.ppgLabel = "NLS_NOM_PULS_OXIM_PLETH";
            cfg.eeg1Label = "NLS_EEG_NAMES_EEG_CHAN1";
            cfg.eeg2Label = "NLS_EEG_NAMES_EEG_CHAN2";
            cfg.eeg3Label = "NLS_EEG_NAMES_EEG_CHAN3";
            cfg.eeg4Label = "NLS_EEG_NAMES_EEG_CHAN4";
            cfg.cvpLabel = "NLS_NOM_PRESS_BLD_VEN_CENT";
            cfg.respLabel = "NLS_NOM_RESP";
            cfg.abpLabel = "NLS_NOM_PRESS_BLD_ART_ABP";
            cfg.artLabel = "NLS_NOM_PRESS_BLD_ART";
            cfg.artPulmLabel = "NLS_NOM_PRESS_BLD_ART_PULM";
        }
    }

    void deriveSubpaths(config_entry& cfg) {
        /*
            All the types of output are stored in their own subfolder of output_path, which are defined in this function
        */
        cfg.snapshot_path = cfg.output_path + "/saved_plot_snapshots/";
        cfg.annealed_data_path = cfg.output_path + "/annealed_output/";
        cfg.noise_data_path = cfg.output_path + "/noise_marking_output/";
        cfg.r_peak_data_path = cfg.output_path + "/r_peak_finding_output/";
        cfg.template_path = cfg.output_path + "/template_outputs/";
        cfg.qtvi_marker_path = cfg.output_path + "/qtvi_marker_path/";

        std::filesystem::create_directories(cfg.snapshot_path);
        std::filesystem::create_directories(cfg.noise_data_path);
        std::filesystem::create_directories(cfg.annealed_data_path);
        std::filesystem::create_directories(cfg.r_peak_data_path);
        std::filesystem::create_directories(cfg.template_path);
        std::filesystem::create_directories(cfg.qtvi_marker_path);

    }

    bool manually_select_folder(config_entry& cfg) {
        /*
            If a folder is not in the config.csv (i.e. its field is empty),
            prompt the user to select it. Fields already populated by
            load_config() are left alone -- no prompt, no overwrite.
        */
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

std::optional<config_entry> load_config(int dataType) {
    /*
    * Loads the config.csv (or whatever is in CONFIG_PATH) and fills up a
    * config_entry struct. The config.csv only has an output_folder; the
    * sub-output paths (template_path, noise_data_path, etc.) are derived
    * from output_path by deriveSubpaths().
    */
    std::ifstream file(CONFIG_PATH);
    if (!file.is_open()) {
        std::cerr << "ERROR: cannot open " << CONFIG_PATH << "\n";
        return std::nullopt;
    }

    std::string input_file_type = (dataType == 1) ? "MESA"
        : (dataType == 2) ? "BITTIUM"
        : (dataType == 3) ? "CHAOS" : "";

    std::string line;
    std::getline(file, line);   // skip header

    while (std::getline(file, line)) {
        auto row = parse_csv_row(line);
        std::string rowType = row[0];
        std::transform(rowType.begin(), rowType.end(), rowType.begin(), ::toupper);
        if (rowType != input_file_type) continue;
        config_entry cfg;

        cfg.dataset_type = rowType;
        cfg.original_file_extention = row[1];
        cfg.sleep_file_extention = row[2];

        // Some CSV cells are intentionally blank (e.g. Bittium has no PPG
        // rate). std::stod throws on empty strings, so guard.
        auto stod_or_zero = [](const std::string& s) -> double {
            if (s.empty()) return 0.0;
            try { return std::stod(s); }
            catch (...) { return 0.0; }
            };
        cfg.ecg_rate = stod_or_zero(row[3]);
        cfg.ppg_rate = stod_or_zero(row[4]);
        cfg.central_venous_pressure_rate = stod_or_zero(row[5]);
        cfg.arterial_blood_pressure_rate = stod_or_zero(row[6]);
        cfg.resp_rate = stod_or_zero(row[7]);
        cfg.target_sampling_rate = stod_or_zero(row[8]);
        cfg.blanking_period = std::stod(row[9]);
		cfg.height_threshold_percent = std::stod(row[10]);
        cfg.bin_length_minutes = stod_or_zero(row[11]);
        cfg.input_path = row[12];
        cfg.output_path = row[13];
        deriveSubpaths(cfg);

        apply_dataset_specific_channel_labels(cfg);
        if (!manually_select_folder(cfg)) return std::nullopt;
        return cfg;
    }
    return std::nullopt;
}