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
            The ECG and PPG channels exist in different datasets, but their names are
            different so they are set here.
        */
        if (cfg.dataset_type == "MESA") {
            cfg.ecg_1_label = "EKG";
            cfg.ppg_label = "Pleth";
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

    // Column-name -> index from the header (survives reordering / new columns).
    std::string header;
    if (!std::getline(file, header)) return std::nullopt;
    std::vector<std::string> headerFields = parse_csv_row(header);
    std::unordered_map<std::string, int> col;
    for (int i = 0; i < (int)headerFields.size(); ++i) {
        std::string h = headerFields[i];
        std::transform(h.begin(), h.end(), h.begin(), ::tolower);
        col[h] = i;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::vector<std::string> row = parse_csv_row(line);

        auto cell = [&](const std::string& name) -> std::string {
            auto it = col.find(name);
            if (it == col.end() || it->second >= (int)row.size()) return {};
            return row[it->second];
            };

        std::string rowType = cell("data_type");
        std::transform(rowType.begin(), rowType.end(), rowType.begin(), ::toupper);
        if (rowType != input_file_type) continue;

        auto stod_or_zero = [](const std::string& s) -> double {
            if (s.empty()) return 0.0;
            try { return std::stod(s); }
            catch (...) { return 0.0; }
            };

        config_entry cfg;
        cfg.dataset_type = rowType;
        cfg.original_file_extention = cell("main_file_extention");
        cfg.sleep_file_extention = cell("sleep_file_extention");
        cfg.ecg_rate = stod_or_zero(cell("ecg_rate"));
        cfg.ppg_rate = stod_or_zero(cell("ppg_rate"));
        cfg.central_venous_pressure_rate = stod_or_zero(cell("cvp_rate"));
        cfg.arterial_blood_pressure_rate = stod_or_zero(cell("abp_rate"));
        cfg.resp_rate = stod_or_zero(cell("resp_rate"));
        cfg.target_sampling_rate = stod_or_zero(cell("upsampled_rate"));
        cfg.blanking_period = stod_or_zero(cell("blanking_period"));
        cfg.height_threshold_percent = stod_or_zero(cell("threshold"));
        cfg.bin_length_minutes = stod_or_zero(cell("bin_size_minutes"));
        cfg.input_path = cell("original_file_path");
        cfg.output_path = cell("output_folder");
        deriveSubpaths(cfg);

        apply_dataset_specific_channel_labels(cfg);
        if (!manually_select_folder(cfg)) return std::nullopt;
        return cfg;
    }
    return std::nullopt;
}