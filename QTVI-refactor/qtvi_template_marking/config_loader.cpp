/**
 * @file   config_loader.cpp
 * @brief  Loads the merged config.csv and walks the source-file directory.
 *         The marking GUI's main() pulls one .bin at a time from make_binfile
 *         and doesn't need to know whether it was just produced or cached.
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

static const std::string CONFIG_PATH = "config.csv";

namespace {

    std::vector<std::string> parseCsvRow(const std::string& line) {
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

    // Channel labels per dataset.
    void apply_dataset_specific_channel_labels(config_entry& cfg) {
        /*
            Each dataset uses different names for the same channel. For example,
            the EKG channel in MESA is equivalent to the NLS_NOM_ECG_ELEC_pOTL_I
            label in CHAOS and the ECG_1 label in Bittium. Each channel is a
            feature of the config_entry cfg, and they are assigned here.
        */
        if (cfg.dataType == "MESA") {
            cfg.ecg1Label = "EKG";
            cfg.ppgLabel = "Pleth";
        }
        else if (cfg.dataType == "BITTIUM") {
            cfg.ecg1Label = "ECG_1";
            cfg.ecg2Label = "ECG_2";
            cfg.ecg3Label = "ECG_3";
            cfg.accelXLabel = "Accelerometer_X";
            cfg.accelYLabel = "Accelerometer_Y";
            cfg.accelZLabel = "Accelerometer_Z";
        }
        else if (cfg.dataType == "CHAOS") {
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

    // Recompute every output subfolder from cfg.output_path. Called once
    // from load_config (so the defaults are populated even if the prompt
    // is skipped) and once from promptForMissingPaths (so a user-selected
    // output_path actually propagates to its derived paths). Without the
    // second call, the subpaths stay stuck at whatever load_config first
    // computed -- which is the wrong value (or empty) whenever the CSV's
    // output_folder column was blank.
    void deriveSubpaths(config_entry& cfg) {
        cfg.bin_file_path = cfg.output_path + "/input_binfile/";
        cfg.noise_data_path = cfg.output_path + "/noise_marking_output/";
        cfg.annealed_data_path = cfg.output_path + "/annealed_output/";
        cfg.r_peak_data_path = cfg.output_path + "/r_peak_finding_output/";
        cfg.template_path = cfg.output_path + "/template_path/";
        cfg.qtvi_marker_path = cfg.output_path + "/qtvi_marker_path/";
        cfg.qtvi_data_path = cfg.output_path + "/qtvi_data_path/";
    }

}   // anonymous namespace


bool load_config(int dataType, config_entry& out) {
    /*
    * Loads the config.csv (or whatever is in CONFIG_PATH) and fills up a
    * config_entry struct. The config.csv only has an output_folder; the
    * sub-output paths (template_path, noise_data_path, etc.) are derived
    * from output_path by deriveSubpaths().
    */
    std::ifstream file(CONFIG_PATH);
    if (!file.is_open()) {
        std::cerr << "ERROR: cannot open " << CONFIG_PATH << "\n";
        return false;
    }

    std::string input_file_type = (dataType == 1) ? "MESA"
        : (dataType == 2) ? "BITTIUM"
        : (dataType == 3) ? "CHAOS" : "";

    std::string line;
    std::getline(file, line);   // skip header

    while (std::getline(file, line)) {
        auto row = parseCsvRow(line);
        if (row.size() < 9) continue;

        std::string rowType = row[0];
        std::transform(rowType.begin(), rowType.end(), rowType.begin(), ::toupper);
        if (rowType != input_file_type) continue;

        out.dataType = rowType;
        out.mainExt = row[1];
        out.sleepExt = row[2];

        // Some CSV cells are intentionally blank (e.g. Bittium has no PPG
        // rate). std::stod throws on empty strings, so guard.
        auto stod_or_zero = [](const std::string& s) -> double {
            if (s.empty()) return 0.0;
            try { return std::stod(s); }
            catch (...) { return 0.0; }
            };
        out.ecgRate = stod_or_zero(row[3]);
        out.ppgRate = stod_or_zero(row[4]);
        out.finalSamplingRate = stod_or_zero(row[5]);
        out.bin_length_minutes = stod_or_zero(row[6]);

        out.output_path = row[8];
        deriveSubpaths(out);

        apply_dataset_specific_channel_labels(out);
        return true;
    }
    return false;
}

bool promptForMissingPaths(config_entry& cfg) {
    // Each entry is (label shown in the dialog title, pointer to the field
    // to fill). Only empty fields prompt the user.
    const std::vector<std::pair<const char*, std::string*>> fields = {
        { "Output Folder of the noise marking script", &cfg.output_path },
    };

    for (const auto& [label, fieldPtr] : fields) {
        if (!fieldPtr->empty()) continue;

        QString title = QString("Select folder for %1 (%2)")
            .arg(label, QString::fromStdString(cfg.dataType));
        QString chosen = QFileDialog::getExistingDirectory(
            nullptr, title, QString(),
            QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
        *fieldPtr = chosen.toStdString();
    }

    // Any change to output_path above means every derived subpath is now
    // stale. Recompute them so callers see the correct, current values.
    deriveSubpaths(cfg);
    return true;
}