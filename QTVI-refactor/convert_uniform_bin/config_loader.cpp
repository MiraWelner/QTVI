/**
 * @file   config_loader.cpp
 * @brief  Loads the config.csv and extracts the feilds that are useful for creating the bin file
 */

#include "config_loader.hpp"
#include "file_to_bin.hpp"

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


std::vector<std::string> parse_config_row(const std::string& line) {
    /*
        The config file is a csv - this just is a util for loading a .csv row
    */
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


void apply_dataset_specific_channel_labels(config_entry& cfg) {
    /*
        Each dataset uses different names for the same channel. For example, the EKG channel in mesa is equivalent
        to the NLS_NOM_ECG_ELEC_pOTL_I label in CHAOS and the ECG_1 label in bittium. Each channel is a feature of the config_entry
        cfg, and they are assigned here.
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
        cfg.ecg1Label = "NLS_NOM_ECG_ELEC_POTL_I";
        cfg.ecg2Label = "NLS_NOM_ECG_ELEC_POTL_II";
        cfg.ecg3Label = "NLS_NOM_ECG_ELEC_POTL_III";
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


bool load_config(int dataType, config_entry& out) {
    /*
    * Loads the config.csv (or whatever is in CONFIG_PATH) and fills up a config_entry
    * struct with all the data. Note that this particular function does NOT load the
    * file paths, because those may be empty. Those are handled in build_derived_paths
    */
    std::ifstream file(CONFIG_PATH);
    if (!file.is_open()) {
        std::cerr << "ERROR: cannot open " << CONFIG_PATH << "\n";
        return false;
    }

    std::string input_file_type = (dataType == 1) ? "MESA" : (dataType == 2) ? "BITTIUM" : (dataType == 3) ? "CHAOS" : "";

    std::string line;

    std::getline(file, line);   // skip header

    while (std::getline(file, line)) {
        auto row = parse_config_row(line);
        // 12 columns: data_type, main_ext, sleep_ext, ecg_rate, ppg_rate,
        // cvp_rate, abp_rate, resp_rate, upsampled_rate, bin_size_minutes,
        // original_file_path, output_folder.
        if (row.size() < 12) continue;

        std::string rowType = row[0];
        std::transform(rowType.begin(), rowType.end(), rowType.begin(), ::toupper);
        if (rowType != input_file_type) continue;

        // Some cells are intentionally blank (e.g. MESA has no cvp_rate).
        // std::stod throws on empty strings, so guard.
        auto stod_or_zero = [](const std::string& s) -> double {
            if (s.empty()) return 0.0;
            try { return std::stod(s); }
            catch (...) { return 0.0; }
            };

        out.dataType = rowType;
        out.mainExt = row[1];
        out.sleepExt = row[2];
        out.ecgRate = stod_or_zero(row[3]);
        out.ppgRate = stod_or_zero(row[4]);
        out.cvpRate = stod_or_zero(row[5]);
        out.abpRate = stod_or_zero(row[6]);
        out.respRate = stod_or_zero(row[7]);
        out.finalSamplingRate = stod_or_zero(row[8]);
        out.bin_length_minutes = stod_or_zero(row[9]);
        out.input_path = row[10];
        out.output_path = row[11];

        apply_dataset_specific_channel_labels(out);
        return true;
    }
    return false;
}

bool promptForMissingPaths(config_entry& cfg) {
    /*
            If the path sections of the config file are empty, prompt the user
            to navigate to the paths via the QT GUI which has a built in file navigator
    */
    const std::vector<std::pair<const char*, std::string*>> fields = {
        { "Raw EDF or DAT files", &cfg.input_path },
        { "Output", &cfg.output_path },
    };

    for (const auto& [label, fieldPtr] : fields) {
        if (!fieldPtr->empty()) continue;

        QString title = QString("%1 (%2)")
            .arg(label, QString::fromStdString(cfg.dataType));
        QString chosen = QFileDialog::getExistingDirectory(
            nullptr, title, QString(),
            QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
        if (chosen.isEmpty()) return false;
        *fieldPtr = chosen.toStdString();
    }
    return true;
}

QStringList discoverSourceFiles(const config_entry& cfg) {
    /*
        Take the folder full of .bin, .dat, or .edf (case insensitive) and run through it, return Qstrings of all the files in it
    */
    std::string extention = cfg.mainExt;
    std::transform(extention.begin(), extention.end(), extention.begin(), ::toupper);
    QString extention_qstring = QString::fromStdString(extention);

    QDirIterator it(QString::fromStdString(cfg.input_path),
        { "*" + extention_qstring.toLower(), "*" + extention_qstring, "*.bin" },
        QDir::Files, QDirIterator::Subdirectories);

    QStringList out;
    while (it.hasNext()) out << it.next();
    out.sort();
    return out;
}