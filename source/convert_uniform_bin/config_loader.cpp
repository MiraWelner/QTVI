/**
 * @file   config_loader.cpp
 * @brief  Loads the config.csv and extracts the fields that are useful for creating the bin file
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


bool load_config(int dataType, config_entry& out) {
    std::ifstream file(CONFIG_PATH);
    if (!file.is_open()) {
        std::cerr << "ERROR: cannot open " << CONFIG_PATH << "\n";
        return false;
    }

    std::string user_selected_dataset = (dataType == 1) ? "MESA"
        : (dataType == 2) ? "BITTIUM"
        : (dataType == 3) ? "CHAOS" : "";

    // Build column-name -> index from the header so reordered / added columns
    // in config.csv don't break parsing. Matched case-insensitively.
    std::string header;
    if (!std::getline(file, header)) return false;
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

        // Cell by column name; "" if the column is missing or the row is short.
        auto cell = [&](const std::string& name) -> std::string {
            auto it = col.find(name);
            if (it == col.end() || it->second >= (int)row.size()) return {};
            return row[it->second];
            };

        std::string rowName = cell("data_type");
        std::transform(rowName.begin(), rowName.end(), rowName.begin(), ::toupper);
        if (rowName != user_selected_dataset) continue;

        auto stod_or_zero = [](const std::string& s) -> double {
            if (s.empty()) return 0.0;
            try { return std::stod(s); }
            catch (...) { return 0.0; }
            };

        out.dataset_type = user_selected_dataset;
        out.original_file_extention = cell("main_file_extention");
        out.sleep_file_extention = cell("sleep_file_extention");
        out.ecg_rate = stod_or_zero(cell("ecg_rate"));
        out.ppg_rate = stod_or_zero(cell("ppg_rate"));
        out.central_venous_pressure_rate = stod_or_zero(cell("cvp_rate"));
        out.arterial_blood_pressure_rate = stod_or_zero(cell("abp_rate"));
        out.accel_rate = stod_or_zero(cell("accel_rate"));
        out.temp_rate = stod_or_zero(cell("temp_rate"));
        out.marker_rate = stod_or_zero(cell("marker_rate"));
        out.resp_rate = stod_or_zero(cell("resp_rate"));
		out.pacemaker_event_rate = stod_or_zero(cell("pacemaker_event_rate"));
        out.high_upsample_rate = stod_or_zero(cell("high_upsampled_rate"));
        out.low_sample_rate = stod_or_zero(cell("low_upsampling_rate"));
        out.sleep_state_length = stod_or_zero(cell("sleepstate_length"));
        out.input_path = cell("original_file_path");
        out.output_path = cell("output_folder");

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

        QString title = QString("%1 (%2)").arg(label, QString::fromStdString(cfg.original_file_extention));
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
    std::string extention = cfg.original_file_extention;
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