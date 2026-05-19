/**
 * @file   input_file_handler.cpp
 * @brief  Loads the merged config.csv, walks the source-file directory,
 *         and converts source files to .bin via file_to_bin. The marking
 *         GUI's main() pulls one .bin at a time from make_binfile and
 *         doesn't need to know whether it was just produced or cached.
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
}   // anonymous namespace


bool load_config(int dataType, config_entry& out) {
    /*
    * Loads the config.csv (or whatever is in CONFIG_PATH) and fills up a config_entry 
    * struct with all the data. The config.csv should only have an output_folder param while the 
    * config_entry has all the paths to the subfolders found within output_folder, like noise_file_path, so those
    * are created within this function
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
        auto row = parseCsvRow(line);
        std::string rowType = row[0];
        std::transform(rowType.begin(), rowType.end(), rowType.begin(), ::toupper);
        if (rowType != input_file_type) continue;
        out.dataType = rowType;
        out.mainExt = row[1];
        out.sleepExt = row[2];
        out.ecgRate = std::stod(row[3]);
        out.ppgRate = std::stod(row[4]);
        out.finalSamplingRate = std::stod(row[5]);
        out.bin_length_minutes = std::stod(row[6]);
        out.input_path = row[7];
        out.output_path = row[8];
        out.bin_file_path = out.output_path + "/input_binfile/";
        out.noise_data_path = out.output_path + "/noise_marking_output/";
        out.annealed_data_path = out.output_path + "/annealed_output/";
        out.r_peak_data_path = out.output_path + "/r_peak_finding_output/";
        out.template_path = out.output_path + "/template_path/";

        apply_dataset_specific_channel_labels(out);
        return true;
    }
    return false;
}

bool promptForMissingPaths(config_entry& cfg) {
    // Ordered so the user is asked for inputs first (sources, then where to
    // cache derived data, then output trees). Each entry is (label shown in
    // the dialog title, pointer to the field to fill).
    const std::vector<std::pair<const char*, std::string*>> fields = {
        { "Original .EDF or .Dat files:",    &cfg.input_path },
        { "Output Folder",           &cfg.output_path      },
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
    return true;
}

QStringList discoverSourceFiles(const config_entry& cfg) {
    std::string srcExt = cfg.mainExt;
    std::transform(srcExt.begin(), srcExt.end(), srcExt.begin(), ::toupper);
    QString srcExtQ = QString::fromStdString(srcExt);

    QDirIterator it(QString::fromStdString(cfg.input_path),
        { "*" + srcExtQ.toLower(), "*" + srcExtQ },
        QDir::Files, QDirIterator::Subdirectories);

    QStringList out;
    while (it.hasNext()) out << it.next();
    out.sort();
    return out;
}