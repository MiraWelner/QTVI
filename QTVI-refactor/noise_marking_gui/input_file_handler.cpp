/**
 * @file   input_file_handler.cpp
 * @brief  Loads the merged config.csv, walks the source-file directory,
 *         and converts source files to .bin via file_to_bin. The marking
 *         GUI's main() pulls one .bin at a time from convertToBin and
 *         doesn't need to know whether it was just produced or cached.
 */

#include "input_file_handler.hpp"
#include "file_to_bin.hpp"

#include <QDir>
#include <QDirIterator>
#include <QFileInfo>

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

    // Channel labels per dataset. These were dropped from config.csv during
    // the merge; if they ever need to vary per recording, restore those
    // columns to config.csv and parse them in loadConfig instead.
    //
    // REPLACE these placeholder strings with the actual channel-name
    // substrings from your data. The find() in build_edf_channel_map does
    // case-insensitive substring matching against EDF channel labels.
    void applyDefaultChannelLabels(config_entry& cfg) {
        if (cfg.dataType == "MESA") {
            cfg.ecg1Label = "EKG";
            cfg.ppgLabel = "Pleth";
        }
        else if (cfg.dataType == "Bittium") {
            cfg.ecg1Label = "ECG";
        }
        else if (cfg.dataType == "CHAOS") {
            cfg.ecg1Label = "NLS_NOM_ECG_ELEC_pOTL_I";
            cfg.ecg2Label = "NLS_NOM_ECG_ELEC_pOTL_II";
            cfg.ecg3Label = "NLS_NOM_ECG_ELEC_pOTL_V";
            cfg.ppgLabel = "NLS_NOM_PLETH";
        }
    }

    // Find a sleep-stage XML next to the source file when sleepExt is set.
    std::filesystem::path findSleepXml(const std::filesystem::path& src,
        const std::string& sleepExt) {
        if (sleepExt.empty()) return {};
        std::string want = sleepExt;
        std::transform(want.begin(), want.end(), want.begin(), ::toupper);
        std::string stem = src.stem().string();
        for (const auto& f :
            std::filesystem::directory_iterator(src.parent_path())) {
            std::string e = f.path().extension().string();
            std::transform(e.begin(), e.end(), e.begin(), ::toupper);
            if (e == want && f.path().stem().string().find(stem) != std::string::npos)
                return f.path();
        }
        return {};
    }

}   // anonymous namespace

bool loadConfig(int dataType, config_entry& out) {
    std::ifstream file(CONFIG_PATH);
    if (!file.is_open()) {
        std::cerr << "ERROR: cannot open " << CONFIG_PATH << "\n";
        return false;
    }

    std::string target = (dataType == 1) ? "MESA"
        : (dataType == 2) ? "BITTIUM"
        : (dataType == 3) ? "CHAOS" : "";

    std::string line;
    std::getline(file, line);   // skip header

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        auto row = parseCsvRow(line);
        if (row.size() < 10) continue;

        std::string rowType = row[0];
        std::transform(rowType.begin(), rowType.end(), rowType.begin(), ::toupper);
        if (rowType != target) continue;

        out.dataType = row[0];
        out.mainExt = row[1];
        out.sleepExt = row[2];
        try {
            out.ecgRate = !row[3].empty() ? std::stod(row[3]) : 0.0;
            out.ppgRate = !row[4].empty() ? std::stod(row[4]) : 0.0;
            out.finalSamplingRate = !row[5].empty() ? std::stod(row[5]) : 1000.0;
        }
        catch (...) {
            out.finalSamplingRate = 1000.0;
        }
        out.bin_length_minutes = std::stod(row[6]);

        out.original_file_path = row[7];
        out.bin_file_path = row[8];
        out.noise_data_path = row[9];
        out.annealed_data_path = row[10];
        out.r_peak_data_path = row[11];
        out.template_path = row[12];

        applyDefaultChannelLabels(out);
        return true;
    }
    return false;
}

QStringList discoverSourceFiles(const config_entry& cfg) {
    std::string srcExt = cfg.mainExt;
    std::transform(srcExt.begin(), srcExt.end(), srcExt.begin(), ::toupper);
    QString srcExtQ = QString::fromStdString(srcExt);

    QDirIterator it(QString::fromStdString(cfg.original_file_path),
        { "*" + srcExtQ.toLower(), "*" + srcExtQ },
        QDir::Files, QDirIterator::Subdirectories);

    QStringList out;
    while (it.hasNext()) out << it.next();
    out.sort();
    return out;
}

std::filesystem::path convertToBin(const std::filesystem::path& src,
    const config_entry& cfg) {
    std::filesystem::path out =
        std::filesystem::path(cfg.bin_file_path) /
        (src.stem().string() + "_" +
            std::to_string((int)cfg.finalSamplingRate) + "_" + std::format("{:03d}", static_cast<int>(cfg.bin_length_minutes)) + ".bin");

    if (std::filesystem::exists(out)) {
        std::cout << "  Using existing bin file: " << out.filename().string() << "\n";
        return out;
    }

    std::filesystem::create_directories(cfg.bin_file_path);

    std::string ext = src.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), ::toupper);
    if (ext == ".EDF") {
        make_binfile_edf(src, findSleepXml(src, cfg.sleepExt), cfg);
    }
    else if (ext == ".DAT" || ext == ".CSV") {
        make_binfile_dat(src, cfg);
    }
    else {
        std::cerr << "ERROR: unsupported extension " << ext
            << " for " << src << "\n";
        return {};
    }
    return out;
}