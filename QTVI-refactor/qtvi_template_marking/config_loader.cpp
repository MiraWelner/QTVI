/**
 * @file   config_loader.cpp
 * @brief  Loads the rate / bin-length cells for the user-selected dataset
 *         from config.csv and prompts for the noise-marking GUI's output
 *         folder. QTVI_markings/ is created inside that folder as the
 *         destination for this GUI's *_template_markings.bin output.
 */

#include "config_loader.hpp"

#include <QFileDialog>
#include <QString>

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

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
            size_t first = f.find_first_not_of(" \t\r\n");
            if (first == std::string::npos) { f.clear(); continue; }
            size_t last = f.find_last_not_of(" \t\r\n");
            f = f.substr(first, last - first + 1);
        }
        return fields;
    }

    double stod_or_zero(const std::string& s) {
        if (s.empty()) return 0.0;
        try { return std::stod(s); }
        catch (...) { return 0.0; }
    }

    // Inside the user-chosen root, this GUI's deliverables live in a single
    // subfolder. template_path is wherever the prior pipeline step
    // (template generation) wrote its templates; we read from there.
    void derive_paths(config_entry& cfg) {
        cfg.template_path = cfg.output_path + "/template_path/";
        cfg.qtvi_marker_path = cfg.output_path + "/QTVI_markings/";
        std::filesystem::create_directories(cfg.qtvi_marker_path);
    }

    // Prompt the user once for the noise-marking GUI's output folder.
    bool prompt_for_output_folder(config_entry& cfg, const std::string& dataType) {
        QString title = QString("Noise-marking GUI output folder (%1)")
            .arg(QString::fromStdString(dataType));
        QString chosen = QFileDialog::getExistingDirectory(
            nullptr, title, QString(),
            QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
        if (chosen.isEmpty()) return false;
        cfg.output_path = chosen.toStdString();
        derive_paths(cfg);
        return true;
    }

}  // anonymous namespace

std::optional<config_entry> load_config(int dataType) {
    std::ifstream file(CONFIG_PATH);
    if (!file.is_open()) {
        std::cerr << "ERROR: cannot open " << CONFIG_PATH << "\n";
        return std::nullopt;
    }

    const std::string wanted_type =
        (dataType == 1) ? "MESA" :
        (dataType == 2) ? "BITTIUM" :
        (dataType == 3) ? "CHAOS" : "";

    std::string line;
    std::getline(file, line);   // skip header

    while (std::getline(file, line)) {
        auto row = parse_csv_row(line);
        // We only need columns 0 (type), 8 (final rate), and 9 (bin minutes).
        // Other columns belong to sibling programs in the pipeline and are
        // ignored here.
        if (row.size() < 10) continue;

        std::string rowType = row[0];
        std::transform(rowType.begin(), rowType.end(), rowType.begin(), ::toupper);
        if (rowType != wanted_type) continue;

        config_entry cfg;
        cfg.finalSamplingRate = stod_or_zero(row[8]);
        cfg.bin_length_minutes = stod_or_zero(row[9]);

        if (!prompt_for_output_folder(cfg, rowType)) return std::nullopt;
        return cfg;
    }
    return std::nullopt;
}