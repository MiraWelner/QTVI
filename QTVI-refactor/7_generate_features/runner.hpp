#pragma once

#include "generate_features.hpp"
#include "bin_io.hpp"
#include <string>
#include <vector>
#include <set>
#include <fstream>
#include <sstream>
#include <iostream>
#include <chrono>
#include <filesystem>

/**
 * @file runner.hpp
 * @brief CLI entry point for the PPG feature pipeline.
 *
 * Reads config.csv, prompts for dataset selection, discovers .bin files,
 * runs multithreaded feature extraction, writes output .bin.
 */

namespace ppg {

    namespace fs = std::filesystem;
    struct DatasetConfig {
        std::string data_type;              ///< col 0:  MESA / Bittium / CHAOS
        std::string annealed_path;          ///< col 6:  annealedPath
        std::string r_mark_path;            ///< col 7:  wavePath (R-peak / wave markings)
        std::string template_path;          ///< col 8:  templatePath (step 5 output)
        std::string template_marking_path;  ///< col 9:  marking_path (step 6 output)
        std::string output_path;
    };

    /**
     * @brief Parse config.csv and return one DatasetConfig per row.
     *

     * @param path  Path to config.csv.
     * @return Vector of DatasetConfig.
     */
    inline std::vector<DatasetConfig> read_config_csv(const std::string& path) {
        std::ifstream file(path);
        if (!file) throw std::runtime_error("Cannot open config: " + path);

        std::vector<DatasetConfig> configs;
        std::string line;
        std::getline(file, line); // skip header

        while (std::getline(file, line)) {
            if (line.empty()) continue;

            std::stringstream ss(line);
            std::vector<std::string> cols;
            std::string field;
            while (std::getline(ss, field, ',')) {
                auto lo = field.find_first_not_of(" \t\r\n");
                auto hi = field.find_last_not_of(" \t\r\n");
                cols.push_back((lo == std::string::npos) ? "" : field.substr(lo, hi - lo + 1));
            }

            if (cols.size() < 11) continue;

            configs.push_back({
                cols[0],   // data_type
                cols[6],   // annealed_path
                cols[7],   // r_mark_path
                cols[8],   // template_path
                cols[9],   // template_marking_path
                cols[10],  // output_path
                });
        }
        return configs;
    }

    // ─── File discovery ─────────────────────────────────────────────────────────

    /**
     * @brief Scan directory for .bin files and extract subject IDs.
     * @param directory  Path to scan.
     * @param suffix     Suffix before .bin (e.g. "_wave_markings"). Empty = stem is the ID.
     * @return Vector of subject ID strings.
     */
    inline std::vector<std::string> parse_ids(
        const std::string& directory,
        const std::string& suffix)
    {
        std::vector<std::string> ids;
        if (!fs::exists(directory)) return ids;

        for (auto& entry : fs::directory_iterator(directory)) {
            if (!entry.is_regular_file()) continue;
            std::string ext = entry.path().extension().string();
            if (ext != ".bin" && ext != ".BIN") continue;

            std::string stem = entry.path().stem().string();
            if (stem.find("_noise_markings") != std::string::npos) continue;

            if (suffix.empty()) {
                ids.push_back(stem);
            }
            else {
                auto pos = stem.find(suffix);
                if (pos != std::string::npos)
                    ids.push_back(stem.substr(0, pos));
            }
        }
        return ids;
    }

    /**
     * @brief Intersect multiple ID lists, preserving order of the first.
     */
    inline std::vector<std::string> intersect_ids(
        const std::vector<std::vector<std::string>>& lists)
    {
        if (lists.empty()) return {};
        std::vector<std::string> result = lists[0];
        for (size_t i = 1; i < lists.size(); ++i) {
            std::set<std::string> s(lists[i].begin(), lists[i].end());
            std::vector<std::string> tmp;
            for (auto& id : result)
                if (s.count(id)) tmp.push_back(id);
            result = std::move(tmp);
        }
        return result;
    }

    // ─── Pipeline ───────────────────────────────────────────────────────────────

    /**
     * @brief Run the feature-extraction pipeline.
     * @param config_csv_path  Path to config.csv.
     * @param skip_existing    Skip subjects whose output .bin already exists.
     */
    inline void run_pipeline(
        const std::string& config_csv_path,
        bool skip_existing = true)
    {
        auto configs = read_config_csv(config_csv_path);
        if (configs.empty()) throw std::runtime_error("No datasets in config.csv");

        // Dataset selection
        std::cout << "Select Dataset:\n";
        for (size_t i = 0; i < configs.size(); ++i)
            std::cout << "  " << (i + 1) << ") " << configs[i].data_type << "\n";

        int choice = 0;
        while (choice < 1 || choice > static_cast<int>(configs.size())) {
            std::cout << "> ";
            if (!(std::cin >> choice)) { std::cin.clear(); std::cin.ignore(10000, '\n'); choice = 0; }
        }
        const auto& cfg = configs[choice - 1];

        std::cout << "\nAnnealed:   " << cfg.annealed_path
            << "\nWave:       " << cfg.r_mark_path
            << "\nTemplates:  " << cfg.template_path
            << "\nMarkings:   " << cfg.template_marking_path
            << "\nOutput:     " << cfg.output_path << "\n\n";

        if (!fs::exists(cfg.output_path))
            fs::create_directories(cfg.output_path);

        // Discover subjects present in all input directories
        auto anneal_ids = parse_ids(cfg.annealed_path, "");
        auto wave_ids = parse_ids(cfg.r_mark_path, "_wave_markings");

        bool has_templates = !cfg.template_path.empty() && fs::exists(cfg.template_path);
        std::vector<std::string> subject_ids;

        if (has_templates) {
            auto tmpl_ids = parse_ids(cfg.template_path, "_template_info");
            subject_ids = intersect_ids({ anneal_ids, wave_ids, tmpl_ids });
        }
        else {
            subject_ids = intersect_ids({ anneal_ids, wave_ids });
        }

        bool has_markings = !cfg.template_marking_path.empty()
            && fs::exists(cfg.template_marking_path);

        std::cout << "Matched " << subject_ids.size() << " subjects.\n";
        std::cout << std::string(70, '*') << "\n\n";

        double total_time = 0.0;
        int success = 0, fail = 0;

        for (size_t i = 0; i < subject_ids.size(); ++i) {
            const auto& id = subject_ids[i];

            std::string out_file = (fs::path(cfg.output_path) /
                (id + "_feature_output.bin")).string();
            if (skip_existing && fs::exists(out_file)) {
                std::cout << id << "_feature_output.bin exists — skipping.\n";
                continue;
            }

            auto t0 = std::chrono::steady_clock::now();
            double avg = (i > 0) ? total_time / i : 0.0;
            double est = avg * (subject_ids.size() - i) / 60.0;

            std::cout << "Beginning " << id
                << " | " << (i + 1) << "/" << subject_ids.size()
                << " | Avg " << avg << "s"
                << " | Est " << est << " min remaining\n";

            try {
                // Load annealed segments
                auto segments = io::read_annealed_bin(
                    (fs::path(cfg.annealed_path) / (id + ".bin")).string());

                // Load wave markings (R-peaks + PPG pairs)
                auto wave_data = io::read_wave_data_bin(
                    (fs::path(cfg.r_mark_path) / (id + "_wave_markings.bin")).string());

                // Load templates (step 5) + markings overlay (step 6)
                std::vector<TemplateInfo> templates;
                if (has_templates) {
                    templates = io::read_template_info_bin(
                        (fs::path(cfg.template_path) / (id + "_template_info.bin")).string());

                    if (has_markings) {
                        io::read_template_markings_bin(
                            (fs::path(cfg.template_marking_path) /
                                (id + "_template_markings.bin")).string(),
                            templates);
                    }
                }

                // Run (multithreaded per bin)
                auto result = generate_features(
                    segments, wave_data, templates,
                    /*window_length=*/ 30,
                    /*sqi_threshold=*/ kInf,
                    /*num_threads=*/   0);

                // Save
                std::cout << "Saving " << result.beats.size() << " beats...\n";
                io::write_feature_output_bin(out_file, result);
                ++success;

            }
            catch (const std::exception& e) {
                std::cerr << "ERROR " << id << ": " << e.what() << "\n";
                ++fail;
            }

            double elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t0).count();
            total_time += elapsed;
            std::cout << "Done in " << elapsed << "s\n"
                << std::string(100, '_') << "\n\n";
        }

        std::cout << "Complete. Success: " << success << " | Failed: " << fail << "\n";
    }

} // namespace ppg