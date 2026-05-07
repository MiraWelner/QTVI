/**
 * @file   post_process.hpp
 * @brief  Per-file pipeline: anneal each input bin, then immediately
 *         compute R-peak / wave-markings on the annealed output.
 *         Header-only because it's only used from gui_handler.cpp and
 *         noise_marking_gui.cpp.
 */

#pragma once

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "annealing_to_bin/anneal_handler.hpp"
#include "config_entry.hpp"
#include "peak_finding/binfile_handling.hpp"
#include "peak_finding/create_ecg_ppg_pairs.hpp"

namespace post_process_detail {

    // Runs the full pipeline for a single file. Returns true on success.
    inline bool processOneFile(const config_entry& cfg,
        const std::filesystem::path& binPath)
    {
        const std::string stem = binPath.stem().string();
        const std::filesystem::path noisePath =
            std::filesystem::path(cfg.noise_data_path) / (stem + "_noise_markings.bin");
        const std::filesystem::path annealedPath =
            std::filesystem::path(cfg.annealed_data_path) / (stem + ".bin");
        const std::filesystem::path rPeakPath =
            std::filesystem::path(cfg.r_peak_data_path) / (stem + "_wave_markings.bin");

        // ---- Step 1: Anneal (with cache) ----
        bool annealedFresh = std::filesystem::exists(annealedPath) &&
            std::filesystem::last_write_time(annealedPath) >=
            std::filesystem::last_write_time(binPath);
        if (annealedFresh && std::filesystem::exists(noisePath)) {
            annealedFresh = std::filesystem::last_write_time(annealedPath) >=
                std::filesystem::last_write_time(noisePath);
        }

        if (annealedFresh) {
            std::cerr << "  [cache] anneal: " << annealedPath.filename() << " up to date\n";
        }
        else {
            std::cerr << "  Annealing: " << stem << "\n";
            if (!annealOneFile(binPath, noisePath, annealedPath, cfg.bin_length_minutes)) {
                std::cerr << "  Skipping R-peaks for " << stem << "\n";
                return false;
            }
        }

        // ---- Step 2: R-peaks / wave-markings (with cache) ----
        bool waveFresh = std::filesystem::exists(rPeakPath) &&
            std::filesystem::last_write_time(rPeakPath) >=
            std::filesystem::last_write_time(annealedPath);

        if (waveFresh) {
            std::cerr << "  [cache] r-peaks: " << rPeakPath.filename() << " up to date\n";
        }
        else {
            std::cerr << "  R-peaks: " << stem << "\n";
            try {
                AnnealedData annealedData = read_input_binfile(annealedPath.string());
                std::vector<output_binfile_data> results = create_ecg_ppg_pairs(
                    std::move(annealedData.bins), 0, true, stem);
                write_output_binfile(rPeakPath.string(), results);
                std::cerr << "  -> " << results.size() << " bins -> "
                    << rPeakPath.filename() << "\n";
            }
            catch (const std::exception& e) {
                std::cerr << "  ERROR (r-peaks) " << stem << ": " << e.what() << "\n";
                return false;
            }
        }
        return true;
    }

}  // namespace post_process_detail

inline int processDataset(const config_entry& cfg) {
    if (cfg.annealed_data_path.empty()) {
        std::cerr << "ERROR: annealed_data_path is empty in config\n";
        return 0;
    }
    if (cfg.r_peak_data_path.empty()) {
        std::cerr << "ERROR: r_peak_data_path is empty in config\n";
        return 0;
    }

    std::filesystem::create_directories(cfg.annealed_data_path);
    std::filesystem::create_directories(cfg.r_peak_data_path);

    int processed = 0;
    int total = 0;
    for (const auto& entry :
        std::filesystem::recursive_directory_iterator(cfg.bin_file_path))
    {
        if (!entry.is_regular_file()) continue;

        const std::string name = entry.path().filename().string();
        if (name.find("_noise_markings") != std::string::npos) continue;

        std::string ext = entry.path().extension().string();
        std::transform(ext.begin(), ext.end(), ext.begin(), ::toupper);
        if (ext != ".BIN") continue;

        ++total;
        std::cerr << "=== " << entry.path().stem().string() << " ===\n";
        if (post_process_detail::processOneFile(cfg, entry.path())) ++processed;
    }

    if (total == 0) {
        std::cerr << "No .bin files found in " << cfg.bin_file_path << "\n";
    }
    else {
        std::cerr << "Done. Processed " << processed << " of " << total << " files.\n";
    }
    return processed;
}