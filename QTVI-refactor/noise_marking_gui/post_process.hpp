#pragma once

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "annealing_to_bin/anneal_handler.hpp"
#include "peak_finding/peakfinding_io.hpp"
#include "config_entry.hpp"
#include "peak_finding/create_ecg_ppg_pairs.hpp"
#include "template_generation/build_templates.hpp"
#include "peak_finding/run_find_r_peaks.hpp"
#include "template_generation/template_io.hpp"

namespace post_process_detail {

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
        const std::filesystem::path templatePath =
            std::filesystem::path(cfg.template_path) /
            (stem + "_" + std::to_string(static_cast<int>(cfg.bin_length_minutes)) +
                "_templates.bin");

        // ---- Step 1: Anneal ----
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
                std::cerr << "  Skipping rest of pipeline for " << stem << "\n";
                return false;
            }
        }

        // ---- Step 2: R-peaks / wave-markings ----
        bool waveFresh = std::filesystem::exists(rPeakPath) &&
            std::filesystem::last_write_time(rPeakPath) >=
            std::filesystem::last_write_time(annealedPath);

        std::vector<output_binfile_data> peakResults;
        bool peakResultsInMemory = false;

        if (waveFresh) {
            std::cerr << "  [cache] r-peaks: " << rPeakPath.filename() << " up to date\n";
        }
        else {
            std::cerr << "  R-peaks: " << stem << "\n";
            try {
                AnnealedData annealedData = read_input_binfile(annealedPath.string());
                peakResults = create_ecg_ppg_pairs(
                    std::move(annealedData.bins), 0, true, stem, cfg.finalSamplingRate, cfg.finalSamplingRate);
                write_output_binfile(rPeakPath.string(), peakResults);
                std::cerr << "  -> " << peakResults.size() << " bins -> "
                    << rPeakPath.filename() << "\n";
                peakResultsInMemory = true;
            }
            catch (const std::exception& e) {
                std::cerr << "  ERROR (r-peaks) " << stem << ": " << e.what() << "\n";
                return false;
            }
        }

        // ---- Step 3: Templates + SAECG + beats ----
        bool templatesFresh = std::filesystem::exists(templatePath) &&
            std::filesystem::last_write_time(templatePath) >=
            std::filesystem::last_write_time(rPeakPath);

        if (templatesFresh) {
            std::cerr << "  [cache] templates: " << templatePath.filename() << " up to date\n";
        }
        else {
            std::cerr << "  Templates: " << stem << "\n";
            try {
                if (!peakResultsInMemory) {
                    // wave_markings holds the R-peak indices, PPG event
                    // indices, preprocessed (squared/absval) ECG channels,
                    // noise flags and pairs. The raw ECG/PPG signals and
                    // the bin-index ranges come from the annealed .bin,
                    // so pass both paths into the re-hydrating overload.
                    std::cerr << "    loading r-peak data from " << rPeakPath.filename()
                        << " (+ raw signals from " << annealedPath.filename() << ")\n";
                    peakResults = read_output_binfile(rPeakPath.string(),
                        annealedPath.string());
                }

                auto [templateData, beatsData] =
                    buildTemplatesAndBeatsFromPeakResults(peakResults);

                std::filesystem::create_directories(cfg.template_path);

                template_io::write_template_binfile(templatePath.string(), templateData);
                std::cerr << "    wrote " << templateData.bins.size()
                    << " bin templates + SAECG -> "
                    << templatePath.filename() << "\n";

                std::filesystem::path beatsPath =
                    std::filesystem::path(cfg.template_path) / (stem + "_beats.bin");
                template_io::write_beats_binfile(beatsPath.string(), beatsData);
                std::cerr << "    wrote ch1-raw beats -> "
                    << beatsPath.filename() << "\n";
            }
            catch (const std::exception& e) {
                std::cerr << "  ERROR (templates) " << stem << ": " << e.what() << "\n";
                return false;
            }
        }
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
    if (cfg.template_path.empty()) {
        std::cerr << "ERROR: template_path is empty in config\n";
        return 0;
    }

    std::filesystem::create_directories(cfg.annealed_data_path);
    std::filesystem::create_directories(cfg.r_peak_data_path);
    std::filesystem::create_directories(cfg.template_path);

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