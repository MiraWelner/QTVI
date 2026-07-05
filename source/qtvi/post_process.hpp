#pragma once

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "annealing/anneal_handler.hpp"
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
            (stem + "_" + std::to_string(static_cast<int>(cfg.bin_size_minutes)) +
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
            if (!annealOneFile(binPath, noisePath, annealedPath, cfg.bin_size_minutes)) {
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
                peakResults = create_ecg_ppg_pairs(std::move(annealedData.bins), 0, true, stem, cfg.ecg_upsample_rate, cfg.ppg_upsample_rate);
                write_output_binfile(rPeakPath.string(), peakResults);
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

        return true;
    }

    // ---------------------------------------------------------------------
    // GUI fast/slow split.
    //
    // prepareViewerJob() does the minimum needed to open the viewer: anneal
    // (cached), raw R-peak detection, and raw/unfiltered/PPG templates. It
    // writes a *provisional* template file (..._templates.partial.bin) for
    // the viewer to open and returns the in-memory state.
    //
    // finalizeViewerJob() does the deferred squared/absval R-peak detection
    // and templating, then writes the CANONICAL wave_markings / templates /
    // beats files. Run it on a worker thread concurrently with manual
    // template marking, then join.
    //
    // Only finalizeViewerJob writes the canonical files, and only once they
    // are fully populated -- so processOneFile's timestamp cache stays
    // correct even if the app is killed mid-marking (canonical files just
    // won't exist yet and regenerate next run). The worker touches only the
    // canonical paths; the viewer reads only the provisional path; the two
    // never collide.
    // ---------------------------------------------------------------------
    struct ViewerJob {
        std::filesystem::path viewerTemplatePath;   // what the viewer opens
        bool needsFinalize = false;                 // false => everything already cached

        // Carried fast -> slow (only meaningful when needsFinalize):
        std::string stem;
        std::string fileID;
        double samplingRate = 0.0;
        std::filesystem::path rPeakPath, templatePath, beatsPath, provisionalPath;
        bool needSqabsDetection = false;            // false when wave_markings already had them

        std::vector<output_binfile_data> peakResults;
        template_io::TemplateFile tmpl;
        template_io::BeatsFile beats;
        std::vector<TemplateInfo> info;

        std::string error;                          // set by finalizeViewerJob on failure
    };

    inline std::optional<ViewerJob> prepareViewerJob(const config_entry& cfg,
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
            (stem + "_" + std::to_string(static_cast<int>(cfg.bin_size_minutes)) +
                "_templates.bin");
        const std::filesystem::path beatsPath =
            std::filesystem::path(cfg.template_path) / (stem + "_beats.bin");
        const std::filesystem::path provisionalPath =
            std::filesystem::path(cfg.template_path) /
            (stem + "_" + std::to_string(static_cast<int>(cfg.bin_size_minutes)) +
                "_templates.partial.bin");

        // ---- Step 1: Anneal (same freshness logic as processOneFile) ----
        bool annealedFresh = std::filesystem::exists(annealedPath) &&
            std::filesystem::last_write_time(annealedPath) >=
            std::filesystem::last_write_time(binPath);
        if (annealedFresh && std::filesystem::exists(noisePath)) {
            annealedFresh = std::filesystem::last_write_time(annealedPath) >=
                std::filesystem::last_write_time(noisePath);
        }
        if (!annealedFresh) {
            std::cerr << "  Annealing: " << stem << "\n";
            if (!annealOneFile(binPath, noisePath, annealedPath, cfg.bin_size_minutes)) {
                std::cerr << "  Skipping rest of pipeline for " << stem << "\n";
                return std::nullopt;
            }
        }
        else {
            std::cerr << "  [cache] anneal: " << annealedPath.filename() << " up to date\n";
        }

        ViewerJob job;
        job.stem = stem;
        job.fileID = stem;
        job.samplingRate = cfg.ecg_upsample_rate;
        job.rPeakPath = rPeakPath;
        job.templatePath = templatePath;
        job.beatsPath = beatsPath;
        job.provisionalPath = provisionalPath;

        bool waveFresh = std::filesystem::exists(rPeakPath) &&
            std::filesystem::last_write_time(rPeakPath) >=
            std::filesystem::last_write_time(annealedPath);
        bool templatesFresh = std::filesystem::exists(templatePath) &&
            std::filesystem::last_write_time(templatePath) >=
            std::filesystem::last_write_time(rPeakPath);

        // ---- Everything already cached: open the canonical file, no worker.
        if (waveFresh && templatesFresh) {
            std::cerr << "  [cache] templates: " << templatePath.filename() << " up to date\n";
            job.viewerTemplatePath = templatePath;
            job.needsFinalize = false;
            return job;
        }

        try {
            if (waveFresh) {
                // wave_markings already holds all three methods + preprocessed
                // signals; only the (deferred) templating remains.
                std::cerr << "  [cache] r-peaks: " << rPeakPath.filename()
                    << " up to date (loading)\n";
                job.peakResults = read_output_binfile(rPeakPath.string(), annealedPath.string());
                job.needSqabsDetection = false;
            }
            else {
                std::cerr << "  R-peaks (raw, fast): " << stem << "\n";
                AnnealedData annealedData = read_input_binfile(annealedPath.string());
                job.peakResults = create_ecg_ppg_pairs_raw(
                    std::move(annealedData.bins), true, stem,
                    cfg.ecg_upsample_rate, cfg.ppg_upsample_rate);
                job.needSqabsDetection = true;
            }

            std::cerr << "  Templates (raw/unfiltered/ppg, fast): " << stem << "\n";
            FastTemplateBuild fast = buildTemplatesAndBeatsFast(job.peakResults);
            job.tmpl = std::move(fast.tmpl);
            job.beats = std::move(fast.beats);
            job.info = std::move(fast.info);

            std::filesystem::create_directories(cfg.template_path);
            template_io::write_template_binfile(provisionalPath.string(), job.tmpl);
            std::cerr << "    wrote provisional viewer file -> "
                << provisionalPath.filename() << "\n";

            job.viewerTemplatePath = provisionalPath;
            job.needsFinalize = true;
            return job;
        }
        catch (const std::exception& e) {
            std::cerr << "  ERROR (fast prep) " << stem << ": " << e.what() << "\n";
            return std::nullopt;
        }
    }

    // Runs on a worker thread. Must not touch Qt. Stores any error in
    // job.error rather than throwing across the thread boundary.
    inline void finalizeViewerJob(ViewerJob& job)
    {
        try {
            if (job.needSqabsDetection) {
                augment_ecg_ppg_pairs_sqabs(job.peakResults, true, job.fileID, job.samplingRate);
                write_output_binfile(job.rPeakPath.string(), job.peakResults);
            }

            mergeTemplatesSlow(job.peakResults, job.tmpl, job.info);

            template_io::write_template_binfile(job.templatePath.string(), job.tmpl);
            template_io::write_beats_binfile(job.beatsPath.string(), job.beats);
        }
        catch (const std::exception& e) {
            job.error = e.what();
        }
        catch (...) {
            job.error = "unknown exception in finalizeViewerJob";
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