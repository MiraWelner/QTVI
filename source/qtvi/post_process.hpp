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
    // are fully populated -- so the freshness checks in prepareViewerJob
    // (annealed / wave / templates up-to-date) stay correct even if the app
    // is killed mid-marking (canonical files just won't exist yet and
    // regenerate next run). The worker touches only the canonical paths; the
    // viewer reads only the provisional path; the two never collide.
    // ---------------------------------------------------------------------
    struct ViewerJob {
        std::filesystem::path viewerTemplatePath;   // what the viewer opens
        bool needsFinalize = false;                 // false => everything already cached

        // Carried fast -> slow (only meaningful when needsFinalize):
        std::string stem;
        std::string fileID;
        double samplingRate = 0.0;   // ECG rate (used by non-template callers below)
        SignalRates rates;           // full per-channel rate set for template pipeline
        std::filesystem::path rPeakPath, templatePath, beatsPath, provisionalPath;
        bool needSqabsDetection = false;            // false when wave_markings already had them

        std::vector<output_binfile_data> peakResults;
        template_io::TemplateFile tmpl;
        template_io::BeatsFile beats;
        std::vector<TemplateInfo> info;

        std::string error;                          // set by finalizeViewerJob on failure
    };

    // All analysis CSVs land in ONE shared folder, a sibling of the template
    // output directory: <template_path>/../csv_for_analysis. Change this one
    // function (or point it at a config field) to relocate every CSV at once.
    inline std::filesystem::path analysisCsvDir(const std::filesystem::path& templatePath) {
        return templatePath.parent_path().parent_path() / "csv_for_analysis";
    }

    inline std::optional<ViewerJob> prepareViewerJob(const config_entry& cfg, const std::filesystem::path& binPath)
    {
        const std::string stem = binPath.stem().string();
        const std::filesystem::path noisePath = std::filesystem::path(cfg.noise_data_path) / (stem + "_noise_markings.bin");
        const std::filesystem::path annealedPath = std::filesystem::path(cfg.annealed_data_path) / (stem + "annealed.bin");
        const std::filesystem::path rPeakPath = std::filesystem::path(cfg.r_peak_data_path) / (stem + "_peak_locations_all_beats.bin");
        const std::filesystem::path templatePath = std::filesystem::path(cfg.template_path) / (stem + "_templates.bin");
        const std::filesystem::path provisionalPath = std::filesystem::path(cfg.template_path) / (stem + "_templates.partial.bin");

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
        // Per-channel rates, forwarded to the template-generation pipeline.
        // A rate of 0 means the channel is absent for this dataset (the
        // slicer skips it silently).
        job.rates = SignalRates{
            cfg.ecg_upsample_rate,
            cfg.ppg_upsample_rate,
            cfg.abp_upsample_rate,
            cfg.art_upsample_rate,
            cfg.art_pulm_upsample_rate
        };
        job.rPeakPath = rPeakPath;
        job.templatePath = templatePath;
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

                // Grab arterial pass-through slots BEFORE the move consumes the bins.
                // Slots match file_to_bin / gui_handler: CH_ABP=33, CH_ART=34, CH_ART_PULM=35.
                const size_t nAnnealed = annealedData.bins.size();
                std::vector<std::vector<double>> abpSlots(nAnnealed), artSlots(nAnnealed), artpSlots(nAnnealed);
                for (size_t i = 0; i < nAnnealed; ++i) {
                    auto& up = annealedData.bins[i].all_upsampled;
                    if (33 < up.size()) abpSlots[i] = up[33];
                    if (34 < up.size()) artSlots[i] = up[34];
                    if (35 < up.size()) artpSlots[i] = up[35];
                }

                job.peakResults = create_ecg_ppg_pairs_raw(
                    std::move(annealedData.bins), true, stem,
                    cfg.ecg_upsample_rate, cfg.ppg_upsample_rate);
                job.needSqabsDetection = true;

                // create_ecg_ppg_pairs_raw doesn't carry the arterial pass-through
                // channels, so attach them here (parallel by bin index).
                for (size_t i = 0; i < job.peakResults.size() && i < nAnnealed; ++i) {
                    job.peakResults[i].abpSignal = std::move(abpSlots[i]);
                    job.peakResults[i].artSignal = std::move(artSlots[i]);
                    job.peakResults[i].artPulmSignal = std::move(artpSlots[i]);
                }
            }

            std::cerr << "  Templates (raw/unfiltered/ppg, fast): " << stem << "\n";
            FastTemplateBuild fast = buildTemplatesAndBeatsFast(job.peakResults, job.rates);
            if (fast.tmpl.bins.empty()) {
                std::cerr << "  no bins for " << stem
                    << " (recording shorter than one bin?); skipping.\n";
                return std::nullopt;   // main.cpp prints "prep failed or skipped"
            }
            job.tmpl = std::move(fast.tmpl);
            job.beats = std::move(fast.beats);
            job.info = std::move(fast.info);

            template_io::write_template_binfile(provisionalPath.string(), job.tmpl);

            job.viewerTemplatePath = provisionalPath;
            job.needsFinalize = true;
            return job;
        }
        catch (const std::exception& e) {
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

                // CSV mirror of wave_markings.bin: one row per (bin, channel,
                // method, peak) with sample index and time in seconds. Written
                // next to the .bin, same stem.
                const std::filesystem::path csvDir = analysisCsvDir(job.templatePath);
                std::filesystem::create_directories(csvDir);
                const std::filesystem::path rPeakCsv =
                    csvDir / (job.stem + "_peak_locations_all_beats.csv");
                write_output_csvfile(rPeakCsv.string(), job.peakResults, job.fileID, job.samplingRate);
            }
            mergeTemplatesSlow(job.peakResults, job.tmpl, job.info, job.rates);
            template_io::write_template_binfile(job.templatePath.string(), job.tmpl);
            std::filesystem::remove(job.provisionalPath);

            // Individual retained beats ("snips"): one column per snip.
            const std::filesystem::path csvDir2 = analysisCsvDir(job.templatePath);
            std::filesystem::create_directories(csvDir2);
            const std::filesystem::path snipsCsv =
                csvDir2 / (job.stem + "_template_snips.csv");
            template_io::write_snips_csv(snipsCsv.string(), job.beats);
        }
        catch (const std::exception& e) {
            job.error = e.what();
        }
        catch (...) {
            job.error = "unknown exception in finalizeViewerJob";
        }
    }

}  // namespace post_process_detail