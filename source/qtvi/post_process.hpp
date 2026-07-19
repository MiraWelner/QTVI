#pragma once

#include <algorithm>
#include <chrono>
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
#include "template_marking_gui/alignment.hpp"   // find_q_column (Q-align)

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
        std::filesystem::path annealedPath;   // for reloading peakResults on the Q-align pass
        std::filesystem::path qAlignPath;      // Q-aligned provisional file (2nd pass)
        bool needSqabsDetection = false;            // false when wave_markings already had them

        std::vector<output_binfile_data> peakResults;
        template_io::TemplateFile tmpl;
        template_io::BeatsFile beats;
        std::vector<TemplateInfo> info;

        std::string error;                          // set by finalizeViewerJob on failure
    };

    // Forward decl: full-rebuild fallback for the Q-align pass (defined below,
    // used by regenerateWithQAlign when cached beats/templates are absent).
    inline bool regenerateWithQAlignFull(ViewerJob& job);

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
        const std::filesystem::path annealedPath = std::filesystem::path(cfg.annealed_data_path) / (stem + "_annealed.bin");
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
        job.annealedPath = annealedPath;

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
            // Squared/absval are unused by the viewer, so we ignore them
            // entirely: no augment_ecg_ppg_pairs_sqabs (detection) and no
            // mergeTemplatesSlow (squared/absval templating). This keeps this
            // step to plain file writes -- it no longer runs the template
            // generator, so it can't race the Q-align rebuild, and the rebuild
            // never waits on heavy squared/absval work.
            if (job.needSqabsDetection) {
                // Fresh raw detection: persist the canonical (raw) r-peaks + CSV.
                write_output_binfile(job.rPeakPath.string(), job.peakResults);

                const std::filesystem::path csvDir = analysisCsvDir(job.templatePath);
                std::filesystem::create_directories(csvDir);
                const std::filesystem::path rPeakCsv =
                    csvDir / (job.stem + "_peak_locations_all_beats.csv");
                write_output_csvfile(rPeakCsv.string(), job.peakResults, job.fileID, job.samplingRate);
            }
            // Canonical templates = the fast (raw/unfilt/ppg) build as-is.
            template_io::write_template_binfile(job.templatePath.string(), job.tmpl);
            std::filesystem::remove(job.provisionalPath);

            // (No snips.csv: serializing every retained beat for every channel
            // dominated finalize time and nothing downstream consumes it.)
        }
        catch (const std::exception& e) {
            job.error = e.what();
        }
        catch (...) {
            job.error = "unknown exception in finalizeViewerJob";
        }
    }

    // Called by the controller when the viewer emits requestQAlignReload()
    // Called by the controller when the viewer emits requestQAlignReload()
    // (first "Finish and Next"). Q-aligns the cached R-pass templates in place
    // (no rebuild), writes the result to a SEPARATE provisional file, and
    // points the job's viewer path at it. Sets job.tmpl to the Q-aligned copy;
    // call it after joining the finalize worker so the worker's canonical
    // (R-pass) write has already completed. Returns false if nothing could be
    // built.
    //
    // NOTE: the Q pass no longer touches the template generator at all -- it
    // Q-shifts the cached R-aligned beats in place. So there is nothing for
    // the background finalize worker (plain file writes) to contend with;
    // joining it first is only belt-and-suspenders.
    inline bool regenerateWithQAlign(ViewerJob& job)
    {
        using clock = std::chrono::steady_clock;
        const auto t0 = clock::now();
        try {
            // Reuse the R-pass PPG/arterial templates as-is and rebuild each
            // ECG raw template by Q-aligning its cached snippets. Needs both
            // the R-pass templates and the cached beats.
            if (job.tmpl.bins.empty() || job.beats.per_channel_beats.empty()) {
                std::cerr << "  [q-align] no cached beats/templates -> full rebuild fallback\n";
                return regenerateWithQAlignFull(job);
            }

            template_io::TemplateFile qtmpl = job.tmpl;   // copy R-pass templates
            qAlignTemplatesFromCache(qtmpl, job.beats, job.rates);

            const std::filesystem::path qPath =
                job.provisionalPath.parent_path() /
                (job.stem + "_templates.qalign.partial.bin");
            template_io::write_template_binfile(qPath.string(), qtmpl);
            job.tmpl = std::move(qtmpl);
            job.qAlignPath = qPath;
            job.viewerTemplatePath = qPath;

            const double ms = std::chrono::duration<double, std::milli>(
                clock::now() - t0).count();
            std::cerr << "  [q-align] reused R-pass templates, Q-shifted ECG in "
                << ms << " ms\n";
            return true;
        }
        catch (const std::exception& e) {
            job.error = e.what();
            return false;
        }
    }

    // Full fallback: used only when the cached beats/templates are unavailable
    // (e.g. a subject opened cold on the Q pass, where prepareViewerJob short-
    // circuited on fresh canonical files and never built job.beats in memory).
    // We can't Q-align without the R-aligned beat matrices, and those are only
    // produced by the fast build -- so re-run the fast (R-only) build from the
    // R-peaks to repopulate the cache, then Q-align from it exactly as the
    // primary path does.
    inline bool regenerateWithQAlignFull(ViewerJob& job)
    {
        using clock = std::chrono::steady_clock;
        const auto t0 = clock::now();
        try {
            if (job.peakResults.empty()) {
                if (job.rPeakPath.empty() || job.annealedPath.empty()) return false;
                std::cerr << "  [q-align] fallback: reloading peakResults from disk\n";
                job.peakResults = read_output_binfile(
                    job.rPeakPath.string(), job.annealedPath.string());
            }
            if (job.peakResults.empty()) return false;

            // Rebuild the R-pass templates from the R-peaks.
            std::cerr << "  [q-align] fallback: full buildTemplatesAndBeatsFast rebuild\n";
            FastTemplateBuild fast = buildTemplatesAndBeatsFast(job.peakResults, job.rates);
            if (fast.tmpl.bins.empty()) return false;

            // Q-shift the ECG raw templates so Q lands on a common column.
            qAlignTemplatesFromCache(fast, job.rates);

            const std::filesystem::path qPath =
                job.provisionalPath.parent_path() /
                (job.stem + "_templates.qalign.partial.bin");
            template_io::write_template_binfile(qPath.string(), fast.tmpl);
            job.tmpl = std::move(fast.tmpl);
            job.beats = std::move(fast.beats);
            job.info = std::move(fast.info);
            job.qAlignPath = qPath;
            job.viewerTemplatePath = qPath;
            const double ms = std::chrono::duration<double, std::milli>(
                clock::now() - t0).count();
            std::cerr << "  [q-align] fallback rebuild done in " << ms << " ms\n";
            return true;
        }
        catch (const std::exception& e) {
            job.error = e.what();
            return false;
        }
    }

}  // namespace post_process_detail