#pragma once

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <omp.h>

#include "annealing/anneal_handler.hpp"
#include "peak_finding/peakfinding_io.hpp"
#include "config_entry.hpp"
#include "peak_finding/create_ecg_ppg_pairs.hpp"
#include "template_generation/build_templates.hpp"
#include "peak_finding/run_find_r_peaks.hpp"
#include "template_generation/template_io.hpp"
#include "template_marking_gui/alignment.hpp"   // find_q_column (Q-align)
#include "quality_check/sqi_ecg.hpp"   // computeEcgSQI / writeEcgSQICsv -> cfg.quality_metric

namespace post_process_detail {

    // Second-pass anchor sequence, one applied per "Finish and Next".
    // R is the primary build (already in job.tmpl), so it is NOT in this list;
    // the cycle is R (shown first) -> each of these in order.
    inline const std::vector<AnchorType>& anchorSequence() {
        static const std::vector<AnchorType> seq = {
            AnchorType::Q_ONSET, AnchorType::J_POINT, AnchorType::T_PEAK,
            AnchorType::P_PEAK,  AnchorType::P_ONSET
        };
        return seq;
    }

    // Short label for the corner readout / CSV tags.
    inline const char* anchorName(AnchorType a) {
        switch (a) {
        case AnchorType::P_ONSET: return "P_ONSET";
        case AnchorType::P_PEAK:  return "P_PEAK";
        case AnchorType::Q_ONSET: return "Q_ONSET";
        case AnchorType::R_PEAK:  return "R_PEAK";
        case AnchorType::J_POINT: return "J_POINT";
        case AnchorType::T_PEAK:  return "T_PEAK";
        }
        return "?";
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
    // are fully populated -- so the freshness checks in prepareViewerJob
    // (annealed / wave / templates up-to-date) stay correct even if the app
    // is killed mid-marking (canonical files just won't exist yet and
    // regenerate next run). The worker touches only the canonical paths; the
    // viewer reads only the provisional path; the two never collide.
    //
    // Callers MUST join the finalize worker before invoking
    // regenerateWithAnchor -- both mutate job.tmpl and would race if run
    // concurrently. (Also documented at each function's declaration below.)
    // ---------------------------------------------------------------------
    struct ViewerJob {
        std::filesystem::path viewerTemplatePath;   // what the viewer opens
        bool needsFinalize = false;                 // false => everything already cached

        // Anchor cycle cursor. 0 = next reload builds anchorSequence()[0].
        // Incremented by the controller after each successful reload.
        size_t anchorStep = 0;

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

        // Carried into the finalize step so it can rerun the squared/absval
        // R-peak detection on the worker thread. Cheaper than re-plumbing
        // the whole cfg into the worker; these are the only fields
        // augment_ecg_ppg_pairs_sqabs actually reads.
        config_entry cfg{};
        bool use_R_algorithm = true;    // = cfg.use_consensus_rpeak in prepareViewerJob
        bool ecg1_inverted = false;
        bool ecg2_inverted = false;
        bool ecg3_inverted = false;

        std::string error;                          // set by finalizeViewerJob on failure
    };

    // All analysis CSVs land in ONE shared folder, a sibling of the template
    // output directory: <template_path>/../csv_for_analysis. Change this one
    // function (or point it at a config field) to relocate every CSV at once.
    inline std::filesystem::path analysisCsvDir(const std::filesystem::path& templatePath) {
        return templatePath.parent_path().parent_path() / "csv_for_analysis";
    }

    inline std::optional<ViewerJob> prepareViewerJob(const config_entry& cfg, const std::filesystem::path& binPath,
        bool ecg1_inverted, bool ecg2_inverted, bool ecg3_inverted)
    {
        const std::string stem = binPath.stem().string();
        const std::filesystem::path noisePath = std::filesystem::path(cfg.noise_data_path) / (stem + "_noise_markings.bin");
        const std::filesystem::path annealedPath = std::filesystem::path(cfg.annealed_data_path) / (stem + "_annealed.bin");
        const std::filesystem::path rPeakPath = std::filesystem::path(cfg.r_peak_data_path) / (stem + "_peak_locations_all_beats.bin");
        const std::filesystem::path templatePath = std::filesystem::path(cfg.template_path) / (stem + "_templates.bin");
        const std::filesystem::path provisionalPath = std::filesystem::path(cfg.template_path) / (stem + "_templates.partial.bin");

        // ---- Step 1: Anneal (always -- freshness check removed) ----
        annealOneFile(binPath, noisePath, annealedPath, cfg.bin_size_minutes, ecg1_inverted, ecg2_inverted, ecg3_inverted);

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
        // Everything finalize needs from cfg. Copy it once here rather than
        // wiring individual fields piecewise later.
        job.cfg = cfg;
        job.use_R_algorithm = cfg.use_consensus_rpeak;
        job.ecg1_inverted = ecg1_inverted;
        job.ecg2_inverted = ecg2_inverted;
        job.ecg3_inverted = ecg3_inverted;
        job.rPeakPath = rPeakPath;
        job.templatePath = templatePath;
        job.provisionalPath = provisionalPath;
        job.annealedPath = annealedPath;

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

        auto t0 = std::chrono::steady_clock::now();
        job.peakResults = create_ecg_ppg_pairs_raw(std::move(annealedData.bins), true, stem, cfg,
            annealedData.ecg1_inverted, annealedData.ecg2_inverted, annealedData.ecg3_inverted);
        auto t1 = std::chrono::steady_clock::now();
        std::cerr << "  [timing] create_ecg_ppg_pairs_raw: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
            << " ms  (use_consensus_rpeak=" << cfg.use_consensus_rpeak << ")\n";
        job.needSqabsDetection = true;

        // create_ecg_ppg_pairs_raw doesn't carry the arterial pass-through
        // channels, so attach them here (parallel by bin index).
        for (size_t i = 0; i < job.peakResults.size() && i < nAnnealed; ++i) {
            job.peakResults[i].abpSignal = std::move(abpSlots[i]);
            job.peakResults[i].artSignal = std::move(artSlots[i]);
            job.peakResults[i].artPulmSignal = std::move(artpSlots[i]);
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

    // Runs on a worker thread. Must not touch Qt. Stores any error in
    // job.error rather than throwing across the thread boundary.
    //
    // Callers MUST join this worker before invoking regenerateWithAnchor:
    // both operate on job.tmpl and would race if allowed to run
    // concurrently. That discipline is already documented at
    // regenerateWithAnchor's declaration.
    inline void finalizeViewerJob(ViewerJob& job)
    {
        // Cap OpenMP so the Qt UI thread always has at least one core to
        // schedule on. Without this, augment_ecg_ppg_pairs_sqabs and the
        // downstream template builders each grab omp_get_max_threads() cores
        // by default; on an N-core machine that saturates all N cores at
        // 100%, the Qt main thread gets starved by the scheduler, and the
        // marking UI freezes -- even though it's technically off the Qt
        // thread. Leaving one core for Qt keeps the GUI responsive.
        // omp_set_nested(0) also blocks any internal parallel-for from
        // spawning a nested team inside the outer one (the same nested-
        // parallelism pattern that made create_arterial_templates slow
        // before -- see its own omp_set_nested call).
        const int hw = std::max(1, static_cast<int>(std::thread::hardware_concurrency()));
        const int workerThreads = std::max(1, hw - 1);
        omp_set_num_threads(workerThreads);
        omp_set_nested(0);

        try {
            if (job.needSqabsDetection) {
                // Fresh raw detection: persist the canonical (raw) r-peaks + CSV
                // first, before squared/absval overwrite the beat lists.
                write_output_binfile(job.rPeakPath.string(), job.peakResults);

                const std::filesystem::path csvDir = analysisCsvDir(job.templatePath);
                std::filesystem::create_directories(csvDir);
                const std::filesystem::path rPeakCsv =
                    csvDir / (job.stem + "_peak_locations_all_beats.csv");
                write_output_csvfile(rPeakCsv.string(), job.peakResults, job.fileID, job.samplingRate);
            }

            // Squared/absval R-peak detection on ECG channels, then slow
            // templating to pack the two extra per-bin blocks (squared,
            // absval) plus their SAECG entries into job.tmpl. Runs entirely
            // on this worker; no Qt access. The worker is expected to have
            // been joined before regenerateWithAnchor is called (see
            // comment above), so no race with the anchor path here.
            auto t_sq0 = std::chrono::steady_clock::now();
            augment_ecg_ppg_pairs_sqabs(job.peakResults, job.use_R_algorithm,
                job.fileID, job.samplingRate, job.cfg,
                job.ecg1_inverted, job.ecg2_inverted, job.ecg3_inverted);
            mergeTemplatesSlow(job.peakResults, job.tmpl, job.info, job.rates);
            auto t_sq1 = std::chrono::steady_clock::now();
            std::cerr << "  [timing] squared/absval build for " << job.stem << ": "
                << std::chrono::duration_cast<std::chrono::milliseconds>(t_sq1 - t_sq0).count()
                << " ms\n";

            // Canonical templates now include squared/absval alongside
            // raw/unfilt/ppg. Write once, then drop the provisional file
            // the viewer has been reading from.
            template_io::write_template_binfile(job.templatePath.string(), job.tmpl);
            std::filesystem::remove(job.provisionalPath);

            // ---- SQI: score every kept beat against this file's own,
            //      just-finalized templates. Written alongside the other
            //      per-file outputs, into cfg.quality_metric. ----
            writeEcgSQICsv(job.cfg, job.stem, job.tmpl, job.beats, job.samplingRate);

            std::cout << "Squared/absval slow processing done for " << job.stem << "\n";

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
    // (first "Finish and Next"). Q-aligns the cached R-pass templates in place
    // (no rebuild), writes the result to a SEPARATE provisional file, and
    // points the job's viewer path at it. Sets job.tmpl to the Q-aligned copy;
    // call it after joining the finalize worker so the worker's canonical
    // (R-pass) write has already completed. Returns false if nothing could be
    // built.
    //
    inline bool regenerateWithAnchorFull(ViewerJob& job, AnchorType anchor);   // fwd decl

    inline bool regenerateWithAnchor(ViewerJob& job, AnchorType anchor)
    {
        try {
            // If the R pass was served from the on-disk cache (waveFresh &&
            // templatesFresh in prepareViewerJob), job.tmpl/job.beats were
            // never built in memory -- both empty. Aligning that would write
            // a 0-bin file ("No bins loaded" on reload), so rebuild from
            // the R-peaks instead.
            if (job.tmpl.bins.empty() || job.beats.per_channel_beats.empty())
                return regenerateWithAnchorFull(job, anchor);

            template_io::TemplateFile atmpl = job.tmpl;   // copy R-pass templates
            alignTemplatesFromCache(atmpl, job.beats, job.rates, anchor);

            const std::filesystem::path aPath =
                job.provisionalPath.parent_path() /
                (job.stem + "_templates.anchor.partial.bin");
            template_io::write_template_binfile(aPath.string(), atmpl);
            job.tmpl = std::move(atmpl);
            job.qAlignPath = aPath;
            job.viewerTemplatePath = aPath;
            return true;
        }
        catch (const std::exception& e) {
            job.error = e.what();
            return false;
        }
    }

    // Fallback for when the R-pass templates/beats weren't held in memory
    // (subject opened from the fresh on-disk cache). Reloads the R-peaks,
    // rebuilds the R-pass templates + beats, then Q-aligns exactly as the
    // primary path does.
    inline bool regenerateWithAnchorFull(ViewerJob& job, AnchorType anchor)
    {
        try {
            if (job.peakResults.empty()) {
                if (job.rPeakPath.empty() || job.annealedPath.empty()) {
                    job.error = "anchor rebuild: no cached beats and no r-peak path";
                    return false;
                }
                job.peakResults = read_output_binfile(
                    job.rPeakPath.string(), job.annealedPath.string());
            }
            if (job.peakResults.empty()) {
                job.error = "anchor rebuild: peakResults empty";
                return false;
            }

            FastTemplateBuild fast = buildTemplatesAndBeatsFast(job.peakResults, job.rates);
            if (fast.tmpl.bins.empty()) {
                job.error = "anchor rebuild: produced 0 bins";
                return false;
            }
            job.tmpl = std::move(fast.tmpl);
            job.beats = std::move(fast.beats);
            job.info = std::move(fast.info);

            template_io::TemplateFile atmpl = job.tmpl;
            alignTemplatesFromCache(atmpl, job.beats, job.rates, anchor);

            const std::filesystem::path aPath =
                job.provisionalPath.parent_path() /
                (job.stem + "_templates.anchor.partial.bin");
            template_io::write_template_binfile(aPath.string(), atmpl);
            job.tmpl = std::move(atmpl);
            job.qAlignPath = aPath;
            job.viewerTemplatePath = aPath;
            return true;
        }
        catch (const std::exception& e) {
            job.error = e.what();
            return false;
        }
    }

}  // namespace post_process_detail