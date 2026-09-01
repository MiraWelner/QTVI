#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <omp.h>

#include "peak_finding/channel_offset.hpp"
#include "peak_finding/run_find_r_peaks.hpp"
#include "peak_finding/create_ecg_ppg_pairs.hpp"
#include "peak_finding/peakfinding_io.hpp"

#include "template_generation/bin_archive.hpp"
#include "template_generation/template_io.hpp"
#include "template_generation/pulse_matched_filter.hpp"
#include "template_generation/build_templates.hpp"
#include "template_generation/premark_beats.hpp"

#include "annealing/anneal_handler.hpp"
#include "config_file_handling/config_entry.hpp"
#include "template_marking_gui/alignment.hpp"
#include "logging/sqi_ecg.hpp"
#include "template_morphology_grouping/morphology_csv.hpp"

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

        template_io::TemplateFile tmplR;      //  R-pass template to be reused by re-alignment

        // Per-anchor raw templates accumulated by regenerateWithAnchor, kept
        // SEPARATE from job.tmpl so the anchor cycle and the finalize worker
        // (which mutates job.tmpl to pack squared/absval) never touch the same
        // object -- no data race, so finalize no longer has to be joined before
        // the first anchor. Merged into job.tmpl.raw_anchors once, at the final
        // promote (after the worker is joined). Key = AnchorType tag.
        std::map<int, std::vector<std::array<template_io::ChannelMethodTemplate, 3>>> anchorAccum;
    };


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
        //
        // TIMED AND ANNOUNCED, because this step prints nothing of its own and
        // sits between "Saved Noise Markings" and the first [timing] line. A
        // stall here was indistinguishable from a hang: no output, no progress,
        // and the two existing instrumentation lines both live downstream of it.
        const auto _a0 = std::chrono::steady_clock::now();
        std::cerr << "[stage] anneal " << stem << " ...\n" << std::flush;
        annealOneFile(binPath, noisePath, annealedPath, cfg.bin_size_minutes, ecg1_inverted, ecg2_inverted, ecg3_inverted);
        const auto _a1 = std::chrono::steady_clock::now();
        std::cerr << "[stage] anneal " << stem << " done in "
            << std::chrono::duration<double>(_a1 - _a0).count() << " s\n" << std::flush;

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

        // ---- shift PPG, and separately ABP/ART/ART_PULM, due to ----------
        // ---- hardware lag -------------------------------------------------
        // The lag needs R peaks and foot events, both produced from a probe
        // pass -- but the shift has to land BEFORE the REAL
        // create_ecg_ppg_pairs_raw call, so SegmentPPG, ppgMinAmps,
        // R-pairing and every template derive from the shifted signal.
        // Nothing then needs renumbering and nothing can fall out of sync.
        // So: measure on a truncated COPY, shift the real segments, and let
        // the normal pass proceed.
        //
        // PPG and the arterial group (ABP/ART/ART_PULM) are measured and
        // applied INDEPENDENTLY -- they are different acquisition paths
        // with different hardware delays, not assumed to share a lag.
        // Within the arterial group, ABP is the representative channel
        // measured (it has no precomputed foot-index field the way PPG's
        // ppgMinAmps does, so its pulse locations are detected fresh here
        // via the same derivative-upstroke detector create_arterial_
        // templates.hpp uses); the resulting ONE lag is applied to all
        // three arterial channels together, since they share one
        // acquisition path.
        //
        // CHAOS ONLY: the hardware lag this measures and corrects is a
        // property of CHAOS's acquisition path specifically. MESA,
        // BITTIUM and SHHS don't have it, so measuring and shifting on
        // them would apply a correction for a lag that isn't actually
        // there. (SHHS additionally has no PPG or arterial channel to
        // measure a lag against, so there is nothing here for it to do.)
        channel_offset::set(cfg.quality_metric, stem);
        channel_offset::Result chOffPpg, chOffArt;
        const bool wantChannelOffset = (cfg.dataset_type == "CHAOS");
        if (wantChannelOffset) {
            {
                auto _cot0 = std::chrono::steady_clock::now();
                auto probeSegs = channel_offset::make_probe(
                    annealedData.bins, cfg.ecg_upsample_rate, cfg.ppg_upsample_rate);
                if (!probeSegs.empty()) {
                    auto probeResults = create_ecg_ppg_pairs_raw(
                        std::move(probeSegs), true, stem, cfg,
                        annealedData.ecg1_inverted, annealedData.ecg2_inverted,
                        annealedData.ecg3_inverted);

                    // PPG group: foot events are already on the bin
                    // (ppgMinAmps, filled in by create_ecg_ppg_pairs_raw's
                    // SegmentPPG call above).
                    std::vector<std::vector<std::size_t>> ppgFeet;
                    ppgFeet.reserve(probeResults.size());
                    for (const auto& b : probeResults) ppgFeet.push_back(b.ppgMinAmps);
                    chOffPpg = channel_offset::measure(probeResults,
                        cfg.ecg_upsample_rate, cfg.ppg_upsample_rate, ppgFeet);

                    // Arterial group: create_ecg_ppg_pairs_raw does NOT
                    // populate b.abpSignal on this in-memory probe path (it
                    // only carries all_upsampled through wholesale; abpSignal
                    // only gets rehydrated from disk by read_output_binfile's
                    // two-arg overload, which this probe deliberately
                    // bypasses to stay fast). So ABP's raw samples are read
                    // straight out of the pass-through slot (33), and its
                    // pulse locations detected fresh -- same detector
                    // create_arterial_templates.hpp's foot-anchored builder
                    // uses -- purely as a repeatable per-pulse fiducial for
                    // this correlation, not a true "foot".
                    std::vector<std::vector<std::size_t>> abpFeet;
                    abpFeet.reserve(probeResults.size());
                    if (cfg.abp_upsample_rate > 0.0) {
                        const int minSep = std::max(1,
                            static_cast<int>(std::llround(0.25 * cfg.abp_upsample_rate)));
                        for (const auto& b : probeResults) {
                            std::vector<std::size_t> feet;
                            if (b.all_upsampled.size() > 33 && !b.all_upsampled[33].empty()) {
                                const std::vector<int> locs =
                                    pulse_matched_filter::derivativePulseLocations(
                                        b.all_upsampled[33], minSep);
                                feet.assign(locs.begin(), locs.end());
                            }
                            abpFeet.push_back(std::move(feet));
                        }
                        chOffArt = channel_offset::measure(probeResults,
                            cfg.ecg_upsample_rate, cfg.abp_upsample_rate, abpFeet);
                    }
                }
                auto _cot1 = std::chrono::steady_clock::now();
                std::cerr << "  [timing] channel_offset probe: "
                    << std::chrono::duration_cast<std::chrono::milliseconds>(_cot1 - _cot0).count()
                    << " ms\n";
            }
            channel_offset::apply(annealedData.bins,
                chOffPpg, cfg.ppg_upsample_rate,
                chOffArt, cfg.abp_upsample_rate,
                cfg.art_upsample_rate, cfg.art_pulm_upsample_rate);
            channel_offset::write_log(chOffPpg, "PPG", /*append=*/false);
            channel_offset::write_log(chOffArt, "ARTERIAL", /*append=*/true);

            if (chOffPpg.ambiguous) {
                std::cerr << "  [channel_offset] " << stem
                    << ": ambiguous, PPG NOT shifted (ratio=" << chOffPpg.ratio
                    << ", lag would have been " << chOffPpg.lag_ms << " ms)\n";
            }
            else {
                std::cerr << "  [channel_offset] " << stem << ": PPG shifted "
                    << chOffPpg.lag_ms << " ms (ratio=" << chOffPpg.ratio << ")\n";
            }
            if (chOffArt.ambiguous) {
                std::cerr << "  [channel_offset] " << stem
                    << ": ambiguous, ABP/ART/ART_PULM NOT shifted (ratio=" << chOffArt.ratio
                    << ", lag would have been " << chOffArt.lag_ms << " ms)\n";
            }
            else {
                std::cerr << "  [channel_offset] " << stem << ": ABP/ART/ART_PULM shifted "
                    << chOffArt.lag_ms << " ms (ratio=" << chOffArt.ratio << ")\n";
            }
        }
        else {
            std::cerr << "  [channel_offset] " << stem << ": dataset_type="
                << cfg.dataset_type
                << " -- correction only runs for CHAOS; skipped\n";
        }

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

        std::cerr << "  Processing Raw Templates (fast stage): " << stem << "\n";
        ecg_move_log::set(cfg.quality_metric, stem);   // per-beat vertical move log
        morphology_csv::set(cfg.template_path, stem);
        FastTemplateBuild fast = buildTemplatesAndBeatsFast(job.peakResults, job.rates);
        if (fast.tmpl.bins.empty()) {
            std::cerr << "  no bins for " << stem  << " (recording shorter than one bin?); skipping.\n";
            return std::nullopt;   // main.cpp prints "prep failed or skipped"
        }
        job.tmpl = std::move(fast.tmpl);
        job.beats = std::move(fast.beats);
        job.info = std::move(fast.info);
        job.tmplR = job.tmpl;      // snapshot R frame (one copy, at prep time)

        if (!cfg.bin_archive_path.empty()) {
            const bool ok = bin_archive::writeBinFeatureArchive(
                cfg.bin_archive_path, stem, job.tmpl.bins,
                cfg.ecg_upsample_rate, "R", &job.beats);
            if (!ok)
                std::cerr << "  [bin_archive] " << stem
                << ": could not write checkpoint to " << cfg.bin_archive_path << "\n";
            else
                std::cerr << "  [bin_archive] " << stem << ": wrote checkpoint ("
                << job.tmpl.bins.size() << " bins)\n";

            // Section 5.5 length/area/volume time series, from the SAME
            // pre-deformation R-pass data. job.peakResults still holds the
            // raw per-channel ECG + detected R-peaks the proportional
            // segmenter needs (job.tmpl has only averaged templates, which
            // can't be re-segmented), so it is computed here rather than in
            // the viewer.
            const std::string ftsPath =
                cfg.bin_archive_path + "/" + stem + "_feature_timeseries.csv";
            const bool okf = normalize_features::writeFeatureTimeSeriesCsv(
                ftsPath, stem, job.peakResults, cfg.ecg_upsample_rate);
            if (!okf)
                std::cerr << "  [feature_ts] " << stem
                << ": could not write " << ftsPath << "\n";
            else
                std::cerr << "  [feature_ts] " << stem << ": wrote length/area/volume series\n";
        }
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
                auto t_io0 = std::chrono::steady_clock::now();
                write_output_binfile(job.rPeakPath.string(), job.peakResults);

                const std::filesystem::path csvDir = job.cfg.r_peak_data_path;
                std::filesystem::create_directories(csvDir);
                const std::filesystem::path rPeakCsv = csvDir / (job.stem + "_peak_locations_all_beats.csv");
                write_output_csvfile(rPeakCsv.string(), job.peakResults, job.fileID, job.samplingRate);
                auto t_io1 = std::chrono::steady_clock::now();

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
            auto t_aug = std::chrono::steady_clock::now();
            mergeTemplatesSlow(job.peakResults, job.tmpl, job.info, job.rates);
            auto t_sq1 = std::chrono::steady_clock::now();

            // Section 4.7 morphology envelope + band scores. Deliberately on
            // THIS thread, not in prepareViewerJob: prepare blocks the marking
            // UI from opening, and premark is pure compute plus two CSV writes
            // that nothing reads during marking, so it belongs on the worker
            // that already runs concurrently with the operator.
            //
            // Safe here: premark mutates neither job.beats nor job.tmpl, so it
            // cannot race the anchor path (which touches job.anchorAccum) nor
            // the mergeTemplatesSlow above (already complete). It reads the RAW
            // per-bin blocks, which augment/merge leave intact.
            //
            // Explicit dir/stem rather than premark::set(): main.cpp parks
            // finalize workers and advances to the next file, so several
            // finalizes can be in flight together and the file-scope
            // g_dir/g_stem would be a data race between them.
            {
                auto t_pm0 = std::chrono::steady_clock::now();
                premark::runAll(job.beats, job.tmpl, job.rates.ecg,
                    job.cfg.quality_metric, job.stem);
                auto t_pm1 = std::chrono::steady_clock::now();
            }

            // Base (R) templates now include squared/absval alongside
            // raw/unfilt/ppg. Write them to the PROVISIONAL file -- the
            // canonical _templates.bin is created only at the END of the
            // anchor cycle (final "Finish"), once all anchors have been
            // accumulated into it. Until then only the .partial.bin exists.
            auto t_tw0 = std::chrono::steady_clock::now();
            template_io::write_template_binfile(job.provisionalPath.string(), job.tmpl);
            auto t_tw1 = std::chrono::steady_clock::now();
            job.viewerTemplatePath = job.provisionalPath;
            std::cout << "Processing Squared and Absolute Value Templates (slow) for " << job.stem << "\n";

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
    // Called by the controller when the viewer emits requestQAlignReload()
    // (each "Finish and Next"). Aligns the R-base templates to `anchor` and
    // accumulates the result into job.anchorAccum (SEPARATE from job.tmpl, so
    // it can run concurrently with the finalize worker). Each step aligns FROM
    // the pristine R snapshot (job.tmplR / job.beatsR), never the previous
    // anchor. At the final anchor the accumulated set is folded into job.tmpl
    // and promoted to the canonical _templates.bin.
    inline bool regenerateWithAnchor(ViewerJob& job, AnchorType anchor)
    {
        try {
            // Align from the pristine R snapshot. alignTemplatesFromCache
            // reads the R base + beats (beats passed by reference, NOT copied,
            // and not mutated in the default non-scoring mode) and writes this
            // anchor's block into atmpl.raw_anchors. We then move that block
            // into the persistent job.tmpl so anchors accumulate across steps
            // into the one growing provisional file.
            //
            // NOTE: per-anchor SQI is NOT written here -- it is deferred to the
            // final anchor (below), so the interactive steps stay fast.
            auto _t0 = std::chrono::steady_clock::now();
            template_io::TemplateFile atmpl = job.tmplR;   // R base to align from (cheap vs beats)
            alignTemplatesFromCache(atmpl, job.beats, job.rates, anchor);
            auto _t1 = std::chrono::steady_clock::now();

            const int anchorTag = static_cast<int>(anchor);
            auto it = atmpl.raw_anchors.find(anchorTag);
            if (it != atmpl.raw_anchors.end())
                job.anchorAccum[anchorTag] = std::move(it->second);   // NOT job.tmpl -- avoids racing finalize

            // Per-step write is TRIMMED: the viewer only reads the CURRENT
            // anchor (readTemplateInfoBin(path, m_currentAnchor)), so the
            // provisional it opens needs just the R scalar base + this one
            // anchor's block -- NOT the whole accumulating set. This keeps the
            // per-step write constant-size instead of growing each anchor.
            // The full accumulated job.tmpl (all anchors) is written once, at
            // the final promote below.
            const auto& seq = anchorSequence();
            const bool finalAnchor =
                (!seq.empty() && anchor == seq.back());

            if (!finalAnchor) {
                // Build the viewer's provisional from the SMALL R-only base
                // (job.tmplR), not job.tmpl -- copying job.tmpl would deep-copy
                // the whole growing accumulated set every step. tmplR carries
                // the R scalar raw the viewer displays; we add just this
                // anchor's block. (Squared/absval aren't in tmplR, but the
                // viewer only draws raw, so the interactive provisional is
                // fine without them.)
                template_io::TemplateFile view = job.tmplR;   // constant-size, R only
                view.raw_anchors.clear();
                auto jt = job.anchorAccum.find(anchorTag);
                if (jt != job.anchorAccum.end())
                    view.raw_anchors[anchorTag] = jt->second;   // just this anchor
                template_io::write_template_binfile(job.provisionalPath.string(), view);
                job.viewerTemplatePath = job.provisionalPath;
                auto _t2 = std::chrono::steady_clock::now();
                std::cerr << "  [timing] anchor " << anchorName(anchor)
                    << ": align="
                    << std::chrono::duration_cast<std::chrono::milliseconds>(_t1 - _t0).count()
                    << "ms write="
                    << std::chrono::duration_cast<std::chrono::milliseconds>(_t2 - _t1).count()
                    << "ms\n";
                return true;   // interactive step done -- fast path, no QC
            }

            // Final anchor. By now the finalize worker has been joined (the
            // controller joins it before the final step -- see main.cpp), so
            // job.tmpl carries the squared/absval scalars AND is safe to touch.
            // Fold the separately-accumulated anchors into it, write the FULL
            // file (all anchors + absval), promote, then run deferred QC.
            for (auto& kv : job.anchorAccum)
                job.tmpl.raw_anchors[kv.first] = std::move(kv.second);
            job.anchorAccum.clear();

            template_io::write_template_binfile(job.provisionalPath.string(), job.tmpl);
            job.viewerTemplatePath = job.provisionalPath;
            {
                std::error_code ec;
                std::filesystem::remove(job.templatePath, ec);
                std::filesystem::rename(job.provisionalPath, job.templatePath, ec);
                if (ec) {
                    std::cerr << "[templates] WARNING: could not promote "
                        << job.provisionalPath.string() << " -> "
                        << job.templatePath.string() << ": " << ec.message() << "\n";
                }
                else {
                    job.viewerTemplatePath = job.templatePath;
                }

                // ---- deferred QC over every anchor (R + the sequence) ----
                // For each anchor: take a fresh copy of the FINALIZED base
                // (job.tmpl -- carries absval in its scalars) and a fresh copy
                // of the R beats, then align in scoring mode so the raw scalars
                // + beats are put in the anchor frame. absval scalars ride
                // along unchanged. Score with writeEcgSQICsv.
                writeEcgSQICsv(job.cfg, job.stem + "_R_PEAK",
                    job.tmpl, job.beats, job.samplingRate);

                for (AnchorType a : seq) {
                    template_io::TemplateFile scoreT = job.tmpl;    // has absval in scalars
                    template_io::BeatsFile scoreB = job.beats;   // R beats, to be shifted
                    alignTemplatesFromCache(scoreT, scoreB, job.rates, a, /*forScoring=*/true);
                    writeEcgSQICsv(job.cfg, job.stem + "_" + anchorName(a),
                        scoreT, scoreB, job.samplingRate);
                }
            }
            return true;
        }
        catch (const std::exception& e) { job.error = e.what(); return false; }
    }
}  // namespace post_process_detail