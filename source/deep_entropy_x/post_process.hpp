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
#include "template_morphology_grouping/envelope_report.hpp"

namespace post_process_detail {

    // The non-R alignments. R is the primary build (already in job.tmpl), so
    // it is NOT in this list.
    //
    // NO LONGER A CYCLE. These used to be applied one per "Finish and Next",
    // each step regenerating a provisional templates file and reloading the
    // viewer on a different alignment -- four windows to mark one subject.
    // buildAllAnchors() now aligns all of them up front and writes ONE file
    // carrying every block, and the viewer holds all four at once (see
    // anchor_view.hpp). The order here is only the order they are built in.
    inline const std::vector<AnchorType>& anchorSequence() {
        static const std::vector<AnchorType> seq = {
            AnchorType::P_ONSET, AnchorType::Q_ONSET, AnchorType::J_POINT,
        };
        return seq;
    }

    // Short label for the corner readout / CSV tags.
    inline const char* anchorName(AnchorType a) {
        switch (a) {
        case AnchorType::P_ONSET: return "P_ONSET";
        case AnchorType::Q_ONSET: return "Q_ONSET";
        case AnchorType::R_PEAK:  return "R_PEAK";
        case AnchorType::J_POINT: return "J_POINT";
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
    // Callers MUST join the finalize worker before invoking buildAllAnchors --
    // both mutate job.tmpl and would race if run concurrently. (Also documented
    // at buildAllAnchors' declaration below.)
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
        // (qAlignPath removed: the second-pass provisional it named was never
        //  assigned by anything, so main.cpp's two cleanup blocks for it could
        //  never fire. buildAllAnchors writes provisionalPath and renames it
        //  over templatePath.)
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
        // THE R-PASS BEATS, and they are not optional. augment_ecg_ppg_pairs_sqabs
        // (in finalizeViewerJob) OVERWRITES per_channel_beats with the
        // squared/absval detection -- the same fact that forces the bin_archive
        // "R" checkpoint to run above it. buildAllAnchors runs after that
        // worker is joined, so aligning from job.beats finds no CH1/CH2/CH3
        // beats at all: every bin is skipped, every anchor's store is written
        // all-empty, chFor falls back to the R base for all four alignments,
        // and the close-up shows one waveform under every bar. Cleared by
        // buildAllAnchors once the alignments exist.
        template_io::BeatsFile beatsR;        //  R-pass beats, pristine

        // (anchorAccum removed with the cycle. It staged each anchor's blocks
        //  away from job.tmpl so the interactive steps could run while the
        //  finalize worker was still mutating it. buildAllAnchors runs once,
        //  after that worker is joined, so it folds straight into job.tmpl.)
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

        // ---- SECTION 4.6 MORPHOLOGY THRESHOLDS, FROM config.csv ----------
        //
        // Applied here, once, before anything partitions. Every spawn and every
        // merge in the record turns on these two numbers, so they must be in
        // force before the first bin is built -- setting them later would leave
        // earlier bins partitioned against the defaults, with nothing on disk
        // saying which bins used which floor.
        //
        // EACH ONE FALLS BACK INDEPENDENTLY. A blank cell reaches here as 0.0,
        // and pairing it with a configured value would fail validation and
        // refuse BOTH -- so setting only the PPG floor, which is the likelier
        // thing to want, would silently do nothing. Each unset floor keeps its
        // own current value instead.
        //
        // AN UNUSABLE VALUE IS REFUSED, NOT CLAMPED. A floor of 0 accepts every
        // beat against every template: one morphology per bin, no ectopy ever
        // separated, and no error anywhere to explain it. Above 1 is the mirror
        // image -- correlation cannot exceed 1, so everything spawns. Either
        // way the defaults stand and the line below says so.
        {
            const bool have_ecg = (cfg.ecg_match_floor != 0.0);
            const bool have_ppg = (cfg.ppg_match_floor != 0.0);
            const double fe = have_ecg ? cfg.ecg_match_floor
                : tbank::matchFloorEcg();
            const double fp = have_ppg ? cfg.ppg_match_floor
                : tbank::matchFloorPpg();

            const char* src_ecg = have_ecg ? "config" : "default";
            const char* src_ppg = have_ppg ? "config" : "default";

            if (!tbank::setMatchFloors(fe, fp)) {
                std::cerr << "  [4.6] REFUSED match floors ECG " << fe
                    << " / PPG " << fp
                    << " -- both must be in (0, 1]. Keeping ECG "
                    << tbank::matchFloorEcg() << " / PPG "
                    << tbank::matchFloorPpg() << "\n";
            }
            else {
                std::cerr << "  [4.6] match floors ECG "
                    << tbank::matchFloorEcg() << " (" << src_ecg
                    << ")  PPG " << tbank::matchFloorPpg()
                    << " (" << src_ppg << ")\n";
            }

            // Minimum beats per template. 0 is a legal, meaningful value --
            // "no minimum" -- so unlike the floors there is nothing to refuse
            // except a negative, and the loader has already clamped that.
            if (tbank::setMinBeats(cfg.min_beats_template_ecg,
                cfg.min_beats_template_ppg)) {
                if (tbank::minBeatsEcg() > 0 || tbank::minBeatsPpg() > 0)
                    std::cerr << "  [4.6] min beats per template: ECG "
                    << tbank::minBeatsEcg() << "  PPG "
                    << tbank::minBeatsPpg()
                    << " -- templates below this are flagged "
                    "too_few_beats and not displayed\n";
                else
                    std::cerr << "  [4.6] no minimum beats per template "
                    "(min_beats_template_ecg/ppg unset)\n";
            }

            // Pulse QC threshold, same treatment: unset keeps the default,
            // unusable is refused rather than clamped.
            if (cfg.ppg_fit_error_pct == 0.0) {
                std::cerr << "  [pulseqc] ppg_fit_error_pct absent from "
                    "config.csv; using default "
                    << 100.0 * pulse_qc::fitErrorFraction() << "%\n";
            }
            else if (!pulse_qc::setFitErrorPct(cfg.ppg_fit_error_pct)) {
                std::cerr << "  [pulseqc] REFUSED ppg_fit_error_pct="
                    << cfg.ppg_fit_error_pct << " -- must be in (0, 100]. "
                    "Keeping " << 100.0 * pulse_qc::fitErrorFraction()
                    << "%\n";
            }
            else {
                std::cerr << "  [pulseqc] pulse fit error threshold "
                    << 100.0 * pulse_qc::fitErrorFraction()
                    << "% (config)\n";
            }
        }
        FastTemplateBuild fast = buildTemplatesAndBeatsFast(job.peakResults, job.rates);
        if (fast.tmpl.bins.empty()) {
            std::cerr << "  no bins for " << stem << " (recording shorter than one bin?); skipping.\n";
            return std::nullopt;   // main.cpp prints "prep failed or skipped"
        }
        job.tmpl = std::move(fast.tmpl);
        job.beats = std::move(fast.beats);
        job.info = std::move(fast.info);
        job.tmplR = job.tmpl;      // snapshot R frame (one copy, at prep time)
        job.beatsR = job.beats;    // and the beats it was built from -- see the field note

        // The R-pass checkpoints (bin archive, feature time series, envelope
        // report) used to run here. They are pure output and cost minutes, so
        // they moved to finalizeViewerJob, which runs concurrently with
        // marking -- see the note at their new home.
        template_io::write_template_binfile(provisionalPath.string(), job.tmpl);
        job.viewerTemplatePath = provisionalPath;
        job.needsFinalize = true;
        return job;
    }

    // Runs on a worker thread. Must not touch Qt. Stores any error in
    // job.error rather than throwing across the thread boundary.
    //
    // Callers MUST join this worker before invoking buildAllAnchors: both
    // operate on job.tmpl and would race if allowed to run concurrently. That
    // discipline is already documented at buildAllAnchors' declaration.
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
            // ---- R-PASS CHECKPOINTS, moved off the critical path ----------
            // These three were in prepareViewerJob, which BLOCKS the marking UI
            // from opening. They are pure output -- nothing downstream in
            // prepare reads them -- and poolBinQuality alone runs computeEcgSQI
            // per beat, per channel, per bin, which on a 40-bin record with
            // ~1000 beats a bin is ~120k evaluations. That was the wait.
            //
            // POSITION IS LOAD-BEARING: above augment_ecg_ppg_pairs_sqabs.
            // All three describe the R pass, and augment overwrites the beat
            // lists. Moving them below it would archive the squared/absval
            // detection under the label "R".
            if (!job.cfg.bin_archive_path.empty()) {
                const bool ok = bin_archive::writeBinFeatureArchive(
                    job.cfg.bin_archive_path, job.stem, job.tmpl.bins,
                    job.rates.ecg, "R", &job.beats);
                if (!ok)
                    std::cerr << "  [bin_archive] " << job.stem
                    << ": could not write checkpoint to "
                    << job.cfg.bin_archive_path << "\n";
                else
                    std::cerr << "  [bin_archive] " << job.stem
                    << ": wrote checkpoint (" << job.tmpl.bins.size() << " bins)\n";

                // Section 5.5 length/area/volume time series, from the SAME
                // pre-deformation R-pass data. job.peakResults still holds the
                // raw per-channel ECG + detected R-peaks the proportional
                // segmenter needs (job.tmpl has only averaged templates, which
                // cannot be re-segmented).
                const std::string ftsPath =
                    job.cfg.bin_archive_path + "/" + job.stem + "_feature_timeseries.csv";
                const bool okf = normalize_features::writeFeatureTimeSeriesCsv(
                    ftsPath, job.stem, job.peakResults, job.rates.ecg);
                if (!okf)
                    std::cerr << "  [feature_ts] " << job.stem
                    << ": could not write " << ftsPath << "\n";
                else
                    std::cerr << "  [feature_ts] " << job.stem
                    << ": wrote length/area/volume series\n";

                // Section 4.7 dynamic envelopes, per beat per segment per
                // channel. Serial by construction (see envelope_report.hpp):
                // the rolling windows are sequential per channel, so this is
                // the one checkpoint here that must not be parallelised over
                // bins. Always written -- a report that only appears when
                // someone remembers a flag is missing from the runs that
                // matter.
                const bool oke = envelope_report::writeEnvelopeReport(
                    job.cfg.bin_archive_path, job.stem, job.tmpl.bins, job.beats,
                    job.rates.ecg);
                if (!oke)
                    std::cerr << "  [envelopes] " << job.stem
                    << ": could not write envelope report to "
                    << job.cfg.bin_archive_path << "\n";
            }


            // Squared/absval R-peak detection on ECG channels, then slow
            // templating to pack the two extra per-bin blocks (squared,
            // absval) plus their SAECG entries into job.tmpl. Runs entirely
            // on this worker; no Qt access. The worker is expected to have
            // been joined before buildAllAnchors is called (see comment
            // above), so no race with the anchor path here.
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
            // Safe here: premark mutates neither job.beats nor job.tmpl, so
            // it cannot race the mergeTemplatesSlow above (already complete).
            // It reads the RAW per-bin blocks, which augment/merge leave
            // intact. buildAllAnchors cannot overlap it either -- it runs only
            // after this whole worker is joined.
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

    // ONE CALL, EVERY ALIGNMENT. Replaces the per-step regenerateWithAnchor
    // cycle: aligns each anchor from the pristine R snapshot, folds every block
    // into job.tmpl, writes the full file once, promotes it, and runs the
    // deferred per-anchor QC.
    //
    // MUST be called AFTER the finalize worker is joined -- it touches
    // job.tmpl, which that worker mutates to pack the squared/absval scalars.
    // The old anchorAccum staging existed purely so the interactive steps could
    // run concurrently with finalize; with one call at one point there was
    // nothing to stage, so anchorAccum and anchorStep are gone from ViewerJob.
    //
    // Each anchor aligns FROM job.tmplR, never from the previous anchor:
    // alignTemplatesFromCache leaves the R base untouched between calls, so the
    // four compose. Returns false if nothing could be built.
    inline bool buildAllAnchors(ViewerJob& job)
    {
        try {
            const auto& seq = anchorSequence();
            auto _tAll = std::chrono::steady_clock::now();

            for (AnchorType anchor : seq) {
                auto _t0 = std::chrono::steady_clock::now();
                // job.beatsR, NOT job.beats: see the field note on beatsR. The
                // template copy is cheap next to the beats, which are passed by
                // reference and not mutated in the default non-scoring mode.
                template_io::TemplateFile atmpl = job.tmplR;
                alignTemplatesFromCache(atmpl, job.beatsR, job.rates, anchor);

                const int anchorTag = static_cast<int>(anchor);
                size_t filled = 0;
                auto it = atmpl.raw_anchors.find(anchorTag);
                if (it != atmpl.raw_anchors.end()) {
                    // COUNT WHAT ACTUALLY ALIGNED. alignTemplatesFromCache
                    // sizes the store to bins up front and leaves skipped
                    // slots default-empty, so the map entry existing says
                    // nothing about whether any bin was aligned -- an all-empty
                    // store round-trips through the file perfectly and then
                    // silently reads back as "use the R base". That is the one
                    // failure mode of this whole feature that looks like
                    // success, so it gets counted and reported rather than
                    // inferred from a timing line.
                    for (const auto& trip : it->second)
                        if (!trip[0].ecgTemplate.empty()
                            || !trip[1].ecgTemplate.empty()
                            || !trip[2].ecgTemplate.empty()) ++filled;
                    job.tmpl.raw_anchors[anchorTag] = std::move(it->second);
                }

                auto _t1 = std::chrono::steady_clock::now();
                std::cerr << "  [timing] anchor " << anchorName(anchor) << ": align="
                    << std::chrono::duration_cast<std::chrono::milliseconds>(_t1 - _t0).count()
                    << "ms, " << filled << "/" << job.tmpl.bins.size() << " bins\n";
                if (filled == 0) {
                    std::cerr << "  [anchors] WARNING: " << anchorName(anchor)
                        << " aligned 0 bins -- the marking window will show the "
                        << "R-aligned template under every landmark for this "
                        << "alignment.\n";
                }
            }

            // ONE WRITE. The old per-step write was trimmed to the R base plus
            // the single anchor the viewer was about to open, because the viewer
            // read exactly one; it now reads all four, so the file has to carry
            // all four.
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

                writeEcgSQICsv(job.cfg, job.stem + "_R_PEAK",
                    job.tmpl, job.beats, job.samplingRate);

                for (AnchorType a : seq) {
                    template_io::TemplateFile scoreT = job.tmpl;
                    template_io::BeatsFile scoreB = job.beats;
                    alignTemplatesFromCache(scoreT, scoreB, job.rates, a, /*forScoring=*/true);
                    writeEcgSQICsv(job.cfg, job.stem + "_" + anchorName(a),
                        scoreT, scoreB, job.samplingRate);
                }
            }
            // The beats snapshot exists only to feed the alignment above; it is
            // a full copy of the R-pass beats and there is no reason to carry it
            // for the length of the marking session.
            job.beatsR = template_io::BeatsFile{};

            std::cerr << "  [timing] all anchors + QC: "
                << std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::steady_clock::now() - _tAll).count() << "ms\n";
            return true;
        }
        catch (const std::exception& e) { job.error = e.what(); return false; }
    }

}  // namespace post_process_detail