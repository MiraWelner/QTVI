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
#include "peak_finding/peakfinding_io.hpp"
#include "peak_finding/make_beats.hpp" 

#include "template_generation/template_structs.hpp"
#include "template_generation/build_bins.hpp"
#include "template_generation/premark_beats.hpp"
#include "template_generation/bank_reload.hpp"
#include "template_generation/template_io.hpp"
#include "template_generation/envelope_report.hpp"
#include "template_generation/beat_substitute.hpp"
#include "annealing/beat_times.hpp"
#include "annealing/record_sleep.hpp"

#include "annealing/anneal_handler.hpp"
#include "config_file_handling/config.hpp"
#include "fiducial_marker_finding/alignment.hpp"
#include "fiducial_marker_finding/ppg_derivative.hpp"

#include "logging/sqi_ecg.hpp"


namespace analysis_job {

    // Short label for the corner readout / CSV tags.
    inline const char* anchorName(AnchorType a) {
        switch (a) {
        case AnchorType::P_ONSET: return "P_ONSET";
        case AnchorType::Q_ONSET: return "Q_ONSET";
        case AnchorType::R_PEAK:  return "R_PEAK";
        case AnchorType::J_POINT: return "J_POINT";
        case AnchorType::T_END:   return "T_END";
        }
        return "?";
    }

    struct AnalysisJob {
        std::string stem;
        std::string fileID;
        double samplingRate = 0.0;
        SignalRates rates;           // full per-channel rate set for template pipeline
        std::filesystem::path rPeakPath;
        std::filesystem::path annealedPath;
        std::vector<output_binfile_data> peakResults;
        template_structs::TemplateFile tmpl;
        template_structs::BeatsFile beats;
        // Per-bin TemplateInfo, carried from the fast build so finalize's
        // mergeTemplatesSlow can pack the squared/absval blocks into tmpl.
        // Non-const by reference there, so it has to live somewhere that
        // outlives prepare.
        std::vector<TemplateInfo> info;
        config_entry cfg{};
        bool ecg1_inverted = false;
        bool ecg2_inverted = false;
        bool ecg3_inverted = false;
        std::string error;                          // set by finalize on failure
        template_structs::TemplateFile r_aligned_template;      //  R-pass template to be reused by re-alignment
        // Two separate things, both built once in prepare and read by the viewer
        // for template_markings.csv. Neither is written to any .bin.
        //   beatTimes: recording seconds of every kept ECG beat, [lead][bin][row].
        //   sleep:     the record's staging; present() is false when it has none.
        beat_times::BeatTimes beatTimes;
        record_sleep::RecordSleep sleep;
        // Every alignment's per-row horizontal shift, in samples, exactly as
        // align_beat_matrix applied it to build that alignment's per-slot
        // averages: [anchor tag][bin][lead][kept row], NaN = row not moved.
        // Copied out of ecg_move_log, which clears its store when it writes the
        // move log -- before the viewer opens. Read by the viewer to draw an
        // individual beat on an anchored average's own frame.
        std::map<int, std::vector<std::array<std::vector<double>, 3>>> ecgRowShifts;
    };


    // highpass_enabled is the operator's checkbox; cfg supplies the cutoff.
    // Applied ahead of anneal_one_file, so ahead of the peak finding, the bins
    // and the templates. Nothing downstream re-filters.
    inline std::optional<AnalysisJob> prepare(const config_entry& cfg, const std::filesystem::path& binPath, bool ecg1_inverted, bool ecg2_inverted, bool ecg3_inverted, bool highpass_enabled)
    {
        const std::string stem = binPath.stem().string();
        const std::filesystem::path noise_bin_path = std::filesystem::path(cfg.noise_data_path) / (stem + "_noise_markings.bin");
        const std::filesystem::path annealedPath = std::filesystem::path(cfg.annealed_data_path) / (stem + "_annealed.bin");
        const std::filesystem::path rPeakPath = std::filesystem::path(cfg.r_peak_data_path) / (stem + "_peak_locations_all_beats.bin");

        const double highpassHz = highpass_enabled ? cfg.waveform_highpass_hz : 0.0;
        std::cerr << "  [highpass] cutoff " << highpassHz
            << " Hz, ECG1/2/3 and PPG, before annealing\n";
        anneal_one_file(binPath, noise_bin_path, annealedPath, cfg.bin_size_minutes, highpassHz, ecg1_inverted, ecg2_inverted, ecg3_inverted);

        AnalysisJob job;
        job.stem = stem;
        job.fileID = stem;
        job.samplingRate = cfg.ecg_upsample_rate;
        // Rate of 0 means absent channel. Built by name in
        // config_entry::signalRates(), not by a positional brace here.
        job.rates = cfg.signalRates();

        job.cfg = cfg;
        job.ecg1_inverted = ecg1_inverted;
        job.ecg2_inverted = ecg2_inverted;
        job.ecg3_inverted = ecg3_inverted;
        job.rPeakPath = rPeakPath;
        job.annealedPath = annealedPath;

        AnnealedData annealedData = read_input_binfile(annealedPath.string());


        channel_offset::set(cfg.training_log, stem);
        channel_offset::Result chOffPpg, chOffArt;
        const bool wantChannelOffset = (cfg.dataset_type == "CHAOS");
        if (wantChannelOffset) {
            {
                auto _cot0 = std::chrono::steady_clock::now();
                auto probeSegs = channel_offset::make_probe(
                    annealedData.bins, cfg.ecg_upsample_rate, cfg.ppg_upsample_rate);
                if (!probeSegs.empty()) {
                    auto probeResults = create_ecg_ppg_pairs_raw(
                        std::move(probeSegs), stem, cfg,
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
                    // pulse locations detected fresh from the channel's own
                    // derivative (ppg_deriv's E-0 census) -- purely as a
                    // repeatable per-pulse fiducial for
                    std::vector<std::vector<std::size_t>> abpFeet;
                    abpFeet.reserve(probeResults.size());
                    if (cfg.abp_upsample_rate > 0.0) {
                        const int minSep = std::max(1,
                            static_cast<int>(std::llround(0.25 * cfg.abp_upsample_rate)));
                        for (const auto& b : probeResults) {
                            std::vector<std::size_t> feet;
                            if (b.all_upsampled.size() > 33 && !b.all_upsampled[33].empty()) {
                                const std::vector<int> locs = ppg_deriv::derivativePulseLocations(b.all_upsampled[33], minSep);
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
            // MEASURED, NOT APPLIED. The signals stay as recorded; see
            // channel_offset.hpp. A confident lag moves the pulse window (via
            // job.rates, below) and the viewer's display, nothing else.
            channel_offset::write_log(chOffPpg, "PPG", /*append=*/false);
            channel_offset::write_log(chOffArt, "ARTERIAL", /*append=*/true);

            auto say = [&](const char* group, const channel_offset::Result& r) {
                std::cerr << "  [channel_offset] " << stem << ": " << group
                    << " lag " << r.lag_ms << " ms (ratio=" << r.ratio << ") -- "
                    << (r.ambiguous ? "ambiguous, not used"
                        : "used for the pulse window and display; signals not shifted")
                    << "\n";
                };
            say("PPG", chOffPpg);
            say("ABP/ART/ART_PULM", chOffArt);
        }

        // ---- THE RECORD'S LAGS ------------------------------------------
        // Measured for CHAOS only; a group with no data is "not measured",
        // which is method NONE, not a lag of 0 ms.
        channel_offset::RecordLag recordLag;
        if (wantChannelOffset) {
            recordLag.ppg = channel_offset::fromResult(chOffPpg,
                chOffPpg.n_r_peaks >= 2 && chOffPpg.n_feet >= 2);
            recordLag.arterial = channel_offset::fromResult(chOffArt,
                chOffArt.n_r_peaks >= 2 && chOffArt.n_feet >= 2);
        }
        job.rates.ppg_lag_ms = recordLag.ppg.usableLagMs();
        job.rates.arterial_lag_ms = recordLag.arterial.usableLagMs();

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
        job.peakResults = create_ecg_ppg_pairs_raw(std::move(annealedData.bins), stem, cfg, annealedData.ecg1_inverted, annealedData.ecg2_inverted, annealedData.ecg3_inverted);


        // create_ecg_ppg_pairs_raw doesn't carry the arterial pass-through
        // channels, so attach them here (parallel by bin index).
        for (size_t i = 0; i < job.peakResults.size() && i < nAnnealed; ++i) {
            job.peakResults[i].abpSignal = std::move(abpSlots[i]);
            job.peakResults[i].artSignal = std::move(artSlots[i]);
            job.peakResults[i].artPulmSignal = std::move(artpSlots[i]);
        }

        std::cerr << "  Processing Raw Templates (fast stage): " << stem << "\n";
        ecg_move_log::set(cfg.training_log, stem);   // per-beat vertical + horizontal move log
        ptt_log::set(cfg.training_log, stem);        // per-beat R->PPG transit times
        templates_io::set(cfg.template_path, stem);
        // AFTER set(), which clears it: the archive this build writes carries
        // this record's lags as a trailer.
        templates_io::setRecordLag(recordLag);
        tbank::setMinBeats(cfg.min_beats_template_ecg, cfg.min_beats_template_ppg);//min beats for displayed templates in the viewer loaded from config

        // MORPHOLOGY SPLIT THRESHOLDS, ONE SIDE AT A TIME. This was a single
        // setMatchFloors(ecg, ppg) that refused BOTH values if either was out
        // of range, so a config with a good ECG threshold and a blank PPG cell
        // silently reverted both. Reported per side now, because a threshold
        // that changes how every bin in the record is partitioned should not
        // be able to revert quietly.
        //
        // A BLANK CELL IS 0.0 AND IS REFUSED, which is the intended outcome:
        // 0.0 is not "no threshold", it is a threshold every beat clears, and
        // accepting it would disable the split while looking configured.
        tbank::setExcludePostEctopic(cfg.exclude_beat_after_ectopic);
        std::cerr << "  [post-ectopic] beat after a PVC/PAC/VT or premature beat: "
            << (cfg.exclude_beat_after_ectopic
                ? "kept in the sinus partition and RR, left out of the average and QT"
                : "treated as an ordinary beat")
            << " (exclude_beat_after_ectopic)\n";
        if (!tbank::setMorphThresholdEcg(cfg.morph_threshold_ecg)) {
            std::cerr << "  [morphology] morph_threshold_ecg="
                << cfg.morph_threshold_ecg << " not usable (need (0, 1]); "
                "keeping " << tbank::morphThresholdEcg() << "\n";
        }
        if (!tbank::setMorphThresholdPpg(cfg.morph_threshold_ppg)) {
            std::cerr << "  [morphology] morph_threshold_ppg="
                << cfg.morph_threshold_ppg << " not usable (need (0, 1]); "
                "keeping " << tbank::morphThresholdPpg() << "\n";
        }
        std::cerr << "  [morphology] split thresholds ecg="
            << tbank::morphThresholdEcg() << ", ppg="
            << tbank::morphThresholdPpg() << "\n";

        // Pulse QC threshold: unset keeps the default, unusable is refused
        // rather than clamped. THE ONLY CALLER of setCorrFloor -- there used
        // to be an unconditional call above as well, which applied the value
        // before this block could refuse it, so "REFUSED ... Keeping X" could
        // print after X had already been replaced.
        // A CORRELATION, NOT A PERCENT. Was ppg_fit_error_pct in (0, 100],
        // printed as 100.0 * the stored fraction; now pulse_qc_corr_floor in
        // (0, 1], printed as itself.
        if (cfg.pulse_qc_corr_floor == 0.0) {
            std::cerr << "  [pulseqc] pulse_qc_corr_floor absent from "
                "config.csv; using default r >= "
                << pulse_qc::corrFloor() << "\n";
        }
        else if (!pulse_qc::setCorrFloor(cfg.pulse_qc_corr_floor)) {
            std::cerr << "  [pulseqc] REFUSED pulse_qc_corr_floor="
                << cfg.pulse_qc_corr_floor << " -- must be in (0, 1]. "
                "Keeping r >= " << pulse_qc::corrFloor() << "\n";
        }
        else {
            std::cerr << "  [pulseqc] pulse correlation floor r >= "
                << pulse_qc::corrFloor() << " (config)\n";
        }

        // THE PRIOR SPLIT IS READ BEFORE THE BUILD. <stem>_templates.bin is
        // rewritten by the build below, inside morphology_csv::writeTemplatesBin,
        // so reading it afterwards would read this run's own output.
        const std::filesystem::path splitPath = std::filesystem::path(cfg.template_path) / (stem + "_templates.bin");
        bank_reload::SplitArchive priorSplit = bank_reload::readSplit(splitPath.string());

        // The prior split is applied INSIDE the build, before the archive is
        // rewritten -- see bank_reload.hpp for the two gates and why.
        FastTemplateBuild fast = buildTemplatesAndBeatsFast(job.peakResults, job.rates, noise_bin_path.string(), &priorSplit);
        if (fast.tmpl.bins.empty()) {
            std::cerr << "  no bins for " << stem << " (recording shorter than one bin?); skipping.\n";
            return std::nullopt;
        }
        job.tmpl = std::move(fast.tmpl);
        job.tmpl.lag = recordLag;   // for the viewer: display and the markings file
        job.beats = std::move(fast.beats);
        job.info = std::move(fast.info);

        // BEAT TIMES, HERE: job.peakResults' R peaks and splices are exactly the
        // ones this build sliced on, and job.info still holds each bin's row ->
        // slice map. The original data .bin's header gives the ECG1 rate and
        // length the annealer cut the bins from; see beat_times.hpp.
        //
        // SLEEP, SEPARATELY, from the same file's stage block. Independent of the
        // beat times: a failure in one leaves the other intact.
        {
            const data_bin_header::Header hdr = data_bin_header::read(binPath.string());
            if (!hdr.ok)
                std::cerr << "  [beat times] could not read the header of " << binPath.string()
                << "; times use ecg_upsample_rate, no sleep staging\n";
            try {
                std::vector<std::array<std::vector<size_t>, 3>> sliceOfRow(job.info.size());
                for (size_t i = 0; i < job.info.size(); ++i)
                    sliceOfRow[i] = job.info[i].ecg_slice_of_row;
                const bool useHdr = hdr.ok && hdr.ecgRateHz() > 0.0;
                job.beatTimes = useHdr
                    ? beat_times::build(job.peakResults, sliceOfRow, hdr.ecgRateHz(), hdr.ecgLength())
                    : beat_times::build(job.peakResults, sliceOfRow, job.rates.ecg);
            }
            catch (const std::exception& e) {
                std::cerr << "  [beat times] failed (" << e.what()
                    << "); template_start_s / template_end_s will be blank\n";
                job.beatTimes = beat_times::BeatTimes{};
            }
            try {
                job.sleep = record_sleep::read(binPath.string(), hdr);
            }
            catch (const std::exception& e) {
                std::cerr << "  [sleep] failed (" << e.what() << ")\n";
                job.sleep = record_sleep::RecordSleep{};
            }
            std::cerr << "  [sleep] " << (job.sleep.present()
                ? "staging found, epoch " + std::to_string(job.sleep.epoch_sec) + " s"
                : std::string("no sleep staging in this record; pct_* columns blank")) << "\n";
        }

        // Applied inside the build above; reported here. The reloaded banks
        // are already in job.tmpl, so the r_aligned_template snapshot below,
        // which every anchor aligns from, carries them.
        bank_reload::printReport(priorSplit.rep);

        job.r_aligned_template = job.tmpl;      // snapshot R frame (one copy, at prep time)

        // Each anchor aligns FROM r_aligned_template, never from the previous
        // one, so the calls compose. R_PEAK is in the list: chFor short-circuits
        // it to the base, but its PER-SLOT averages are not a no-op and R is
        // what the grid draws on Automatic.
        for (AnchorType a : anchor_view::anchor_array) {
            template_structs::TemplateFile atmpl = job.r_aligned_template;
            alignTemplatesFromCache(atmpl, job.beats, job.rates, a);

            const int tag = static_cast<int>(a);

            // This anchor's per-row shifts, before ecg_move_log::write clears them.
            for (const ecg_move_log::HorizontalAnchor& h : ecg_move_log::g_horizontal)
                if (h.label == anchor_view::label(a)) { job.ecgRowShifts[tag] = h.bins; break; }

            auto it = atmpl.raw_anchors.find(tag);
            if (it != atmpl.raw_anchors.end())
                job.tmpl.raw_anchors[tag] = std::move(it->second);

            // The per-slot averages for this anchor. `atmpl` dies at the end of
            // this iteration, so anything left in it is lost.
            auto bit = atmpl.bank_anchors.find(tag);
            if (bit != atmpl.bank_anchors.end())
                job.tmpl.bank_anchors[tag] = std::move(bit->second);
        }

        // <stem>_beat_moves.csv, now that both halves exist: the vertical
        // shifts the fast build stashed, and every anchor's horizontal ones.
        ecg_move_log::write(job.rates.ecg);
        // Every pulse channel's transit times, stashed by the build.
        ptt_log::writeTransit();

        std::cerr.flush();
        return job;
    }

    inline void finalize(AnalysisJob& job)
    {
        // Leave one core for the Qt UI thread: without this, the OpenMP
        // regions below each take omp_get_max_threads() and the marking GUI
        // gets starved even though it is on another thread. omp_set_nested(0)
        // stops an inner parallel-for spawning a nested team inside the
        // outer one.
        const int hw = std::max(1, static_cast<int>(std::thread::hardware_concurrency()));
        const int workerThreads = std::max(1, hw - 1);
        omp_set_num_threads(workerThreads);
        omp_set_nested(0);
        const LeadPolarity pol{ { job.ecg1_inverted, job.ecg2_inverted, job.ecg3_inverted } };
        try {
            if (!job.cfg.template_path.empty()) {
                const std::string ftsPath = job.cfg.template_path + "/" + job.stem + "_pq_and_qrs_data.csv";
                normalize_features::writeFeatureTimeSeriesCsv(ftsPath, job.stem, job.peakResults, job.rates.ecg, pol);
                envelope_report::writeEnvelopeReport(job.cfg.template_path, job.stem, job.tmpl.bins, job.beats, job.rates.ecg, pol);
            }

            std::cout << "Processing Squared and Absolute Value Templates (slow) for " << job.stem << "\n";
            augment_ecg_ppg_pairs_sqabs(job.peakResults, job.fileID, job.samplingRate, job.cfg);

            // Below augment, so all 9 channel x method blocks are populated in
            // both files rather than just the 3 raw ones.
            write_output_binfile(job.rPeakPath.string(), job.peakResults);
            const std::filesystem::path rPeakCsv =
                std::filesystem::path(job.cfg.r_peak_data_path) / (job.stem + "_peak_locations_all_beats.csv");
            write_output_csvfile(rPeakCsv.string(), job.peakResults, job.fileID, job.samplingRate);

            mergeTemplatesSlow(job.peakResults, job.tmpl, job.info, job.rates);
            premark::runAll(job.beats, job.tmpl, job.rates.ecg, pol, job.cfg.training_log, job.stem);
            writeEcgSQICsv(job.cfg, job.stem + "_R_PEAK", job.tmpl, job.beats, job.samplingRate, pol);
        }
        catch (const std::exception& e) {
            job.error = e.what();
        }
        catch (...) {
            job.error = "unknown exception in finalize";
        }
    }


    struct BankSnapshot {
        std::array<tbank::TemplateBank, 3> ecg_bank;
        tbank::TemplateBank                ppg_bank;
    };

    inline bool commit(AnalysisJob& job, const std::vector<BankSnapshot>& banks)
    {
        if (banks.size() < job.tmpl.bins.size()) {
            // Not fatal -- the operator may have closed the viewer before every
            // bin was paged in -- but silence here would look like a record
            // whose later bins were all reviewed and found unconfirmed.
            std::cerr << "  [templates] note: " << banks.size() << " of "
                << job.tmpl.bins.size() << " bins came back from the viewer; "
                "the rest keep their generated banks\n";
        }
        for (size_t i = 0; i < job.tmpl.bins.size() && i < banks.size(); ++i) {
            job.tmpl.bins[i].ecg_bank = banks[i].ecg_bank;
            job.tmpl.bins[i].ppg_bank = banks[i].ppg_bank;
        }
        // ---- WHERE THE VERDICT GOES NOW ------------------------------
        //
        // <stem>_templates.bin, rewritten in place. It was <stem>_bins.bin,
        // whose Sections 3-5 carried the banks via tbank_ser -- the only place
        // confirmed_by_operator ever reached disk. That file also carried the
        // bin-wide averages, a cache of what _beats.bin can rebuild, so it is
        // retired and the sub-templates keep the archive named for them.
        //
        // POINTERS INTO job.tmpl, TAKEN AFTER THE FOLD ABOVE. The loop just
        // above is what puts the operator's banks into job.tmpl.bins, so a view
        // built before it would rewrite the archive with the generated
        // partition and report success.
        std::vector<templates_io::BinBanks> view(job.tmpl.bins.size());
        for (size_t i = 0; i < job.tmpl.bins.size(); ++i) {
            for (int c = 0; c < 3; ++c)
                view[i].chan[c] = &job.tmpl.bins[i].ecg_bank[c];
            view[i].chan[3] = &job.tmpl.bins[i].ppg_bank;
        }

        const std::filesystem::path splitPath =
            std::filesystem::path(job.cfg.template_path) / (job.stem + "_templates.bin");
        const templates_io::ConfirmPatchReport rep =
            templates_io::applyOperatorConfirmations(splitPath.string(), view);
        templates_io::printConfirmPatchReport(rep);

        // A MISSING ARCHIVE IS NOT A FAILED COMMIT. A record whose build
        // produced no templates has nothing to patch; only a file that exists
        // and could not be rewritten is worth failing on, because that is the
        // case where the verdict was lost.
        if (!rep.error.empty()) job.error = rep.error;
        return rep.written || !rep.present;
    }

}  // namespace analysis_job