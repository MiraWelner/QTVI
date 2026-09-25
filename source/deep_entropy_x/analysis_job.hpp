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

#include "template_generation/template_io.hpp"
#include "template_generation/build_bins.hpp"
#include "template_generation/premark_beats.hpp"
#include "template_morphology_grouping/bank_reload.hpp"
#include "fiducial_marker_finding/ppg_derivative.hpp"

#include "annealing/anneal_handler.hpp"
#include "config_file_handling/config.hpp"
#include "fiducial_marker_finding/alignment.hpp"
#include "logging/sqi_ecg.hpp"
#include "template_morphology_grouping/morphology_csv.hpp"
#include "template_morphology_grouping/envelope_report.hpp"
#include "template_morphology_grouping/beat_substitute.hpp"


namespace analysis_job {

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

    struct AnalysisJob {
        std::string stem;
        std::string fileID;
        double samplingRate = 0.0;
        SignalRates rates;           // full per-channel rate set for template pipeline
        std::filesystem::path rPeakPath, binsPath;
        std::filesystem::path annealedPath;
        std::vector<output_binfile_data> peakResults;
        template_io::TemplateFile tmpl;
        template_io::BeatsFile beats;
        // Per-bin TemplateInfo, carried from the fast build so finalize's
        // mergeTemplatesSlow can pack the squared/absval blocks into tmpl.
        // Non-const by reference there, so it has to live somewhere that
        // outlives prepare.
        std::vector<TemplateInfo> info;
        config_entry cfg{};
        bool use_consensus_peakfind_alg = true;    // = cfg.use_consensus_rpeak in prepare
        bool ecg1_inverted = false;
        bool ecg2_inverted = false;
        bool ecg3_inverted = false;
        std::string error;                          // set by finalize on failure
        template_io::TemplateFile r_aligned_template;      //  R-pass template to be reused by re-alignment
    };


    inline std::optional<AnalysisJob> prepare(const config_entry& cfg, const std::filesystem::path& binPath, bool ecg1_inverted, bool ecg2_inverted, bool ecg3_inverted)
    {
        const std::string stem = binPath.stem().string();
        const std::filesystem::path noise_bin_path = std::filesystem::path(cfg.noise_data_path) / (stem + "_noise_markings.bin");
        const std::filesystem::path annealedPath = std::filesystem::path(cfg.annealed_data_path) / (stem + "_annealed.bin");
        const std::filesystem::path rPeakPath = std::filesystem::path(cfg.r_peak_data_path) / (stem + "_peak_locations_all_beats.bin");
        const std::filesystem::path binsPath = std::filesystem::path(cfg.template_path) / (stem + "_bins.bin");

        anneal_one_file(binPath, noise_bin_path, annealedPath, cfg.bin_size_minutes, cfg.waveform_highpass_hz, ecg1_inverted, ecg2_inverted, ecg3_inverted);

        AnalysisJob job;
        job.stem = stem;
        job.fileID = stem;
        job.samplingRate = cfg.ecg_upsample_rate;
        //rate of 0 means absent channel
        job.rates = SignalRates{
            cfg.ecg_upsample_rate,
            cfg.ppg_upsample_rate,
            cfg.abp_upsample_rate,
            cfg.art_upsample_rate,
            cfg.art_pulm_upsample_rate,
            cfg.region_around_Rpeak_for_morphology_split,
            cfg.region_around_PPGPeak_for_morphology_split
        };

        job.cfg = cfg;
        job.use_consensus_peakfind_alg = cfg.use_consensus_rpeak;
        job.ecg1_inverted = ecg1_inverted;
        job.ecg2_inverted = ecg2_inverted;
        job.ecg3_inverted = ecg3_inverted;
        job.rPeakPath = rPeakPath;
        job.binsPath = binsPath;
        job.annealedPath = annealedPath;

        AnnealedData annealedData = read_input_binfile(annealedPath.string());


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
        job.peakResults = create_ecg_ppg_pairs_raw(std::move(annealedData.bins), true, stem, cfg, annealedData.ecg1_inverted, annealedData.ecg2_inverted, annealedData.ecg3_inverted);


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
        tbank::setMinBeats(cfg.min_beats_template_ecg, cfg.min_beats_template_ppg);//min beats for displayed templates in the viewer loaded from config

        // MORPHOLOGY SPLIT FLOORS. setMatchFloors takes both or neither: it
        // returns false and sets NOTHING if either value is outside (0, 1],
        // and a blank config cell arrives here as 0.0. So a refused pair
        // leaves both floors at their defaults, which changes how every bin in
        // the record is partitioned -- reported rather than silent.
        if (!tbank::setMatchFloors(cfg.ecg_match_floor, cfg.ppg_match_floor)) {
            std::cerr << "  [morphology] REFUSED match floors from config.csv (ecg="
                << cfg.ecg_match_floor << ", ppg=" << cfg.ppg_match_floor
                << ") -- both must be in (0, 1]. Keeping defaults ecg="
                << tbank::matchFloorEcg() << ", ppg="
                << tbank::matchFloorPpg() << "\n";
        }

        // Pulse QC threshold: unset keeps the default, unusable is refused
        // rather than clamped. THE ONLY CALLER of setFitErrorPct -- there used
        // to be an unconditional call above as well, which applied the value
        // before this block could refuse it, so "REFUSED ... Keeping X" could
        // print after X had already been replaced.
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

        // THE PRIOR SPLIT IS READ BEFORE THE BUILD. <stem>_templates.bin is
        // rewritten by the build below, inside morphology_csv::writeTemplatesBin,
        // so reading it afterwards would read this run's own output.
        const std::filesystem::path splitPath = std::filesystem::path(cfg.template_path) / (stem + "_templates.bin");
        bank_reload::SplitArchive priorSplit = bank_reload::readSplit(splitPath.string());

        FastTemplateBuild fast = buildTemplatesAndBeatsFast(job.peakResults, job.rates, noise_bin_path.string());
        if (fast.tmpl.bins.empty()) {
            std::cerr << "  no bins for " << stem << " (recording shorter than one bin?); skipping.\n";
            return std::nullopt;
        }
        job.tmpl = std::move(fast.tmpl);
        job.beats = std::move(fast.beats);
        job.info = std::move(fast.info);

        // AFTER the fresh build, because it overwrites what that build
        // partitioned; BEFORE the r_aligned_template snapshot, so the R frame
        // every anchor aligns from carries the reloaded banks.
        {
            const bank_reload::SplitReport rep =
                bank_reload::applySplit(priorSplit, job.tmpl);
            bank_reload::printReport(rep);
        }

        job.r_aligned_template = job.tmpl;      // snapshot R frame (one copy, at prep time)

        // Each anchor aligns FROM r_aligned_template, never from the previous
        // one, so the calls compose. R_PEAK is in the list: chFor short-circuits
        // it to the base, but its PER-SLOT averages are not a no-op and R is
        // what the grid draws on Automatic.
        for (AnchorType a : anchor_view::anchor_array) {
            template_io::TemplateFile atmpl = job.r_aligned_template;
            alignTemplatesFromCache(atmpl, job.beats, job.rates, a);

            const int tag = static_cast<int>(a);

            auto it = atmpl.raw_anchors.find(tag);
            if (it != atmpl.raw_anchors.end())
                job.tmpl.raw_anchors[tag] = std::move(it->second);

            // The per-slot averages for this anchor. `atmpl` dies at the end of
            // this iteration, so anything left in it is lost.
            auto bit = atmpl.bank_anchors.find(tag);
            if (bit != atmpl.bank_anchors.end())
                job.tmpl.bank_anchors[tag] = std::move(bit->second);
        }
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
            // POSITION IS LOAD-BEARING: these writes describe the RAW pass and
            // must stay above augment_ecg_ppg_pairs_sqabs, which overwrites the
            // beat lists in peakResults.
            write_output_binfile(job.rPeakPath.string(), job.peakResults);
            const std::filesystem::path csvDir = job.cfg.r_peak_data_path;
            const std::filesystem::path rPeakCsv = csvDir / (job.stem + "_peak_locations_all_beats.csv");
            write_output_csvfile(rPeakCsv.string(), job.peakResults, job.fileID, job.samplingRate);

            if (!job.cfg.template_path.empty()) {
                const std::string ftsPath = job.cfg.template_path + "/" + job.stem + "_pq_and_qrs_data.csv";
                normalize_features::writeFeatureTimeSeriesCsv(ftsPath, job.stem, job.peakResults, job.rates.ecg, pol);
                envelope_report::writeEnvelopeReport(job.cfg.template_path, job.stem, job.tmpl.bins, job.beats, job.rates.ecg, pol);
            }

            augment_ecg_ppg_pairs_sqabs(job.peakResults, job.use_consensus_peakfind_alg, job.fileID, job.samplingRate, job.cfg, job.ecg1_inverted, job.ecg2_inverted, job.ecg3_inverted);
            mergeTemplatesSlow(job.peakResults, job.tmpl, job.info, job.rates);
            premark::runAll(job.beats, job.tmpl, job.rates.ecg, pol, job.cfg.quality_metric, job.stem);
            writeEcgSQICsv(job.cfg, job.stem + "_R_PEAK", job.tmpl, job.beats, job.samplingRate, pol);
            std::cout << "Processing Squared and Absolute Value Templates (slow) for " << job.stem << "\n";
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
        try {
            template_io::write_template_binfile(job.binsPath.string(), job.tmpl);
            std::cerr << "  [templates] wrote " << job.binsPath.string()
                << " with the operator's confirmations\n";
            return true;
        }
        catch (const std::exception& e) {
            job.error = e.what();
            std::cerr << "  [templates] ERROR: could not write "
                << job.binsPath.string() << ": " << e.what()
                << " -- this record has no templates file\n";
            return false;
        }
    }

}  // namespace analysis_job