/**
 * @file   make_beats.hpp
 * @brief  Find the fiduical markers that mark the start and the end of each beat - ie the r peaks and ppg mins
 *
 *
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-26
 */

#pragma once

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include "peakfinding_structs.hpp"
#include "SegmentPPG.hpp"
#include "JoinedRR.hpp"
#include "config_file_handling/config.hpp"

static inline void run_rpeak_detection(const std::vector<double>& sig, double ecgRate,
    const std::string& fileID, std::vector<std::size_t>& rIndex, config_entry cfg, bool inverted) {

    rIndex.clear();

    // A zero-variance channel is file_to_bin's "channel absent" placeholder.
    // Leaving rIndex empty here is what CreateEcgTemplatesFast's
    // (signal non-empty AND chN.raw non-empty) presence test relies on.
    if (std_dev(sig) == 0.0) return;

    try {
        if (cfg.use_consensus_rpeak) {
            JoinedRRResult jrr = JoinedRR_full(sig, ecgRate, fileID, inverted);
            rIndex = std::move(jrr.peaks);
        }
        else {
            // single-detector fallback: custom rpeakdetect at its default threshold
            RPeakDetectResult r = rpeakdetect(sig, ecgRate, 0.2, inverted);
            rIndex = std::move(r.r_peak_index);
        }
    }
    catch (...) {
        rIndex.clear();
        std::fprintf(stderr, "  [rpeak] detector threw for %s -- no peaks for this "
            "channel/bin\n", fileID.c_str());
        std::fflush(stderr);
    }
}


static inline void detect_channel_raw(ChannelRPeaks& result, const std::vector<double>& signal, double ecgRate,
    const std::string& fileID, config_entry cfg, bool inverted)
{
    if (signal.empty()) return;
    run_rpeak_detection(signal, ecgRate, fileID, result.raw, cfg, inverted);
}


/* Squaring and rectification both map an R-peak to a positive excursion
   whatever the lead polarity, so the peak is a local maximum on both derived
   signals and no inverted-lead flag applies here. */
static inline void detect_channel_sqabs(ChannelRPeaks& result, const std::vector<double>& signal, double ecgRate,
    const std::string& fileID, config_entry cfg)
{
    if (signal.empty()) return;

    /* Method 2: squared signal */
    result.squared_signal.resize(signal.size());
    for (std::size_t i = 0; i < signal.size(); ++i) {
        result.squared_signal[i] = signal[i] * signal[i];
    }
    run_rpeak_detection(result.squared_signal, ecgRate, fileID, result.squared,
        cfg, /*inverted=*/false);

    /* Method 3: absolute value signal */
    result.absval_signal.resize(signal.size());
    for (std::size_t i = 0; i < signal.size(); ++i) {
        result.absval_signal[i] = std::fabs(signal[i]);
    }
    run_rpeak_detection(result.absval_signal, ecgRate, fileID, result.absval,
        cfg, /*inverted=*/false);
}

inline std::vector<output_binfile_data> create_ecg_ppg_pairs_raw(std::vector<AnnealedSegment> annealedSegments,
    std::string fileID, config_entry cfg,
    bool ecg1_inverted, bool ecg2_inverted, bool ecg3_inverted) {

    std::vector<output_binfile_data> data(annealedSegments.size());
    if (cfg.use_consensus_rpeak) {
        std::cout << "Using consensus peak finding method\n";
    }
    else {
        std::cout << "Using only one r peak detection method\n";
    }

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < static_cast<int>(annealedSegments.size()); ++i) {

        AnnealedSegment seg = std::move(annealedSegments[i]);
        auto& d = data[i];

        d.index = static_cast<int>(i);

        d.ppgSignal = seg.ppg_signal;
        d.ecgSignal = seg.ecg_signal_1;
        d.ecgSignal2 = seg.ecg_signal_2;
        d.ecgSignal3 = seg.ecg_signal_3;

        // CHAOS: subtract the whole-signal median ONCE, here, before any
        // detection or alignment runs -- from all three ECG channels AND
        // from PPG. Removes the per-segment DC pedestal globally (cheaper
        // and simpler than doing it per beat, per method, per bin
        // downstream).
        if (cfg.dataset_type == "CHAOS") {
            auto subtract_median = [](std::vector<double>& sig) {
                if (sig.empty()) return;
                std::vector<double> vals;
                vals.reserve(sig.size());
                for (double v : sig) if (!std::isnan(v)) vals.push_back(v);
                if (vals.empty()) return;
                const size_t m = vals.size() / 2;
                std::nth_element(vals.begin(), vals.begin() + m, vals.end());
                double med = vals[m];
                if (vals.size() % 2 == 0) {
                    const double hi = med;
                    med = 0.5 * (*std::max_element(vals.begin(), vals.begin() + m) + hi);
                }
                for (double& v : sig) if (!std::isnan(v)) v -= med;
                };
            subtract_median(d.ecgSignal);
            subtract_median(d.ecgSignal2);
            subtract_median(d.ecgSignal3);
            subtract_median(d.ppgSignal);

        }

        d.ppg_bin_indexs = std::move(seg.ppg_bin_indexs);
        d.ecg_bin_indexs = std::move(seg.ecg_bin_indexs);

        d.all_upsampled = std::move(seg.all_upsampled);
        d.all_raw_pairs_flat = std::move(seg.all_raw_pairs_flat);

        const bool hasPPG = !d.ppgSignal.empty();

        /* Step 1 - PPG pulse segmentation. Supplies the valley list that
           pairRtoPPGBeat reads below. */
        if (hasPPG) {
            try {
                SegmentPPGResult ppgResult = SegmentPPG(d.ppgSignal, cfg.ppg_upsample_rate);
                d.ppgMinAmps = ppgResult.minAmps;
                d.ppgMaxAmps = ppgResult.maxAmps;
            }
            catch (const std::exception& e) {
                // REPORTED, not just absorbed. The valleys are no longer a
                // precondition for pulse TEMPLATES (see the has_ppg note in
                // make_averaged_templates.hpp), but they are still what
                // pairRtoPPGBeat needs below, so losing them costs the R-to-pulse
                // pairing and has to be visible rather than inferred from a
                // missing column three stages later.
                std::fprintf(stderr,
                    "  [ppg-seg] SegmentPPG failed (%s) -- no PPG valleys for "
                    "this bin; R-to-pulse pairing skipped, pulse templates "
                    "unaffected\n", e.what());
                std::fflush(stderr);
                d.ppgMinAmps.clear();
                d.ppgMaxAmps.clear();
            }
            catch (...) {
                std::fprintf(stderr,
                    "  [ppg-seg] SegmentPPG failed (unknown exception) -- no PPG"
                    " valleys for this bin\n");
                std::fflush(stderr);
                d.ppgMinAmps.clear();
                d.ppgMaxAmps.clear();
            }
        }

        /* Step 2 - ECG R-peak detection, RAW method only. detect_channel_raw
           returns on an empty signal, and an emptiness test would not mean
           much anyway: file_to_bin fills absent channels with placeholder
           vectors, so a non-empty ecgSignal2 does not prove CH2 exists. The
           zero-variance guard inside the detector leaves .raw empty for a
           constant placeholder; the authoritative presence test is
           CreateEcgTemplatesFast's (signal non-empty AND chN.raw non-empty). */
        detect_channel_raw(d.ch1, d.ecgSignal, cfg.ecg_upsample_rate, fileID, cfg, ecg1_inverted);
        detect_channel_raw(d.ch2, d.ecgSignal2, cfg.ecg_upsample_rate, fileID, cfg, ecg2_inverted);
        detect_channel_raw(d.ch3, d.ecgSignal3, cfg.ecg_upsample_rate, fileID, cfg, ecg3_inverted);

    }

    return data;
}


inline void augment_ecg_ppg_pairs_sqabs(std::vector<output_binfile_data>& data, std::string fileID, double ecgRate, config_entry cfg)
{
    //runs peakfinding on the squared and abs val lines
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < static_cast<int>(data.size()); ++i) {
        auto& d = data[i];
        detect_channel_sqabs(d.ch1, d.ecgSignal, ecgRate, fileID, cfg);
        detect_channel_sqabs(d.ch2, d.ecgSignal2, ecgRate, fileID, cfg);
        detect_channel_sqabs(d.ch3, d.ecgSignal3, ecgRate, fileID, cfg);
    }
}