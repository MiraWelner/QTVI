/**
 * @file   create_ecg_ppg_pairs.hpp
 * @brief  Identify R-peaks on ECG channels using three preprocessing methods
 *         (raw, squared, abs-value), pair channel-1 raw R-peaks with PPG valleys.
 *         Original, preprocessed, and pre-bandpass signals are stored in the output.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-26
 */

#pragma once

#include <vector>
#include <string>
#include <cmath>
#include "peakfinding_io.hpp"
#include "SegmentPPG.hpp"
#include "JoinedRR.hpp"
#include "pairRtoPPGBeat.hpp"
#include "StatsUtils.hpp"

 /**
  * @brief  Run R-peak detection on three versions of a signal: raw,
  *         squared, and absolute value. Each is independently noise-checked.
  *         The squared and absval signals are stored in the result, along
  *         with the pre-bandpass (detrended) version of each.
  *
  * @param[in]  signal    Original ECG channel.
  * @param[in]  fileID    Study identifier forwarded to JoinedRR.
  * @param[in]  hasPPG    Whether PPG data is available for cross-validation.
  * @param[in]  ppgCount  Number of PPG valleys (ignored when hasPPG is false).
  * @return     ChannelRPeaks with results and signals from all three methods.
  */
static inline ChannelRPeaks detect_channel_3way(
    const std::vector<double>& signal,
    const std::string& fileID,
    bool hasPPG,
    std::size_t ppgCount)
{
    ChannelRPeaks result;

    if (signal.empty()) return result;

    auto run_detection = [&](const std::vector<double>& sig,
        std::vector<std::size_t>& rIndex,
        bool& noisy) {
            noisy = false;
            rIndex.clear();

            if (std_dev(sig) == 0.0) {
                noisy = true;
                return;
            }

            try {
                JoinedRRResult jrr = JoinedRR_full(sig, fileID);
                rIndex = std::move(jrr.peaks);

                if (hasPPG && ppgCount > 0) {
                    double r = static_cast<double>(rIndex.size());
                    double p = static_cast<double>(ppgCount);
                    if (r < p / 2.0 || p * 1.5 < r) {
                        noisy = true;
                    }
                }
            }
            catch (...) {
                rIndex.clear();
                noisy = true;
            }
        };

    /* Method 1: raw signal */
    run_detection(signal, result.raw, result.raw_noisy);

    /* Method 2: squared signal */
    result.squared_signal.resize(signal.size());
    for (std::size_t i = 0; i < signal.size(); ++i) {
        result.squared_signal[i] = signal[i] * signal[i];
    }
    run_detection(result.squared_signal, result.squared, result.squared_noisy);

    /* Method 3: absolute value signal */
    result.absval_signal.resize(signal.size());
    for (std::size_t i = 0; i < signal.size(); ++i) {
        result.absval_signal[i] = std::fabs(signal[i]);
    }
    run_detection(result.absval_signal, result.absval, result.absval_noisy);

    return result;
}

/**
 * @brief  Build placeholder pairs when no valid ECG R-peaks are available.
 *
 * @param[in]  ppgMinAmps  PPG valley indices to preserve.
 * @return     Pair matrix with [ppg_valley_idx, -1.0] per row.
 */
static inline std::vector<std::vector<double>> build_unpaired(
    const std::vector<std::size_t>& ppgMinAmps)
{
    std::vector<std::vector<double>> pairs;
    pairs.reserve(ppgMinAmps.size());
    for (std::size_t k = 0; k < ppgMinAmps.size(); ++k) {
        pairs.push_back({ static_cast<double>(ppgMinAmps[k]), -1.0 });
    }
    return pairs;
}

/**
 * @brief  Main processing pipeline: detect R-peaks using three methods on
 *         up to three ECG channels, segment PPG, and pair channel-1 raw
 *         R-peaks with PPG valleys.
 *
 * @param[in] annealedSegments  Per-segment ECG/PPG data (consumed via move).
 * @param[in] dbg_plot          Debug plotting flag (reserved, currently unused).
 * @param[in] use_R_algorithm   When false, R-peak detection is skipped entirely.
 * @param[in] fileID            Study identifier forwarded to sub-algorithms.
 * @return    One output_binfile_data per segment.
 */
inline std::vector<output_binfile_data> create_ecg_ppg_pairs(
    std::vector<AnnealedSegment> annealedSegments,
    int dbg_plot, bool use_R_algorithm, std::string fileID)
{
    std::vector<output_binfile_data> data(annealedSegments.size());

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < static_cast<int>(annealedSegments.size()); ++i) {

        AnnealedSegment seg = std::move(annealedSegments[i]);
        auto& d = data[i];

        d.index = static_cast<int>(i);

        /* Store original signals unmodified */
        d.ppgSignal = seg.ppg_signal;
        d.ecgSignal = seg.ecg_signal_1;
        d.ecgSignal2 = seg.ecg_signal_2;
        d.ecgSignal3 = seg.ecg_signal_3;
        d.ppg_bin_indexs = std::move(seg.ppg_bin_indexs);
        d.ecg_bin_indexs = std::move(seg.ecg_bin_indexs);

        /* Pass-through: route the full 41-channel set from input to output
           untouched. Move-semantics since the input segment is consumed
           in this loop iteration. The peakfinding algorithm below does
           not read or modify these vectors. */
        d.all_upsampled = std::move(seg.all_upsampled);
        d.all_raw_pairs_flat = std::move(seg.all_raw_pairs_flat);

        bool hasPPG = !seg.ppg_signal.empty();
        bool hasEcg2 = !seg.ecg_signal_2.empty();
        bool hasEcg3 = !seg.ecg_signal_3.empty();

        /* ----------------------------------------------------------------
         * Step 1 - PPG pulse segmentation
         * ---------------------------------------------------------------- */
        if (hasPPG) {
            try {
                SegmentPPGResult ppgResult = SegmentPPG(seg.ppg_signal);
                d.ppgMinAmps = ppgResult.minAmps;
                d.ppgMaxAmps = ppgResult.maxAmps;
            }
            catch (...) {
                d.ppgMinAmps.clear();
                d.ppgMaxAmps.clear();
            }
        }

        /* ----------------------------------------------------------------
         * Step 2 - ECG R-peak detection (3 methods per channel)
         *          Pre-bandpass signals are now captured automatically
         *          by detect_channel_3way via JoinedRR_full.
         * ---------------------------------------------------------------- */
        if (use_R_algorithm && !seg.ecg_signal_1.empty()) {
            d.ch1 = detect_channel_3way(seg.ecg_signal_1, fileID, hasPPG, d.ppgMinAmps.size());

            if (hasEcg2) {
                d.ch2 = detect_channel_3way(seg.ecg_signal_2, fileID, hasPPG, d.ppgMinAmps.size());
            }
            if (hasEcg3) {
                d.ch3 = detect_channel_3way(seg.ecg_signal_3, fileID, hasPPG, d.ppgMinAmps.size());
            }
        }

        /* ----------------------------------------------------------------
         * Step 3 - Pair channel-1 raw R-peaks with PPG valleys
         * ---------------------------------------------------------------- */
        d.bad_segment = false;
        if (!d.ch1.raw.empty() && !d.ppgMinAmps.empty()) {
            try {
                d.pairs = pairRtoPPGBeat(d.ecgSignal, d.ppgSignal, d.ch1.raw, d.ppgMinAmps);
            }
            catch (...) {
                d.pairs = build_unpaired(d.ppgMinAmps);
                d.bad_segment = true;
            }
        }
        else if (!d.ppgMinAmps.empty()) {
            d.pairs = build_unpaired(d.ppgMinAmps);
            d.bad_segment = true;
        }
    }

    return data;
}