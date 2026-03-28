/**
 * @file   create_ecg_ppg_pairs.hpp
 * @brief  Identify R-peaks on ECG channels using three preprocessing methods
 *         (raw, squared, abs-value), pair channel-1 raw R-peaks with PPG valleys.
 *         Original signals are always stored unmodified in the output.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-26
 */

#pragma once

#include <vector>
#include <string>
#include <cmath>
#include "binfile_handling.h"
#include "SegmentPPG.h"
#include "JoinedRR.h"
#include "pairRtoPPGBeat.h"
#include "StatsUtils.h"

 /**
  * @brief  Run R-peak detection on three versions of a signal: raw,
  *         squared, and absolute value. Each is independently noise-checked.
  *
  * @param[in]  signal    Original ECG channel.
  * @param[in]  fileID    Study identifier forwarded to JoinedRR.
  * @param[in]  hasPPG    Whether PPG data is available for cross-validation.
  * @param[in]  ppgCount  Number of PPG valleys (ignored when hasPPG is false).
  * @return     ChannelRPeaks with results from all three methods.
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
                rIndex = JoinedRR(sig, fileID);

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
    std::vector<double> sq(signal.size());
    for (std::size_t i = 0; i < signal.size(); ++i) {
        sq[i] = signal[i] * signal[i];
    }
    run_detection(sq, result.squared, result.squared_noisy);

    /* Method 3: absolute value signal */
    std::vector<double> ab(signal.size());
    for (std::size_t i = 0; i < signal.size(); ++i) {
        ab[i] = std::fabs(signal[i]);
    }
    run_detection(ab, result.absval, result.absval_noisy);

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

        bool hasPPG = !seg.ppg_signal.empty();
        bool hasEcg2 = !seg.ecg_signal_2.empty();
        bool hasEcg3 = !seg.ecg_signal_3.empty();

        /* ----------------------------------------------------------------
         * Step 1 — PPG pulse segmentation
         * ---------------------------------------------------------------- */
        if (hasPPG) {
            try {
                SegmentPPGResult ppgResult = SegmentPPG(seg.ppg_signal);
                d.ppgMinAmps = std::move(ppgResult.ppgMinAmps);
                d.ppgMaxAmps = std::move(ppgResult.maxAmps);
                d.bad_segment = false;
            }
            catch (...) {
                d.ppgMinAmps.clear();
                d.ppgMaxAmps.clear();
                d.bad_segment = true;
            }
        }
        else {
            d.bad_segment = false;
        }

        /* ----------------------------------------------------------------
         * Step 2 — R-peak detection (3 methods × up to 3 channels)
         * ---------------------------------------------------------------- */
        if (use_R_algorithm) {
            d.ch1 = detect_channel_3way(seg.ecg_signal_1, fileID, hasPPG,
                d.ppgMinAmps.size());

            if (hasEcg2) {
                d.ch2 = detect_channel_3way(seg.ecg_signal_2, fileID, hasPPG,
                    d.ppgMinAmps.size());
            }

            if (hasEcg3) {
                d.ch3 = detect_channel_3way(seg.ecg_signal_3, fileID, hasPPG,
                    d.ppgMinAmps.size());
            }
        }

        /* Free detection-phase memory early */
        seg.ecg_signal_2.clear();
        seg.ecg_signal_2.shrink_to_fit();
        seg.ecg_signal_3.clear();
        seg.ecg_signal_3.shrink_to_fit();
        seg.sleep_state_signal.clear();
        seg.sleep_state_signal.shrink_to_fit();

        /* ----------------------------------------------------------------
         * Step 3 — Pair channel-1 raw R-peaks with PPG valleys
         * ---------------------------------------------------------------- */
        bool rIsNoise = d.ch1.raw_noisy;

        if (!rIsNoise && !d.ch1.raw.empty() && hasPPG) {
            try {
                d.pairs = pairRtoPPGBeat(
                    seg.ecg_signal_1, seg.ppg_signal,
                    d.ch1.raw, d.ppgMinAmps);
            }
            catch (...) {
                if (!d.bad_segment) {
                    d.ch1.raw.clear();
                    d.pairs = build_unpaired(d.ppgMinAmps);
                }
                else {
                    d.pairs.clear();
                }
            }
        }
        else if (!rIsNoise && !d.ch1.raw.empty() && !hasPPG) {
            d.pairs.clear();
        }
        else {
            if (!d.bad_segment && hasPPG) {
                d.ch1.raw.clear();
                d.pairs = build_unpaired(d.ppgMinAmps);
            }
            else {
                d.ch1.raw.clear();
                d.pairs.clear();
            }
        }
    }

    return data;
}