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
  // Shared single-signal R-peak detection (was the run_detection lambda
  // inside detect_channel_3way). Factored out so the raw-only and
  // squared/absval halves can be invoked independently.
static inline void run_rpeak_detection(
    const std::vector<double>& sig,
    double ecgRate,
    const std::string& fileID,
    bool hasPPG,
    std::size_t ppgCount,
    std::vector<std::size_t>& rIndex,
    bool& noisy)
{
    noisy = false;
    rIndex.clear();

    if (std_dev(sig) == 0.0) {
        noisy = true;
        return;
    }

    try {
        JoinedRRResult jrr = JoinedRR_full(sig, ecgRate, fileID);
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
}

// FAST half: raw-method detection only. Fills result.raw / raw_noisy.
// This is everything the PPG pairing and the viewer-displayed templates
// depend on.
static inline void detect_channel_raw(
    ChannelRPeaks& result,
    const std::vector<double>& signal,
    double ecgRate,
    const std::string& fileID,
    bool hasPPG,
    std::size_t ppgCount)
{
    if (signal.empty()) return;
    run_rpeak_detection(signal, ecgRate, fileID, hasPPG, ppgCount,
        result.raw, result.raw_noisy);
}

// SLOW half: squared + absolute-value preprocessing and detection.
// Fills result.squared(_signal)/absval(_signal) and their noise flags.
// Deferred off the critical path -- nothing the viewer reads depends on
// these.
static inline void detect_channel_sqabs(
    ChannelRPeaks& result,
    const std::vector<double>& signal,
    double ecgRate,
    const std::string& fileID,
    bool hasPPG,
    std::size_t ppgCount)
{
    if (signal.empty()) return;

    /* Method 2: squared signal */
    result.squared_signal.resize(signal.size());
    for (std::size_t i = 0; i < signal.size(); ++i) {
        result.squared_signal[i] = signal[i] * signal[i];
    }
    run_rpeak_detection(result.squared_signal, ecgRate, fileID, hasPPG, ppgCount,
        result.squared, result.squared_noisy);

    /* Method 3: absolute value signal */
    result.absval_signal.resize(signal.size());
    for (std::size_t i = 0; i < signal.size(); ++i) {
        result.absval_signal[i] = std::fabs(signal[i]);
    }
    run_rpeak_detection(result.absval_signal, ecgRate, fileID, hasPPG, ppgCount,
        result.absval, result.absval_noisy);
}

// Original 3-method entry point, preserved by composition (raw then
// squared/absval) so existing callers are unchanged.
static inline ChannelRPeaks detect_channel_3way(
    const std::vector<double>& signal,
    double ecgRate,
    const std::string& fileID,
    bool hasPPG,
    std::size_t ppgCount)
{
    ChannelRPeaks result;
    if (signal.empty()) return result;
    detect_channel_raw(result, signal, ecgRate, fileID, hasPPG, ppgCount);
    detect_channel_sqabs(result, signal, ecgRate, fileID, hasPPG, ppgCount);
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
 * @param[in] use_R_algorithm   When false, R-peak detection is skipped entirely.
 * @param[in] fileID            Study identifier forwarded to sub-algorithms.
 * @param[in] ecgRate           ECG sample rate in Hz.
 * @param[in] ppgRate           PPG sample rate in Hz.
 * @return    One output_binfile_data per segment.
 */
inline std::vector<output_binfile_data> create_ecg_ppg_pairs_raw( std::vector<AnnealedSegment> annealedSegments, bool use_R_algorithm, std::string fileID, double ecgRate, double ppgRate){
    std::vector<output_binfile_data> data(annealedSegments.size());

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < static_cast<int>(annealedSegments.size()); ++i) {

        AnnealedSegment seg = std::move(annealedSegments[i]);
        auto& d = data[i];

        d.index = static_cast<int>(i);

        /* Store original signals unmodified. Kept in-memory because the
           squared/absval augment below (which may run later, on another
           thread) re-derives its preprocessed signals from these. */
        d.ppgSignal = seg.ppg_signal;
        d.ecgSignal = seg.ecg_signal_1;
        d.ecgSignal2 = seg.ecg_signal_2;
        d.ecgSignal3 = seg.ecg_signal_3;
        d.ppg_bin_indexs = std::move(seg.ppg_bin_indexs);
        d.ecg_bin_indexs = std::move(seg.ecg_bin_indexs);

        d.all_upsampled = std::move(seg.all_upsampled);
        d.all_raw_pairs_flat = std::move(seg.all_raw_pairs_flat);

        bool hasPPG = !d.ppgSignal.empty();
        bool hasEcg2 = !d.ecgSignal2.empty();
        bool hasEcg3 = !d.ecgSignal3.empty();

        /* Step 1 - PPG pulse segmentation */
        if (hasPPG) {
            try {
                SegmentPPGResult ppgResult = SegmentPPG(d.ppgSignal, ppgRate);
                d.ppgMinAmps = ppgResult.minAmps;
                d.ppgMaxAmps = ppgResult.maxAmps;
            }
            catch (...) {
                d.ppgMinAmps.clear();
                d.ppgMaxAmps.clear();
            }
        }

        /* Step 2 - ECG R-peak detection, RAW method only. */
        if (use_R_algorithm && !d.ecgSignal.empty()) {
            detect_channel_raw(d.ch1, d.ecgSignal, ecgRate, fileID, hasPPG, d.ppgMinAmps.size());
            if (hasEcg2)
                detect_channel_raw(d.ch2, d.ecgSignal2, ecgRate, fileID, hasPPG, d.ppgMinAmps.size());
            if (hasEcg3)
                detect_channel_raw(d.ch3, d.ecgSignal3, ecgRate, fileID, hasPPG, d.ppgMinAmps.size());
        }

        /* Step 3 - Pair channel-1 raw R-peaks with PPG valleys */
        d.bad_segment = false;
        if (!d.ch1.raw.empty() && !d.ppgMinAmps.empty()) {
            try {
                d.pairs = pairRtoPPGBeat(d.ecgSignal, d.ppgSignal, d.ch1.raw, d.ppgMinAmps, ecgRate, ppgRate);
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

/**
 * @brief  Fill the squared/absval R-peaks (and preprocessed signals) onto
 *         an existing peakResults produced by create_ecg_ppg_pairs_raw().
 *         Reads each channel's original ECG from d.ecgSignal/2/3, so the
 *         signals must still be present (they are, on both the in-memory
 *         and the re-hydrated-from-disk paths).
 *
 *         Touches only the squared/absval fields of each ChannelRPeaks --
 *         disjoint from everything the raw pass and the viewer read -- so
 *         this is safe to run on a worker thread while the raw results are
 *         consumed elsewhere, provided the raw pass has fully finished.
 */
inline void augment_ecg_ppg_pairs_sqabs(
    std::vector<output_binfile_data>& data, bool use_R_algorithm,
    std::string fileID, double ecgRate)
{
    if (!use_R_algorithm) return;

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < static_cast<int>(data.size()); ++i) {
        auto& d = data[i];
        bool hasPPG = !d.ppgSignal.empty();
        if (!d.ecgSignal.empty())
            detect_channel_sqabs(d.ch1, d.ecgSignal, ecgRate, fileID, hasPPG, d.ppgMinAmps.size());
        if (!d.ecgSignal2.empty())
            detect_channel_sqabs(d.ch2, d.ecgSignal2, ecgRate, fileID, hasPPG, d.ppgMinAmps.size());
        if (!d.ecgSignal3.empty())
            detect_channel_sqabs(d.ch3, d.ecgSignal3, ecgRate, fileID, hasPPG, d.ppgMinAmps.size());
    }
}

inline std::vector<output_binfile_data> create_ecg_ppg_pairs(std::vector<AnnealedSegment> annealedSegments, int dbg_plot, bool use_R_algorithm, std::string fileID,
    double ecgRate, double ppgRate)
{
    (void)dbg_plot;  // reserved, currently unused
    // Preserve the original all-methods behaviour by composition.
    std::vector<output_binfile_data> data =
        create_ecg_ppg_pairs_raw(std::move(annealedSegments), use_R_algorithm,
            fileID, ecgRate, ppgRate);
    augment_ecg_ppg_pairs_sqabs(data, use_R_algorithm, fileID, ecgRate);
    return data;
}