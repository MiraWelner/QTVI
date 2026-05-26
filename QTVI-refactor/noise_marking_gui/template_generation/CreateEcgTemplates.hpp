/**
 * @file   CreateEcgTemplates.hpp
 * @brief  Create ECG templates for each bin using EnsembleTemplate.
 *         Builds templates from 3 channels x 4 preprocessing methods
 *         (raw, squared, absval, unfiltered).
 *
 *         The "unfiltered" method uses the original ECG signal (ecgSignal,
 *         ecgSignal2, ecgSignal3) with the raw R-peaks to build a template
 *         from the signal before any preprocessing or filtering.
 *
 *         For the "raw" method on ch1 only, we additionally capture the
 *         surviving aligned beats that contributed to the template, so
 *         downstream code can write them out for QC visualization.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-27
 */
#pragma once

#include "TemplateTypes.hpp"
#include "EnsembleTemplate.hpp"
#include <atomic>

#ifdef _OPENMP
#include <omp.h>
#endif

struct SingleMethodResult {
    vector<double> ecgTemplate;
    double ppg_alignment_point;
    double avg_r_expand;
};

static inline SingleMethodResult build_ecg_template_for_method(
    const vector<double>& ecgSignal,
    const vector<size_t>& rpeaks,
    const vector<vector<double>>& pairs,
    double std_multiplier,
    vector<vector<double>>* out_kept_beats = nullptr)
{
    SingleMethodResult res;
    res.ecgTemplate = {};
    res.ppg_alignment_point = NaN;
    res.avg_r_expand = 0.0;

    if (rpeaks.size() < 2 || ecgSignal.empty()) return res;

    // MATLAB Daniel: reduce ECG to the first 10% of the bin and drop
    // R-peaks at/after (reduced_size - 1). The trailing -1 is a 1-sample
    // safety margin so the very last beat's window can't run off the end.
    const size_t reduced_size =
        static_cast<size_t>(std::round(static_cast<double>(ecgSignal.size()) / 10.0));
    if (reduced_size < 2) return res;

    vector<double> ecg_reduced(ecgSignal.begin(),
        ecgSignal.begin() + reduced_size);

    vector<size_t> r;
    r.reserve(rpeaks.size());
    for (auto idx : rpeaks) {
        // MATLAB: r = ecgRIndex(ecgRIndex < reduced_size - 1)
        if (idx + 1 < reduced_size) r.push_back(idx);
    }
    if (r.size() < 2) return res;

    vector<size_t> lens;
    lens.reserve(r.size() - 1);
    for (size_t j = 0; j + 1 < r.size(); ++j)
        lens.push_back(static_cast<size_t>((r[j + 1] - r[j]) / 5));

    {
        vector<double> ld(lens.begin(), lens.end());
        res.avg_r_expand = median(ld);
    }

    try {
        res.ecgTemplate = EnsembleTemplate(
            ecg_reduced, r, std_multiplier, "ecg", lens, out_kept_beats);
    }
    catch (...) {
        res.ecgTemplate = {};
        if (out_kept_beats) out_kept_beats->clear();
        return res;
    }

    if (!pairs.empty()) {
        vector<double> diffs;
        for (const auto& p : pairs) {
            if (p.size() >= 2 && p[0] >= 0 && p[1] >= 0)
                diffs.push_back(p[0] - p[1]);
        }
        if (!diffs.empty()) res.ppg_alignment_point = median(diffs);
    }

    return res;
}

static inline void init_channel_result(EcgChannelResult& cr, size_t n) {
    cr.ecgTemplates_raw.resize(n);
    cr.ecgTemplates_squared.resize(n);
    cr.ecgTemplates_absval.resize(n);
    cr.ecgTemplates_unfiltered.resize(n);

    cr.ppg_alignment_point_raw.resize(n, NaN);
    cr.ppg_alignment_point_squared.resize(n, NaN);
    cr.ppg_alignment_point_absval.resize(n, NaN);
    cr.ppg_alignment_point_unfiltered.resize(n, NaN);

    cr.avg_r_expand_raw.resize(n, 0.0);
    cr.avg_r_expand_squared.resize(n, 0.0);
    cr.avg_r_expand_absval.resize(n, 0.0);
    cr.avg_r_expand_unfiltered.resize(n, 0.0);

    cr.kept_beats_raw.resize(n);
}

/**
 * @brief  Process all 4 methods for one channel in one bin.
 *
 * @param cr             Channel result accumulator
 * @param bins           All bins (for pairs access)
 * @param i              Current bin index
 * @param ecgSignal      The signal used for raw/sq/abs methods (may be preprocessed)
 * @param origSignal     The original unfiltered ECG signal for this channel
 * @param ch             Channel R-peaks struct
 * @param std_multiplier Outlier threshold multiplier
 * @param capture_raw_beats  If true, capture the surviving aligned beats
 *                           from the "raw" method into cr.kept_beats_raw[i].
 */
static inline void process_channel(
    EcgChannelResult& cr,
    const vector<output_binfile_data>& bins,
    size_t i,
    const vector<double>& ecgSignal,
    const vector<double>& origSignal,
    const ChannelRPeaks& ch,
    double std_multiplier,
    bool capture_raw_beats = false)
{
    const auto& bin = bins[i];

    // Method 1: raw (detection signal + raw R-peaks)
    vector<vector<double>>* capture =
        (capture_raw_beats && i < cr.kept_beats_raw.size())
        ? &cr.kept_beats_raw[i] : nullptr;
    auto raw_res = build_ecg_template_for_method(
        ecgSignal, ch.raw, bin.pairs, std_multiplier, capture);
    cr.ecgTemplates_raw[i] = raw_res.ecgTemplate;
    cr.ppg_alignment_point_raw[i] = raw_res.ppg_alignment_point;
    cr.avg_r_expand_raw[i] = raw_res.avg_r_expand;

    // Method 2: squared (squared signal + squared R-peaks)
    const auto& sq_sig = ch.squared_signal.empty() ? ecgSignal : ch.squared_signal;
    auto sq_res = build_ecg_template_for_method(sq_sig, ch.squared, bin.pairs, std_multiplier);
    cr.ecgTemplates_squared[i] = sq_res.ecgTemplate;
    cr.ppg_alignment_point_squared[i] = sq_res.ppg_alignment_point;
    cr.avg_r_expand_squared[i] = sq_res.avg_r_expand;

    // Method 3: absval (abs-value signal + absval R-peaks)
    const auto& abs_sig = ch.absval_signal.empty() ? ecgSignal : ch.absval_signal;
    auto abs_res = build_ecg_template_for_method(abs_sig, ch.absval, bin.pairs, std_multiplier);
    cr.ecgTemplates_absval[i] = abs_res.ecgTemplate;
    cr.ppg_alignment_point_absval[i] = abs_res.ppg_alignment_point;
    cr.avg_r_expand_absval[i] = abs_res.avg_r_expand;

    // Method 4: unfiltered (original ECG signal + raw R-peaks)
    auto unfilt_res = build_ecg_template_for_method(origSignal, ch.raw, bin.pairs, std_multiplier);
    cr.ecgTemplates_unfiltered[i] = unfilt_res.ecgTemplate;
    cr.ppg_alignment_point_unfiltered[i] = unfilt_res.ppg_alignment_point;
    cr.avg_r_expand_unfiltered[i] = unfilt_res.avg_r_expand;
}

inline EcgTemplateResult CreateEcgTemplates(
    const vector<output_binfile_data>& bins,
    double std_multiplier)
{
    size_t n = bins.size();
    EcgTemplateResult res;
    init_channel_result(res.ch1, n);
    init_channel_result(res.ch2, n);
    init_channel_result(res.ch3, n);

    std::atomic<int> done{ 0 };

    int max_threads = std::min(8, (int)n);
#pragma omp parallel for schedule(dynamic) num_threads(max_threads)
    for (int i = 0; i < static_cast<int>(n); ++i) {

        const auto& bin = bins[i];

        // Ch1: capture the surviving raw-method beats for QC output
        process_channel(res.ch1, bins, i, bin.ecgSignal, bin.ecgSignal, bin.ch1,
            std_multiplier, /*capture_raw_beats=*/true);

        // Ch2: ecgSignal2 is both the detection signal and the original
        if (!bin.ecgSignal2.empty())
            process_channel(res.ch2, bins, i, bin.ecgSignal2, bin.ecgSignal2, bin.ch2,
                std_multiplier);

        // Ch3: ecgSignal3 is both the detection signal and the original
        if (!bin.ecgSignal3.empty())
            process_channel(res.ch3, bins, i, bin.ecgSignal3, bin.ecgSignal3, bin.ch3,
                std_multiplier);

        ++done;
    }
    return res;
}