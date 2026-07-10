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
 *         For the "raw" method we also capture two extra outputs:
 *           - The surviving aligned beats that contributed to the template
 *             (ch1 only -- downstream code writes these out for QC).
 *           - The per-sample std across those beats (all channels), which
 *             the viewer draws as a gray band around the displayed
 *             template. The other three methods (squared/absval/unfiltered)
 *             don't get std computed since the viewer never displays them.
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
    vector<double> ecgTemplate_std;   // empty for methods that don't compute std
    double ppg_alignment_point;
    double avg_r_expand;
};

static inline SingleMethodResult build_ecg_template_for_method(
    const vector<double>& ecgSignal,
    const vector<size_t>& rpeaks,
    const vector<vector<double>>& pairs,
    double std_multiplier,
    vector<vector<double>>* out_kept_beats = nullptr,
    bool compute_std = false)
{
    SingleMethodResult res;
    res.ecgTemplate = {};
    res.ecgTemplate_std = {};
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
    for (size_t j = 0; j + 1 < r.size(); ++j) {
        lens.push_back(static_cast<size_t>((r[j + 1] - r[j]) / 5));
    }
    res.avg_r_expand = median(vector<double>(lens.begin(), lens.end()));


    // Snip P-to-P instead of R-to-R. Cutting at each R put the *next* beat's
    // P at the tail; shifting every boundary back by kEcgPrePFrac*avg_r_expand
    // starts each window before its own P, so the averaged beat reads
    // P,QRS,T with the P at the front. Increased from 1.5 to 2.0 so slower
    // heart rates / longer PR intervals don't clip the P off the front.
    //
    // R therefore sits at index preP within the template. The viewer anchors
    // the PPG to R using this same fraction (see kEcgRAnchorFrac in
    // TemplateViewerWindow.cpp) -- the two MUST stay equal so ECG and PPG
    // remain time-aligned.
    constexpr double kEcgPrePFrac = 2.0;
    const size_t preP = static_cast<size_t>(std::llround(kEcgPrePFrac * res.avg_r_expand));
    vector<size_t> seg(r.size());
    for (size_t j = 0; j < r.size(); ++j)
        seg[j] = (r[j] > preP) ? r[j] - preP : 0;

    try {
        vector<double>* std_out = compute_std ? &res.ecgTemplate_std : nullptr;
        res.ecgTemplate = EnsembleTemplate(
            ecg_reduced, seg, std_multiplier, "ecg", {},   // P-to-P cut, no expand
            out_kept_beats, std_out);
    }
    catch (...) {
        res.ecgTemplate = {};
        res.ecgTemplate_std = {};
        if (out_kept_beats) out_kept_beats->clear();
        return res;
    }

    if (!pairs.empty()) {
        vector<double> diffs;
        for (const auto& p : pairs) {
            // p[0]=PPG foot, p[1]=R. Skip unpaired (<0) and step-8 fallback
            // rows (foot==R), which would inject spurious zero delays.
            if (p.size() >= 2 && p[0] >= 0 && p[1] >= 0 && p[0] != p[1])
                diffs.push_back(p[0] - p[1]);   // foot - R = transit delay (samples)
        }
        if (!diffs.empty()) res.ppg_alignment_point = median(diffs);
    }
    return res;
}

static inline void init_channel_result(EcgChannelResult& cr, size_t n) {
    cr.ecgTemplates_raw.resize(n);
    cr.ecgTemplates_raw_std.resize(n);
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
 *         Only the raw method computes std (the other three are never
 *         displayed in the viewer). Only ch1 captures the surviving raw
 *         beats for QC output.
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
 // FAST methods: raw (the displayed one, with std) + unfiltered. These are
 // everything the viewer renders. Captures the ch1 raw beats for QC.
static inline void process_channel_fast(
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

    // Method 1: raw (detection signal + raw R-peaks). Only method with std.
    vector<vector<double>>* capture =
        (capture_raw_beats && i < cr.kept_beats_raw.size())
        ? &cr.kept_beats_raw[i] : nullptr;
    auto raw_res = build_ecg_template_for_method(
        ecgSignal, ch.raw, bin.pairs, std_multiplier,
        capture, /*compute_std=*/true);
    cr.ecgTemplates_raw[i] = raw_res.ecgTemplate;
    cr.ecgTemplates_raw_std[i] = raw_res.ecgTemplate_std;
    cr.ppg_alignment_point_raw[i] = raw_res.ppg_alignment_point;
    cr.avg_r_expand_raw[i] = raw_res.avg_r_expand;

    // Method 4: unfiltered (original ECG signal + raw R-peaks). No std.
    auto unfilt_res = build_ecg_template_for_method(
        origSignal, ch.raw, bin.pairs, std_multiplier,
        nullptr, /*compute_std=*/false);
    cr.ecgTemplates_unfiltered[i] = unfilt_res.ecgTemplate;
    cr.ppg_alignment_point_unfiltered[i] = unfilt_res.ppg_alignment_point;
    cr.avg_r_expand_unfiltered[i] = unfilt_res.avg_r_expand;
}

// SLOW methods: squared + absval. Not displayed by the viewer; safe to
// compute off the critical path. Writes only the squared/absval fields of
// cr (which process_channel_fast leaves untouched).
static inline void process_channel_slow(
    EcgChannelResult& cr,
    const vector<output_binfile_data>& bins,
    size_t i,
    const vector<double>& ecgSignal,
    const ChannelRPeaks& ch,
    double std_multiplier)
{
    const auto& bin = bins[i];

    // Method 2: squared (squared signal + squared R-peaks). No std.
    const auto& sq_sig = ch.squared_signal.empty() ? ecgSignal : ch.squared_signal;
    auto sq_res = build_ecg_template_for_method(
        sq_sig, ch.squared, bin.pairs, std_multiplier,
        nullptr, /*compute_std=*/false);
    cr.ecgTemplates_squared[i] = sq_res.ecgTemplate;
    cr.ppg_alignment_point_squared[i] = sq_res.ppg_alignment_point;
    cr.avg_r_expand_squared[i] = sq_res.avg_r_expand;

    // Method 3: absval (abs-value signal + absval R-peaks). No std.
    const auto& abs_sig = ch.absval_signal.empty() ? ecgSignal : ch.absval_signal;
    auto abs_res = build_ecg_template_for_method(
        abs_sig, ch.absval, bin.pairs, std_multiplier,
        nullptr, /*compute_std=*/false);
    cr.ecgTemplates_absval[i] = abs_res.ecgTemplate;
    cr.ppg_alignment_point_absval[i] = abs_res.ppg_alignment_point;
    cr.avg_r_expand_absval[i] = abs_res.avg_r_expand;
}

// Original all-4-methods entry point, preserved by composition.
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
    process_channel_fast(cr, bins, i, ecgSignal, origSignal, ch,
        std_multiplier, capture_raw_beats);
    process_channel_slow(cr, bins, i, ecgSignal, ch, std_multiplier);
}

// FAST pass: raw + unfiltered templates (everything the viewer needs).
// Leaves the squared/absval vectors sized-but-empty for CreateEcgTemplatesSlow.
inline EcgTemplateResult CreateEcgTemplatesFast(
    const vector<output_binfile_data>& bins,
    double std_multiplier)
{
    size_t n = bins.size();
    EcgTemplateResult res;
    init_channel_result(res.ch1, n);
    init_channel_result(res.ch2, n);
    init_channel_result(res.ch3, n);

    int max_threads = std::min(8, (int)n);
#pragma omp parallel for schedule(dynamic) num_threads(max_threads)
    for (int i = 0; i < static_cast<int>(n); ++i) {
        const auto& bin = bins[i];
        process_channel_fast(res.ch1, bins, i, bin.ecgSignal, bin.ecgSignal, bin.ch1,
            std_multiplier, /*capture_raw_beats=*/true);
        if (!bin.ecgSignal2.empty())
            process_channel_fast(res.ch2, bins, i, bin.ecgSignal2, bin.ecgSignal2, bin.ch2,
                std_multiplier);
        if (!bin.ecgSignal3.empty())
            process_channel_fast(res.ch3, bins, i, bin.ecgSignal3, bin.ecgSignal3, bin.ch3,
                std_multiplier);
    }
    return res;
}

// SLOW pass: fills the squared/absval templates onto an EcgTemplateResult
// that has already been sized (e.g. by CreateEcgTemplatesFast, or by
// init_channel_result). Touches only squared/absval fields.
inline void CreateEcgTemplatesSlow(
    const vector<output_binfile_data>& bins,
    double std_multiplier,
    EcgTemplateResult& res)
{
    size_t n = bins.size();
    int max_threads = std::min(8, (int)n);
#pragma omp parallel for schedule(dynamic) num_threads(max_threads)
    for (int i = 0; i < static_cast<int>(n); ++i) {
        const auto& bin = bins[i];
        process_channel_slow(res.ch1, bins, i, bin.ecgSignal, bin.ch1, std_multiplier);
        if (!bin.ecgSignal2.empty())
            process_channel_slow(res.ch2, bins, i, bin.ecgSignal2, bin.ch2, std_multiplier);
        if (!bin.ecgSignal3.empty())
            process_channel_slow(res.ch3, bins, i, bin.ecgSignal3, bin.ch3, std_multiplier);
    }
}

inline EcgTemplateResult CreateEcgTemplates(
    const vector<output_binfile_data>& bins,
    double std_multiplier)
{
    EcgTemplateResult res = CreateEcgTemplatesFast(bins, std_multiplier);
    CreateEcgTemplatesSlow(bins, std_multiplier, res);
    return res;
}