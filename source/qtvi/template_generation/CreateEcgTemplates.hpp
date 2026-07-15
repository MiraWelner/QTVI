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
#include "template_marking_gui\alignment.hpp"
#include <atomic>

#ifdef _OPENMP
#include <omp.h>
#endif

struct SingleMethodResult {
    vector<double> ecgTemplate;
    vector<double> ecgTemplate_std;   // empty for methods that don't compute std
    double ppg_alignment_point;
    double avg_r_expand;
    // Number of slices that survived the drop rules (r1<=r0, len<3) and
    // were fed to the column-wise median. Under Patch B this is nearly
    // always n_rpeaks - 1, but any future length/amplitude filter would
    // drop it lower; the widget shows this so bins with poor signal are
    // visible at a glance.
    size_t n_beats = 0;
};

static inline SingleMethodResult build_ecg_template_for_method(
    const vector<double>& ecgSignal,
    const vector<size_t>& rpeaks,
    const vector<vector<double>>& pairs,
    double /*std_multiplier*/,
    double ecgRate,
    vector<vector<double>>* out_kept_beats = nullptr,
    bool compute_std = false)
{
    // NEW SLICING (Phase A of the pair-window refactor):
    //   For every consecutive R-peak pair (R_i, R_{i+1}) in `rpeaks`,
    //   slice ecgSignal at [R_i - 1s, R_{i+1} + 1s] and stack the slices.
    //   All slices start 1s before R_first, so R_first sits at column
    //   `ecgRate` (== 1s worth of samples) in every slice, and therefore
    //   in the column-wise-median template.
    //
    //   Slices have variable length (RR varies beat-to-beat). Shorter
    //   slices contribute NaN past their real end; the NaN-skipping
    //   column median downweights those columns automatically.
    //
    //   The old machinery (10% bin reduction, preP snipping via
    //   kEcgPrePFrac, EnsembleTemplate's foot alignment for PPG) is
    //   gone. The whole template is now deterministically R-anchored.
    //
    //   avg_r_expand is written as ecgRate / 2 so the viewer's
    //   preP = 2 * avg_r_expand math still points R at column
    //   `ecgRate` -- keeps Patch A compatible with the old viewer
    //   until Patch C simplifies it.
    SingleMethodResult res;
    res.ecgTemplate = {};
    res.ecgTemplate_std = {};
    res.ppg_alignment_point = NaN;
    res.avg_r_expand = 0.0;

    if (rpeaks.size() < 2 || ecgSignal.empty() || ecgRate <= 0.0) return res;

    // Alignment (see alignment.hpp): slice each beat as
    // [R_i - 0.25*RR_i, R_i + 0.75*RR_i], then do TWO in-place NaN-padded
    // shifts on a shared axis:
    //   Pass 1: R-align -- every beat's R lands at r_aligned_col
    //   Pass 2: Q-align -- every beat's Q lands at q_aligned_col
    // The returned beats are all the SAME WIDTH; NaN cells outside each
    // beat's real range don't participate in the column-wise median or std.
    const alignment::BeatSet aligned =
        alignment::extract_beats_and_mode(ecgSignal, rpeaks);
    if (aligned.beats.empty() || aligned.mode_length <= 0) return res;

    const size_t maxLen = aligned.beats.front().size();   // shared-axis width
    res.n_beats = aligned.beats.size();

    // avg_r_expand: back to a rate-derived constant (0.15 * ecgRate =
    // samples per fifth of a nominal RR cycle at the ECG rate). This is
    // used by downstream code that assumes a stable per-subject value --
    // PPG marker seeder, arterial seeder, etc. -- so it must NOT depend on
    // per-bin alignment output. Alignment's own R column lives on the
    // BeatSet (r_aligned_col); code that needs the aligned R column reads
    // that directly.
    res.avg_r_expand = 0.15 * ecgRate;

    // Column-wise NaN-skipping median over the aligned beats => template.
    res.ecgTemplate.assign(maxLen, NaN);
    for (size_t c = 0; c < maxLen; ++c) {
        std::vector<double> col;
        col.reserve(aligned.beats.size());
        for (const auto& sl : aligned.beats) {
            const double v = sl[c];
            if (!std::isnan(v)) col.push_back(v);
        }
        if (col.empty()) continue;
        std::sort(col.begin(), col.end());
        const size_t nc = col.size();
        res.ecgTemplate[c] = (nc % 2 == 0)
            ? 0.5 * (col[nc / 2 - 1] + col[nc / 2])
            : col[nc / 2];
    }

    // Optional per-sample std over the same aligned-beat matrix (sample std,
    // NaN skip, ddof=1). Used to draw the gray band under the raw template.
    if (compute_std) {
        res.ecgTemplate_std.assign(maxLen, 0.0);
        for (size_t c = 0; c < maxLen; ++c) {
            double sum = 0.0;
            size_t n = 0;
            for (const auto& sl : aligned.beats)
                if (!std::isnan(sl[c])) { sum += sl[c]; ++n; }
            if (n < 2) continue;
            const double mean = sum / static_cast<double>(n);
            double ss = 0.0;
            for (const auto& sl : aligned.beats)
                if (!std::isnan(sl[c])) {
                    const double d = sl[c] - mean;
                    ss += d * d;
                }
            res.ecgTemplate_std[c] = std::sqrt(ss / static_cast<double>(n - 1));
        }
    }

    // Retain the aligned per-beat slices for the snips CSV.
    if (out_kept_beats) {
        out_kept_beats->clear();
        out_kept_beats->reserve(aligned.beats.size());
        for (const auto& sl : aligned.beats) out_kept_beats->push_back(sl);
    }

    // PPG transit delay (median foot - R across paired beats), unchanged.
    if (!pairs.empty()) {
        std::vector<double> diffs;
        for (const auto& p : pairs) {
            if (p.size() >= 2 && p[0] >= 0 && p[1] >= 0 && p[0] != p[1])
                diffs.push_back(p[0] - p[1]);
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

    cr.n_beats_raw.resize(n, 0);

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
//
// Patch A change: ecgRate is threaded in for pair-window slicing, and
// `masterPeaks` (== bin.ch1.raw) drives the slicing for every channel so
// every channel's template covers the same real-time window.
static inline void process_channel_fast(
    EcgChannelResult& cr,
    const vector<output_binfile_data>& bins,
    size_t i,
    const vector<double>& ecgSignal,
    const vector<double>& origSignal,
    const vector<size_t>& masterPeaks,
    double std_multiplier,
    double ecgRate,
    bool capture_raw_beats = false)
{
    const auto& bin = bins[i];

    // Method 1: raw (detection signal + master R-peaks). Only method with std.
    vector<vector<double>>* capture =
        (capture_raw_beats && i < cr.kept_beats_raw.size())
        ? &cr.kept_beats_raw[i] : nullptr;
    auto raw_res = build_ecg_template_for_method(
        ecgSignal, masterPeaks, bin.pairs, std_multiplier, ecgRate,
        capture, /*compute_std=*/true);
    cr.ecgTemplates_raw[i] = raw_res.ecgTemplate;
    cr.ecgTemplates_raw_std[i] = raw_res.ecgTemplate_std;
    cr.ppg_alignment_point_raw[i] = raw_res.ppg_alignment_point;
    cr.avg_r_expand_raw[i] = raw_res.avg_r_expand;
    cr.n_beats_raw[i] = raw_res.n_beats;

    // Method 4: unfiltered (original ECG signal + master R-peaks). No std.
    auto unfilt_res = build_ecg_template_for_method(
        origSignal, masterPeaks, bin.pairs, std_multiplier, ecgRate,
        nullptr, /*compute_std=*/false);
    cr.ecgTemplates_unfiltered[i] = unfilt_res.ecgTemplate;
    cr.ppg_alignment_point_unfiltered[i] = unfilt_res.ppg_alignment_point;
    cr.avg_r_expand_unfiltered[i] = unfilt_res.avg_r_expand;
}

// SLOW methods: squared + absval. Not displayed by the viewer; safe to
// compute off the critical path. Writes only the squared/absval fields of
// cr (which process_channel_fast leaves untouched).
//
// Patch A change: same as fast pass, master ch1 peaks drive slicing. The
// per-channel preprocessed signals are still used (squared/absval variants
// have different amplitudes), but they're indexed at the master R-peak
// positions since R sample indices map 1:1 across ECG preprocessing.
static inline void process_channel_slow(
    EcgChannelResult& cr,
    const vector<output_binfile_data>& bins,
    size_t i,
    const vector<double>& ecgSignal,
    const ChannelRPeaks& ch,
    const vector<size_t>& masterPeaks,
    double std_multiplier,
    double ecgRate)
{
    const auto& bin = bins[i];

    // Method 2: squared (squared signal + master R-peaks). No std.
    const auto& sq_sig = ch.squared_signal.empty() ? ecgSignal : ch.squared_signal;
    auto sq_res = build_ecg_template_for_method(
        sq_sig, masterPeaks, bin.pairs, std_multiplier, ecgRate,
        nullptr, /*compute_std=*/false);
    cr.ecgTemplates_squared[i] = sq_res.ecgTemplate;
    cr.ppg_alignment_point_squared[i] = sq_res.ppg_alignment_point;
    cr.avg_r_expand_squared[i] = sq_res.avg_r_expand;

    // Method 3: absval (abs-value signal + master R-peaks). No std.
    const auto& abs_sig = ch.absval_signal.empty() ? ecgSignal : ch.absval_signal;
    auto abs_res = build_ecg_template_for_method(
        abs_sig, masterPeaks, bin.pairs, std_multiplier, ecgRate,
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
    const vector<size_t>& masterPeaks,
    double std_multiplier,
    double ecgRate,
    bool capture_raw_beats = false)
{
    process_channel_fast(cr, bins, i, ecgSignal, origSignal, masterPeaks,
        std_multiplier, ecgRate, capture_raw_beats);
    process_channel_slow(cr, bins, i, ecgSignal, ch, masterPeaks,
        std_multiplier, ecgRate);
}

// FAST pass: raw + unfiltered templates (everything the viewer needs).
// Leaves the squared/absval vectors sized-but-empty for CreateEcgTemplatesSlow.
inline EcgTemplateResult CreateEcgTemplatesFast(
    const vector<output_binfile_data>& bins,
    double std_multiplier,
    double ecgRate)
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
        const auto& master = bin.ch1.raw;   // ch1.raw drives every channel's slicing
        process_channel_fast(res.ch1, bins, i, bin.ecgSignal, bin.ecgSignal,
            master, std_multiplier, ecgRate, /*capture_raw_beats=*/true);
        // Only build ch2/ch3 templates when the channel is REAL: both the
        // signal is present AND R-peak detection actually found something.
        // file_to_bin fills absent channels with placeholder vectors (see
        // "0.0 = channel absent" in file_to_bin.hpp), so `.empty()` alone
        // isn't a reliable "channel exists" signal -- but a truly absent
        // channel will never produce R-peaks (zero-variance placeholder
        // trips run_rpeak_detection's std_dev==0 noisy flag). Requiring
        // bin.ch2.raw non-empty catches those.
        if (!bin.ecgSignal2.empty() && !bin.ch2.raw.empty())
            process_channel_fast(res.ch2, bins, i, bin.ecgSignal2, bin.ecgSignal2,
                master, std_multiplier, ecgRate, /*capture_raw_beats=*/true);
        if (!bin.ecgSignal3.empty() && !bin.ch3.raw.empty())
            process_channel_fast(res.ch3, bins, i, bin.ecgSignal3, bin.ecgSignal3,
                master, std_multiplier, ecgRate, /*capture_raw_beats=*/true);
    }
    return res;
}

// SLOW pass: fills the squared/absval templates onto an EcgTemplateResult
// that has already been sized (e.g. by CreateEcgTemplatesFast, or by
// init_channel_result). Touches only squared/absval fields.
inline void CreateEcgTemplatesSlow(
    const vector<output_binfile_data>& bins,
    double std_multiplier,
    double ecgRate,
    EcgTemplateResult& res)
{
    size_t n = bins.size();
    int max_threads = std::min(8, (int)n);
#pragma omp parallel for schedule(dynamic) num_threads(max_threads)
    for (int i = 0; i < static_cast<int>(n); ++i) {
        const auto& bin = bins[i];
        const auto& master = bin.ch1.raw;
        process_channel_slow(res.ch1, bins, i, bin.ecgSignal, bin.ch1, master,
            std_multiplier, ecgRate);
        // Same "real channel" gate as the fast pass (see comment there).
        if (!bin.ecgSignal2.empty() && !bin.ch2.raw.empty())
            process_channel_slow(res.ch2, bins, i, bin.ecgSignal2, bin.ch2, master,
                std_multiplier, ecgRate);
        if (!bin.ecgSignal3.empty() && !bin.ch3.raw.empty())
            process_channel_slow(res.ch3, bins, i, bin.ecgSignal3, bin.ch3, master,
                std_multiplier, ecgRate);
    }
}

inline EcgTemplateResult CreateEcgTemplates(
    const vector<output_binfile_data>& bins,
    double std_multiplier,
    double ecgRate)
{
    EcgTemplateResult res = CreateEcgTemplatesFast(bins, std_multiplier, ecgRate);
    CreateEcgTemplatesSlow(bins, std_multiplier, ecgRate, res);
    return res;
}