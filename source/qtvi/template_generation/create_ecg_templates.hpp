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

#include "template_structs.hpp"
#include "template_marking_gui/alignment.hpp"
#include <atomic>
#include <fstream>
#include <string>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

 // Beat-move log destination. Set once from main/post_process before the build
 // (read-only afterwards, and the writer runs single-threaded post-loop, so no
 // race). Empty dir/stem => no log.
namespace ecg_move_log {
    inline std::string g_dir;
    inline std::string g_stem;
    inline void set(const std::string& dir, const std::string& stem) { g_dir = dir; g_stem = stem; }

    // Write one channel's per-bin/per-beat two-stage vertical shifts. Called
    // single-threaded (after the parallel build loop). `first` truncates +
    // writes the header; later channels append.
    inline void write_channel(const char* channel,
        const std::vector<std::vector<double>>& tp,
        const std::vector<std::vector<double>>& pq, bool first) {
        if (g_dir.empty() || g_stem.empty()) return;
        std::ofstream f(g_dir + "/" + g_stem + "_beat_moves.csv",
            first ? std::ios::trunc : std::ios::app);
        if (!f) return;
        if (first) f << "stem,channel,bin,beat,tp_mv_shift,pq_mv_shift\n";
        for (size_t b = 0; b < tp.size(); ++b)
            for (size_t k = 0; k < tp[b].size(); ++k) {
                const double pqv = (b < pq.size() && k < pq[b].size())
                    ? pq[b][k] : std::numeric_limits<double>::quiet_NaN();
                f << g_stem << ',' << channel << ',' << b << ',' << k
                    << ',' << tp[b][k] << ',' << pqv << '\n';
            }
    }
}

struct SingleMethodResult {
    vector<double> ecgTemplate;
    vector<double> ecg_template_iqr;   // empty for methods that don't compute std
    double ppg_alignment_point;
    int r_col = -1;   // true R column (alignment's r_aligned_col)
    int ref_beat_index = -1;
    size_t n_beats = 0;
    // Per-beat per-stage vertical DC shifts (TP stage, PQ stage) from the
    // two-stage leveling, surfaced for the move log.
    vector<double> tp_shift;
    vector<double> pq_shift;
};

static inline SingleMethodResult build_ecg_template_for_method(const vector<double>& ecgSignal, const vector<size_t>& rpeaks, const vector<vector<double>>& pairs,
    double ecgRate, vector<vector<double>>* out_kept_beats = nullptr, bool compute_iqr = false) {
    SingleMethodResult res;
    res.ecgTemplate = {};
    res.ecg_template_iqr = {};
    res.ppg_alignment_point = NaN;
    res.r_col = -1;

    if (rpeaks.size() < 2 || ecgSignal.empty() || ecgRate <= 0.0) return res;

    const alignment::ecg_beat_set aligned = alignment::extract_beats_and_align(ecgSignal, rpeaks, ecgRate);
    if (aligned.beats.empty() || aligned.median_length <= 0) return res;

    res.tp_shift = aligned.tp_shift;   // surface for the move log
    res.pq_shift = aligned.pq_shift;

    // Beats with baseline_source == NONE had neither a usable TP nor PQ
    // isoelectric reference, so their DC level is untrustworthy -- exclude
    // them from every amplitude-dependent aggregate below (median template,
    // std band, and the surviving-beats QC capture), same as the Tukey/
    // wave-score rejections that already ran upstream in extract_beats_and_
    // align(). NOTE: this changes prior behavior -- previously an
    // unavailable baseline meant "use the beat un-shifted anyway"; now it
    // means "exclude it entirely" per spec. If baseline_source is empty or
    // mismatched in length (e.g. ref_beat_index was invalid so Pass 3 never
    // ran), fall back to using every beat unfiltered rather than silently
    // producing an empty template.
    const bool haveSrc = aligned.baseline_source.size() == aligned.beats.size();
    std::vector<const std::vector<double>*> usable;
    usable.reserve(aligned.beats.size());
    for (size_t i = 0; i < aligned.beats.size(); ++i) {
        if (haveSrc && aligned.baseline_source[i] == alignment::BaselineSource::NONE) continue;
        usable.push_back(&aligned.beats[i]);
    }
    if (usable.empty()) {   // every beat's baseline was NONE -- fail safe, don't zero the template
        usable.reserve(aligned.beats.size());
        for (const auto& b : aligned.beats) usable.push_back(&b);
    }

    const size_t maxLen = usable.front()->size();   // shared-axis width
    res.n_beats = usable.size();

    // R column: the detected-R fiducial the template was built around, straight
    // from alignment (every beat's detected R lands at r_aligned_col). Passed
    // through as-is -- no re-detection (a window search would grab Q or S).
    res.r_col = aligned.r_aligned_col;
    res.ref_beat_index = aligned.ref_beat_index;
    // Column-wise NaN-skipping median over the aligned beats => template.
    res.ecgTemplate.assign(maxLen, NaN);
    for (size_t c = 0; c < maxLen; ++c) {
        std::vector<double> col;
        col.reserve(usable.size());
        for (const auto* sl : usable) {
            const double v = (*sl)[c];
            if (!std::isnan(v)) col.push_back(v);
        }
        if (col.empty()) continue;
        std::sort(col.begin(), col.end());
        const size_t nc = col.size();
        res.ecgTemplate[c] = (nc % 2 == 0)
            ? 0.5 * (col[nc / 2 - 1] + col[nc / 2])
            : col[nc / 2];
    }

    // Per-sample robust spread over the same aligned-beat matrix: per-sample
    // STD (ddof=1) -- changed from IQR per spec step 7. Used to draw the
    // gray band under the raw template.
    // NOTE: field/param names say "iqr" but this now computes the
    // per-sample STD (ddof=1), not the interquartile range -- changed
    // per spec step 7. Renaming ecg_template_iqr/compute_iqr throughout
    // the codebase (TemplateTypes.hpp, BinPlotWidget, TemplateBinIO,
    // template_io, ...) is a separate, larger follow-up; left as-is here
    // to keep this change to the computation only.
    if (compute_iqr) {
        res.ecg_template_iqr.assign(maxLen, 0.0);
        std::vector<double> col;
        col.reserve(usable.size());
        for (size_t c = 0; c < maxLen; ++c) {
            col.clear();
            for (const auto* sl : usable)
                if (!std::isnan((*sl)[c])) col.push_back((*sl)[c]);
            const size_t nc = col.size();
            if (nc < 2) continue;
            double mean = 0.0;
            for (double v : col) mean += v;
            mean /= static_cast<double>(nc);
            double sumsq = 0.0;
            for (double v : col) sumsq += (v - mean) * (v - mean);
            res.ecg_template_iqr[c] = std::sqrt(sumsq / static_cast<double>(nc - 1));   // ddof = 1
        }
    }

    // Retain the aligned per-beat slices for the snips CSV (baseline-
    // unreliable beats excluded, same set used for the template above).
    if (out_kept_beats) {
        out_kept_beats->clear();
        out_kept_beats->reserve(usable.size());
        for (const auto* sl : usable) out_kept_beats->push_back(*sl);
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
    cr.ecgTemplates_raw_iqr.resize(n);
    cr.ecgTemplates_squared.resize(n);
    cr.ecgTemplates_absval.resize(n);
    cr.ecgTemplates_unfiltered.resize(n);
    cr.ref_index_raw.resize(n, -1);

    cr.ppg_alignment_point_raw.resize(n, NaN);
    cr.ppg_alignment_point_squared.resize(n, NaN);
    cr.ppg_alignment_point_absval.resize(n, NaN);
    cr.ppg_alignment_point_unfiltered.resize(n, NaN);

    cr.r_col_raw.resize(n, -1);
    cr.r_col_squared.resize(n, -1);
    cr.r_col_absval.resize(n, -1);
    cr.r_col_unfiltered.resize(n, -1);

    cr.n_beats_raw.resize(n, 0);

    cr.kept_beats_raw.resize(n);
    cr.tp_shift_raw.resize(n);
    cr.pq_shift_raw.resize(n);
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
    double ecgRate,
    bool capture_raw_beats = false)
{
    const auto& bin = bins[i];

    // Method 1: raw (detection signal + master R-peaks). Only method with std.
    vector<vector<double>>* capture =
        (capture_raw_beats && i < cr.kept_beats_raw.size())
        ? &cr.kept_beats_raw[i] : nullptr;
    auto raw_res = build_ecg_template_for_method(
        ecgSignal, masterPeaks, bin.pairs, ecgRate,
        capture, /*compute_iqr=*/true);
    cr.ecgTemplates_raw[i] = raw_res.ecgTemplate;
    cr.ecgTemplates_raw_iqr[i] = raw_res.ecg_template_iqr;
    cr.ppg_alignment_point_raw[i] = raw_res.ppg_alignment_point;
    cr.r_col_raw[i] = raw_res.r_col;
    cr.n_beats_raw[i] = raw_res.n_beats;
    cr.ref_index_raw[i] = raw_res.ref_beat_index;
    if (i < cr.tp_shift_raw.size()) {
        cr.tp_shift_raw[i] = std::move(raw_res.tp_shift);   // distinct i -> race-free
        cr.pq_shift_raw[i] = std::move(raw_res.pq_shift);
    }

    // Method 4: unfiltered (original ECG signal + master R-peaks). No std.
    auto unfilt_res = build_ecg_template_for_method(
        origSignal, masterPeaks, bin.pairs, ecgRate,
        nullptr, /*compute_iqr=*/false);
    cr.ecgTemplates_unfiltered[i] = unfilt_res.ecgTemplate;
    cr.ppg_alignment_point_unfiltered[i] = unfilt_res.ppg_alignment_point;
    cr.r_col_unfiltered[i] = unfilt_res.r_col;
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
    double ecgRate)
{
    const auto& bin = bins[i];

    // Method 2: squared (squared signal + master R-peaks). No std.
    const auto& sq_sig = ch.squared_signal.empty() ? ecgSignal : ch.squared_signal;
    auto sq_res = build_ecg_template_for_method(
        sq_sig, masterPeaks, bin.pairs, ecgRate,
        nullptr, /*compute_iqr=*/false);
    cr.ecgTemplates_squared[i] = sq_res.ecgTemplate;
    cr.ppg_alignment_point_squared[i] = sq_res.ppg_alignment_point;
    cr.r_col_squared[i] = sq_res.r_col;

    // Method 3: absval (abs-value signal + master R-peaks). No std.
    const auto& abs_sig = ch.absval_signal.empty() ? ecgSignal : ch.absval_signal;
    auto abs_res = build_ecg_template_for_method(
        abs_sig, masterPeaks, bin.pairs, ecgRate,
        nullptr, /*compute_iqr=*/false);
    cr.ecgTemplates_absval[i] = abs_res.ecgTemplate;
    cr.ppg_alignment_point_absval[i] = abs_res.ppg_alignment_point;
    cr.r_col_absval[i] = abs_res.r_col;
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
    double ecgRate,
    bool capture_raw_beats = false)
{
    process_channel_fast(cr, bins, i, ecgSignal, origSignal, masterPeaks,
        ecgRate, capture_raw_beats);
    process_channel_slow(cr, bins, i, ecgSignal, ch, masterPeaks,
        ecgRate);
}

// FAST pass: raw + unfiltered templates (everything the viewer needs).
// Leaves the squared/absval vectors sized-but-empty for CreateEcgTemplatesSlow.
inline EcgTemplateResult CreateEcgTemplatesFast(
    const vector<output_binfile_data>& bins,
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
            master, ecgRate, /*capture_raw_beats=*/true);
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
                master, ecgRate, /*capture_raw_beats=*/true);
        if (!bin.ecgSignal3.empty() && !bin.ch3.raw.empty())
            process_channel_fast(res.ch3, bins, i, bin.ecgSignal3, bin.ecgSignal3,
                master, ecgRate, /*capture_raw_beats=*/true);
    }

    // Single-threaded, post-loop: write the per-beat vertical move log
    // (CH1 truncates+headers, CH2/CH3 append).
    ecg_move_log::write_channel("CH1", res.ch1.tp_shift_raw, res.ch1.pq_shift_raw, /*first=*/true);
    ecg_move_log::write_channel("CH2", res.ch2.tp_shift_raw, res.ch2.pq_shift_raw, /*first=*/false);
    ecg_move_log::write_channel("CH3", res.ch3.tp_shift_raw, res.ch3.pq_shift_raw, /*first=*/false);

    return res;
}

// SLOW pass: fills the squared/absval templates onto an EcgTemplateResult
// that has already been sized (e.g. by CreateEcgTemplatesFast, or by
// init_channel_result). Touches only squared/absval fields.
inline void CreateEcgTemplatesSlow(
    const vector<output_binfile_data>& bins,
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
            ecgRate);
        // Same "real channel" gate as the fast pass (see comment there).
        if (!bin.ecgSignal2.empty() && !bin.ch2.raw.empty())
            process_channel_slow(res.ch2, bins, i, bin.ecgSignal2, bin.ch2, master,
                ecgRate);
        if (!bin.ecgSignal3.empty() && !bin.ch3.raw.empty())
            process_channel_slow(res.ch3, bins, i, bin.ecgSignal3, bin.ch3, master,
                ecgRate);
    }
}

inline EcgTemplateResult CreateEcgTemplates(
    const vector<output_binfile_data>& bins,
    double ecgRate)
{
    EcgTemplateResult res = CreateEcgTemplatesFast(bins, ecgRate);
    CreateEcgTemplatesSlow(bins, ecgRate, res);
    return res;
}