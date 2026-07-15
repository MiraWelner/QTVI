/**
 * @file   CreatePPGTemplates.hpp
 * @brief  Build PPG (and arterial-channel) templates by slicing per beat
 *         at the SAME real-time windows the ECG uses:
 *
 *             [t_R_i - pad, t_R_{i+1} + pad]
 *
 *         with pad = 0.25 s and ch1.raw R-peaks (ECG-frame samples) driving
 *         every channel's slicing. Each channel's slice bounds are converted
 *         from ECG samples to that channel's own samples via the rate ratio
 *         (channelRate / ecgRate). Column-wise NaN-skipping median across
 *         beats produces the template; the first R lands at column
 *         `pad * channelRate` (0.25 s in) by construction, matching the
 *         ECG's own R position in its own template.
 *
 *         Under this scheme the PPG template is no longer foot-anchored;
 *         it's R-anchored just like the ECG. No more find_foot_pulseox /
 *         AlignWaves / EnsembleTemplate for PPG -- the slicing math IS the
 *         alignment. ppgStartSample() should therefore be 0 in the viewer
 *         (real-time-aligned by construction). That viewer cleanup is
 *         Patch C.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <vector>
#include "template_marking_gui/alignment.hpp"
#include "TemplateTypes.hpp"

struct PPGTemplatesResult {
    vector<vector<double>> templates;   // [bin][sample]
    vector<vector<double>> stds;        // [bin][sample], same shape as templates
    vector<vector<vector<double>>> kept; // [bin][beat][sample] retained snips
};

/**
 * @brief  Shared "slice a channel by ch1.raw R-pairs" pulse-averager.
 *
 *         For every consecutive R-pair (R_i, R_{i+1}) in `masterPeaksEcg`
 *         (ECG-frame sample indices from bin.ch1.raw), compute the real-time
 *         window [t_R_i - pad, t_R_{i+1} + pad] and pull it out of `signal`
 *         (which lives in this channel's own sample space at `channelRate`).
 *         The R sample index maps into channel space as
 *         `round(r_ecg * channelRate / ecgRate)`.
 *
 *         Slices are aligned so that column 0 is `pad` seconds before the
 *         first R -- so R_first sits at column `padSamplesCh` in every
 *         slice. Variable-length RR => variable slice length; short slices
 *         contribute NaN past their real end. Column-wise NaN-skipping
 *         median => the template.
 *
 *         Reused by CreateArterialTemplates via a member-pointer for the
 *         signal.
 */
static inline void build_pulse_template_pair_windowed(
    const std::vector<double>& signal,
    double channelRate,
    const std::vector<size_t>& masterPeaksEcg,
    double ecgRate,
    double padSeconds,
    std::vector<double>& outTemplate,
    std::vector<double>& outStd,
    std::vector<std::vector<double>>& outKeptBeats)
{
    outTemplate.clear();
    outStd.clear();
    outKeptBeats.clear();

    if (signal.empty() || masterPeaksEcg.size() < 2 ||
        channelRate <= 0.0 || ecgRate <= 0.0) return;

    const double scale = channelRate / ecgRate;

    // Convert ch1.raw R-peaks from ECG samples to this channel's samples.
    std::vector<size_t> peaksCh;
    peaksCh.reserve(masterPeaksEcg.size());
    for (size_t r : masterPeaksEcg)
        peaksCh.push_back(static_cast<size_t>(std::llround(
            static_cast<double>(r) * scale)));

    // Per-bin peak-aligned + foot-vertical-aligned beat matrix.
    const auto aligned = alignment::extract_ppg_beats_and_align(signal, peaksCh);
    if (aligned.beats.empty()) return;

    const size_t maxLen = aligned.beats.front().size();

    // Column-wise NaN-skipping median => template.
    outTemplate.assign(maxLen, NaN);
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
        outTemplate[c] = (nc % 2 == 0)
            ? 0.5 * (col[nc / 2 - 1] + col[nc / 2])
            : col[nc / 2];
    }

    // Per-sample std (sample std, NaN skip, ddof=1).
    outStd.assign(maxLen, 0.0);
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
        outStd[c] = std::sqrt(ss / static_cast<double>(n - 1));
    }

    // Retain aligned per-beat slices for downstream (snips CSV, etc).
    outKeptBeats.reserve(aligned.beats.size());
    for (const auto& sl : aligned.beats) outKeptBeats.push_back(sl);

    // Silence unused warnings — the alignment function computes its own
    // slicing pad from RR intervals, so padSeconds isn't used here.
    (void)padSeconds; (void)channelRate; (void)ecgRate;
}

/**
 * @brief  PPG templates for every bin. Slicing is driven by bin.ch1.raw
 *         (ECG-frame R-peaks), consistent with the ECG templater.
 *
 * @param bins        Input bins.
 * @param ecgRate     ECG sample rate (for R-peak time base).
 * @param ppgRate     PPG sample rate.
 * @param padSeconds  0.25 (matches CreateEcgTemplates).
 */
inline PPGTemplatesResult CreatePPGTemplates(
    const vector<output_binfile_data>& bins,
    double ecgRate,
    double ppgRate,
    double padSeconds = 0.3)
{
    size_t n = bins.size();
    PPGTemplatesResult out;
    out.templates.assign(n, {});
    out.stds.assign(n, {});
    out.kept.assign(n, {});

    int ppg_threads = std::min(8, static_cast<int>(n > 0 ? n : 1));
#pragma omp parallel for schedule(dynamic) num_threads(ppg_threads)
    for (int i = 0; i < static_cast<int>(n); ++i) {
        const auto& b = bins[i];
        if (b.bad_segment || b.ppgSignal.empty() || b.ch1.raw.size() < 2)
            continue;
        try {
            build_pulse_template_pair_windowed(
                b.ppgSignal, ppgRate,
                b.ch1.raw, ecgRate,
                padSeconds,
                out.templates[i], out.stds[i], out.kept[i]);
        }
        catch (...) {
            out.templates[i] = {};
            out.stds[i] = {};
            out.kept[i] = {};
        }
    }

    return out;
}