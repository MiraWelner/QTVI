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
    const int padCh = static_cast<int>(std::llround(padSeconds * channelRate));
    const int signalN = static_cast<int>(signal.size());

    // Convert ch1.raw R-peaks from ECG samples to THIS channel's samples.
    std::vector<int> peaksCh;
    peaksCh.reserve(masterPeaksEcg.size());
    for (size_t r : masterPeaksEcg)
        peaksCh.push_back(static_cast<int>(std::llround(static_cast<double>(r) * scale)));

    struct Slice { std::vector<double> data; };
    std::vector<Slice> slices;
    slices.reserve(peaksCh.size() > 0 ? peaksCh.size() - 1 : 0);
    size_t maxLen = 0;

    for (size_t i = 0; i + 1 < peaksCh.size(); ++i) {
        const int r0 = peaksCh[i];
        const int r1 = peaksCh[i + 1];
        if (r1 <= r0) continue;

        const int startSig = r0 - padCh;
        const int endSig = r1;                 // no trailing pad -- matches ECG slicer
        const int len = endSig - startSig;
        if (len < 3) continue;

        std::vector<double> s(static_cast<size_t>(len), NaN);
        const int copyStart = std::max(0, startSig);
        const int copyEnd = std::min(signalN, endSig);
        for (int k = copyStart; k < copyEnd; ++k)
            s[static_cast<size_t>(k - startSig)] = signal[k];

        if (static_cast<size_t>(len) > maxLen) maxLen = static_cast<size_t>(len);
        slices.push_back({ std::move(s) });
    }
    if (slices.empty()) return;

    for (auto& sl : slices) sl.data.resize(maxLen, NaN);

    // Column-wise NaN-skipping median.
    outTemplate.assign(maxLen, NaN);
    for (size_t c = 0; c < maxLen; ++c) {
        std::vector<double> col;
        col.reserve(slices.size());
        for (const auto& sl : slices) {
            const double v = sl.data[c];
            if (!std::isnan(v)) col.push_back(v);
        }
        if (col.empty()) continue;
        std::sort(col.begin(), col.end());
        const size_t nc = col.size();
        outTemplate[c] = (nc % 2 == 0)
            ? 0.5 * (col[nc / 2 - 1] + col[nc / 2])
            : col[nc / 2];
    }

    // Per-sample std (ddof=1, NaN-skip). Empty columns => 0.
    outStd.assign(maxLen, 0.0);
    for (size_t c = 0; c < maxLen; ++c) {
        double sum = 0.0;
        size_t n = 0;
        for (const auto& sl : slices)
            if (!std::isnan(sl.data[c])) { sum += sl.data[c]; ++n; }
        if (n < 2) continue;
        const double mean = sum / static_cast<double>(n);
        double ss = 0.0;
        for (const auto& sl : slices)
            if (!std::isnan(sl.data[c])) {
                const double d = sl.data[c] - mean;
                ss += d * d;
            }
        outStd[c] = std::sqrt(ss / static_cast<double>(n - 1));
    }

    // Retain aligned per-beat slices.
    outKeptBeats.reserve(slices.size());
    for (auto& sl : slices) outKeptBeats.push_back(std::move(sl.data));
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

    // Pad templates to a common length with NaN, as before. Pad stds with 0
    // so the band collapses to the line in padded regions (no gray bulge at
    // the tail).
    size_t max_len = 0;
    for (const auto& t : out.templates)
        if (t.size() > max_len) max_len = t.size();

    for (auto& t : out.templates) t.resize(max_len, NaN);
    for (auto& s : out.stds)      s.resize(max_len, 0.0);

    return out;
}