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
#include "template_marking_gui/NormalizeFeatures.hpp"

struct PPGTemplatesResult {
    vector<vector<double>> templates;   // [bin][sample]
    vector<vector<double>> iqrs;        // [bin][sample], same shape as templates
    vector<vector<vector<double>>> kept; // [bin][beat][sample] retained snips
    vector<int> peakCol;                // [bin] systolic peak column (R1..R2)
    vector<int> footCol;                // [bin] foot column (R1..peak)
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
    std::vector<double>& outIqr,
    std::vector<std::vector<double>>& outKeptBeats,
    int& outPeakCol,
    int& outFootCol)
{
    outTemplate.clear();
    outIqr.clear();
    outKeptBeats.clear();
    outPeakCol = -1;
    outFootCol = -1;

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

    // ---- Deterministic PPG fiducials from the real R-peaks -------------
    // The template is R-anchored: R1 lands at column padSeconds*channelRate
    // (0.25 s in) by construction. R2 = R1 + one RR interval. We compute the
    // systolic peak as the max in [R1, R2] (exactly one pulse -> no risk of
    // grabbing a later pulse) and the foot as the min in [R1, peak]. Using the
    // true R-pair interval (not a fixed window) makes these exact.
    // Computed BEFORE the spread below, since the spread needs outFootCol.
    {
        const int N = static_cast<int>(outTemplate.size());
        const int r1 = std::clamp(
            static_cast<int>(std::llround(padSeconds * channelRate)), 0,
            std::max(0, N - 1));
        // Median RR in ECG samples -> channel samples.
        std::vector<double> gaps;
        gaps.reserve(masterPeaksEcg.size());
        for (size_t k = 1; k < masterPeaksEcg.size(); ++k)
            gaps.push_back(static_cast<double>(masterPeaksEcg[k] - masterPeaksEcg[k - 1]));
        int rrCh = 0;
        if (!gaps.empty()) {
            std::sort(gaps.begin(), gaps.end());
            const double medGapEcg = gaps[gaps.size() / 2];
            rrCh = static_cast<int>(std::llround(medGapEcg * scale));
        }
        const int r2 = (rrCh > 0) ? std::min(N - 1, r1 + rrCh) : (N - 1);

        if (N > 0 && r2 > r1) {
            int pk = r1; double pmax = -std::numeric_limits<double>::infinity();
            for (int i = r1; i <= r2; ++i)
                if (!std::isnan(outTemplate[i]) && outTemplate[i] > pmax) {
                    pmax = outTemplate[i]; pk = i;
                }
            int ft = r1; double fmin = std::numeric_limits<double>::infinity();
            for (int i = r1; i <= pk; ++i)
                if (!std::isnan(outTemplate[i]) && outTemplate[i] < fmin) {
                    fmin = outTemplate[i]; ft = i;
                }
            outPeakCol = pk;
            outFootCol = ft;
        }
    }

    // Per-sample robust spread, in the SAME units the mean trace will
    // eventually be displayed/exported in (minus only the final /ref
    // division, which normalize_features::scale_array_by_ref applies at
    // display time -- never here). Each beat is first converted to its own
    // local perfusion-index ratio using ITS OWN foot (outFootCol), per the
    // documented algorithm -- never a median/global foot -- then the
    // cross-beat IQR is taken of those local-ratio values. This is NOT the
    // same as taking the IQR of raw amplitudes, because the per-sample
    // transform's slope differs beat-to-beat (each beat has its own foot).
    outIqr = (outFootCol >= 0)
        ? normalize_features::local_ratio_iqr(aligned.beats, outFootCol)
        : std::vector<double>(maxLen, 0.0);

    // Retain aligned per-beat slices for downstream (snips CSV, etc).
    outKeptBeats.reserve(aligned.beats.size());
    for (const auto& sl : aligned.beats) outKeptBeats.push_back(sl);

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
    out.iqrs.assign(n, {});
    out.kept.assign(n, {});
    out.peakCol.assign(n, -1);
    out.footCol.assign(n, -1);

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
                out.templates[i], out.iqrs[i], out.kept[i],
                out.peakCol[i], out.footCol[i]);
        }
        catch (...) {
            out.templates[i] = {};
            out.iqrs[i] = {};
            out.kept[i] = {};
            out.peakCol[i] = -1;
            out.footCol[i] = -1;
        }
    }

    return out;
}