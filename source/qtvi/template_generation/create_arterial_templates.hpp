/**
 * @file   create_arterial_templates.hpp
 * @brief  Two distinct pipelines, per spec:
 *
 *         PPG: R-anchored, same real-time windows the ECG uses --
 *             [t_R_i - pad, t_R_{i+1} + pad], driven by ch1.raw R-peaks
 *             (ECG-frame samples), converted to this channel's own samples
 *             via the rate ratio (channelRate / ecgRate). See
 *             CreatePulseTemplates / build_pulse_template_pair_windowed.
 *
 *         ABP / ART / ART_PULM: FOOT-anchored. Per spec, these are NOT
 *             sliced from borrowed ECG R-peaks -- systolic peaks are
 *             self-detected directly from each channel's own waveform, and
 *             find_foot_pulseox (the intersecting-tangent method) locates
 *             each beat's true foot within its own peak-to-peak segment.
 *             Every beat is re-sliced so its OWN foot lands at a fixed
 *             column. See CreateArterialTemplates /
 *             build_arterial_template_foot_anchored.
 *
 *         Both pipelines share the same matched-filter QC (pulse_matched_
 *         filter.hpp): the first 20 beats seed a provisional template, and
 *         every later beat is rejected once its normalized error against
 *         that template reaches 5%.
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
#include "pulse_matched_filter.hpp"
#include "find_foot_pulseox.hpp"

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
 *         Shared by CreatePulseTemplates for every channel (PPG, ABP, ART,
 *         ART_PULM) via a member-pointer for the signal.
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

    // ---- Matched-filter QC (spec: PPG/arterial pulse filter) -----------
    // Beats are already located (R-anchored above). The provisional
    // template is built from only the first 20 beats (or all of them, if
    // fewer than 20 are available); those seed beats are kept unconditionally.
    // Every remaining beat is then tested by NORMALIZED ERROR vs that
    // template (||beat - templ|| / ||templ||) and rejected once error
    // reaches 5%. Falls back to all beats if the filter would reject
    // everything (degenerate template).
    std::vector<std::vector<double>> filteredBeats;
    {
        const int w = static_cast<int>(aligned.beats.front().size());
        const size_t nSeed = std::min<size_t>(20, aligned.beats.size());
        const std::vector<std::vector<double>> seedBeats(
            aligned.beats.begin(), aligned.beats.begin() + nSeed);
        const std::vector<double> provisional =
            pulse_matched_filter::buildTemplate(seedBeats, w);

        filteredBeats.assign(aligned.beats.begin(), aligned.beats.begin() + nSeed);
        for (size_t i = nSeed; i < aligned.beats.size(); ++i) {
            const double err = pulse_matched_filter::normalizedError(aligned.beats[i], provisional);
            if (err < 0.05) filteredBeats.push_back(aligned.beats[i]);
        }
        if (filteredBeats.empty()) filteredBeats = aligned.beats;   // don't discard all
    }
    const auto& beatsForTemplate = filteredBeats;

    const size_t maxLen = beatsForTemplate.front().size();

    // Column-wise NaN-skipping median => template.
    outTemplate.assign(maxLen, NaN);
    for (size_t c = 0; c < maxLen; ++c) {
        std::vector<double> col;
        col.reserve(beatsForTemplate.size());
        for (const auto& sl : beatsForTemplate) {
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
        ? normalize_features::local_ratio_iqr(beatsForTemplate, outFootCol)
        : std::vector<double>(maxLen, 0.0);

    // Retain aligned per-beat slices for downstream (snips CSV, etc).
    outKeptBeats.reserve(beatsForTemplate.size());
    for (const auto& sl : beatsForTemplate) outKeptBeats.push_back(sl);

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
inline PPGTemplatesResult CreatePulseTemplates(
    const vector<output_binfile_data>& bins,
    std::vector<double> output_binfile_data::* sigMember,
    double ecgRate,
    double channelRate,
    double padSeconds = 0.3)
{
    size_t n = bins.size();
    PPGTemplatesResult out;
    out.templates.assign(n, {});
    out.iqrs.assign(n, {});
    out.kept.assign(n, {});
    out.peakCol.assign(n, -1);
    out.footCol.assign(n, -1);

    if (channelRate <= 0.0) return out;   // channel absent from this dataset

    int ppg_threads = std::min(8, static_cast<int>(n > 0 ? n : 1));
#pragma omp parallel for schedule(dynamic) num_threads(ppg_threads)
    for (int i = 0; i < static_cast<int>(n); ++i) {
        const auto& b = bins[i];
        if (b.bad_segment || (b.*sigMember).empty() || b.ch1.raw.size() < 2)
            continue;
        try {
            build_pulse_template_pair_windowed(
                b.*sigMember, channelRate,
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

/**
 * @brief  Foot-anchored pulse template for ONE arterial channel of ONE bin,
 *         per spec. Unlike PPG above, this is entirely self-contained --
 *         no borrowed ECG R-peaks. Pipeline:
 *           (1) self-detect systolic peaks directly on `signal` (local
 *               maxima, refractory window sized off channelRate);
 *           (2) slice [prevPeak, thisPeak] segments and batch them through
 *               find_foot_pulseox to locate each beat's true foot;
 *           (3) re-slice on a shared axis so every beat's OWN foot lands at
 *               a fixed column (padSamples), one foot-to-foot interval
 *               (median-length fallback for the last beat) plus trailing
 *               pad wide, NaN-padding short beats;
 *           (4) matched-filter QC (pulse_matched_filter): first 20 beats
 *               seed a provisional template, later beats rejected at >=5%
 *               normalized error;
 *           (5) column-wise NaN-skipping median -> template.
 */
static inline void build_arterial_template_foot_anchored(
    const std::vector<double>& signal,
    double channelRate,
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
    if (signal.empty() || channelRate <= 0.0) return;

    const int n = static_cast<int>(signal.size());

    // ---- (1) self-detect systolic peaks: local maxima, refractory window
    // generous enough for up to ~240 bpm (0.25 s minimum separation). -----
    const int minSep = std::max(1, static_cast<int>(std::llround(0.25 * channelRate)));
    std::vector<int> peaks;
    for (int i = 1; i + 1 < n; ++i) {
        if (std::isnan(signal[i - 1]) || std::isnan(signal[i]) || std::isnan(signal[i + 1])) continue;
        if (signal[i] > signal[i - 1] && signal[i] >= signal[i + 1]) {
            if (!peaks.empty() && (i - peaks.back()) < minSep) {
                if (signal[i] > signal[peaks.back()]) peaks.back() = i;   // stronger peak wins
            }
            else {
                peaks.push_back(i);
            }
        }
    }
    if (peaks.size() < 2) return;

    // ---- (2) [prevPeak, thisPeak] segments -> find_foot_pulseox --------
    const size_t nBeatsRaw = peaks.size() - 1;
    std::vector<std::vector<double>> segments(nBeatsRaw);
    for (size_t k = 0; k < nBeatsRaw; ++k) {
        const int a = peaks[k], b = peaks[k + 1];
        segments[k].assign(signal.begin() + a, signal.begin() + b + 1);
    }
    const FootResult feet = find_foot_pulseox(segments);

    std::vector<int> footAbs(nBeatsRaw);
    for (size_t k = 0; k < nBeatsRaw; ++k)
        footAbs[k] = peaks[k] + static_cast<int>(feet.idx[k]);

    // ---- (3) re-slice: every beat's OWN foot lands at column padSamples.
    // Width = one foot-to-foot interval (median-length fallback for the
    // last beat, which has no "next foot") plus lead/trail pad. -----------
    const int padSamples = std::max(0, static_cast<int>(std::llround(padSeconds * channelRate)));
    std::vector<int> gaps;
    gaps.reserve(nBeatsRaw);
    for (size_t k = 0; k + 1 < nBeatsRaw; ++k) gaps.push_back(footAbs[k + 1] - footAbs[k]);
    int medGap = std::max(1, padSamples * 4);
    if (!gaps.empty()) {
        std::vector<int> g = gaps;
        std::sort(g.begin(), g.end());
        medGap = std::max(1, g[g.size() / 2]);
    }
    const int width = padSamples + medGap + padSamples;

    std::vector<std::vector<double>> beats;
    beats.reserve(nBeatsRaw);
    for (size_t k = 0; k < nBeatsRaw; ++k) {
        const int foot = footAbs[k];
        const int start = foot - padSamples;
        std::vector<double> beat(width, NaN);
        for (int c = 0; c < width; ++c) {
            const int idx = start + c;
            if (idx >= 0 && idx < n) beat[c] = signal[idx];
        }
        beats.push_back(std::move(beat));
    }
    if (beats.empty()) return;

    // ---- (4) matched-filter QC: first 20 beats seed the template, the
    // rest tested by 5% normalized error (mirrors the PPG QC step). ------
    std::vector<std::vector<double>> filteredBeats;
    {
        const size_t nSeed = std::min<size_t>(20, beats.size());
        const std::vector<std::vector<double>> seedBeats(beats.begin(), beats.begin() + nSeed);
        const std::vector<double> provisional =
            pulse_matched_filter::buildTemplate(seedBeats, width);

        filteredBeats.assign(beats.begin(), beats.begin() + nSeed);
        for (size_t i = nSeed; i < beats.size(); ++i) {
            const double err = pulse_matched_filter::normalizedError(beats[i], provisional);
            if (err < 0.05) filteredBeats.push_back(beats[i]);
        }
        if (filteredBeats.empty()) filteredBeats = beats;   // don't discard all
    }

    // ---- (5) column-wise NaN-skipping median -> template. --------------
    outTemplate.assign(width, NaN);
    for (int c = 0; c < width; ++c) {
        std::vector<double> col;
        col.reserve(filteredBeats.size());
        for (const auto& sl : filteredBeats) {
            const double v = sl[c];
            if (!std::isnan(v)) col.push_back(v);
        }
        if (col.empty()) continue;
        const size_t nc = col.size();
        const size_t mid = nc / 2;
        std::nth_element(col.begin(), col.begin() + mid, col.end());
        const double hi = col[mid];
        if (nc % 2) {
            outTemplate[c] = hi;
        }
        else {
            // lower median = max of the left partition nth_element already
            // produced -- no second full sort needed.
            const double lo = *std::max_element(col.begin(), col.begin() + mid);
            outTemplate[c] = 0.5 * (lo + hi);
        }
    }

    // The foot is fixed by construction -- every beat's own foot was
    // shifted to this column when it was re-sliced in step (3).
    outFootCol = padSamples;
    {
        int pk = outFootCol; double pmax = -Inf;
        for (int i = outFootCol; i < width; ++i)
            if (!std::isnan(outTemplate[i]) && outTemplate[i] > pmax) { pmax = outTemplate[i]; pk = i; }
        outPeakCol = pk;
    }

    outIqr = (outFootCol >= 0)
        ? normalize_features::local_ratio_iqr(filteredBeats, outFootCol)
        : std::vector<double>(width, 0.0);

    outKeptBeats.reserve(filteredBeats.size());
    for (const auto& sl : filteredBeats) outKeptBeats.push_back(sl);
}

/**
 * @brief  Foot-anchored templates for every bin of ONE arterial channel
 *         (ABP / ART / ART_PULM), per spec. No ECG R-peaks involved --
 *         self-detected and self-anchored, see
 *         build_arterial_template_foot_anchored.
 *
 * @param bins        Input bins.
 * @param sigMember   Member-pointer selecting the channel (abpSignal /
 *                     artSignal / artPulmSignal).
 * @param channelRate Sample rate of that channel.
 * @param padSeconds  Lead-in before the foot / trail-out pad. 0.25 matches
 *                     the PPG/ECG convention.
 */
inline PPGTemplatesResult CreateArterialTemplates(
    const vector<output_binfile_data>& bins,
    std::vector<double> output_binfile_data::* sigMember,
    double channelRate,
    double padSeconds = 0.25)
{
    size_t n = bins.size();
    PPGTemplatesResult out;
    out.templates.assign(n, {});
    out.iqrs.assign(n, {});
    out.kept.assign(n, {});
    out.peakCol.assign(n, -1);
    out.footCol.assign(n, -1);

    if (channelRate <= 0.0) return out;   // channel absent from this dataset

    // find_foot_pulseox (called per-bin below) has its OWN internal
    // #pragma omp parallel for over beat rows. Without this cap, that
    // becomes a nested parallel region inside the per-bin loop just below --
    // either real thread oversubscription (if nesting is enabled somewhere
    // else in this process) or repeated fork/join overhead on every single
    // bin x channel call (if it isn't). omp_set_nested(0) forces the inner
    // region to always collapse to the calling thread; parallelism stays at
    // the per-bin level, where it's actually worth it.
    // (omp_set_max_active_levels is OpenMP 3.0+ and isn't available under
    // MSVC's default /openmp flag, which only implements OpenMP 2.0 --
    // omp_set_nested is the 2.0-era equivalent and does the same thing.)
    omp_set_nested(0);

    int threads = std::min(8, static_cast<int>(n > 0 ? n : 1));
#pragma omp parallel for schedule(dynamic) num_threads(threads)
    for (int i = 0; i < static_cast<int>(n); ++i) {
        const auto& b = bins[i];
        if (b.bad_segment || (b.*sigMember).empty()) continue;
        try {
            build_arterial_template_foot_anchored(
                b.*sigMember, channelRate, padSeconds,
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