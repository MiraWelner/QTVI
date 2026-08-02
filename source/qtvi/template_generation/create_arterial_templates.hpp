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
 *         filter.hpp): two-pass median-based rejection. Build a reference
 *         template as the median of ALL candidate beats, then reject any
 *         candidate whose normalized error against that reference exceeds
 *         5%. The final template is rebuilt from the survivors.
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
#include "ppg_matched_filter.hpp"   // derivativePulseLocations for arterial stage-1 census
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

    // ---- Matched-filter QC, two-pass, per spec:
    //   (a) build a REFERENCE template as the column-wise NaN-skipping
    //       median across ALL candidate beats. The median is robust to
    //       outliers without needing to pick a fixed "seed" count.
    //   (b) score every candidate against the reference by normalized
    //       error ||beat - ref|| / ||ref||; keep beats whose error is
    //       below 5%.
    // The final template below is then rebuilt from the survivors,
    // giving the two-pass: median-of-all -> reject high-error ->
    // re-median. Same wave-score pruning logic as ECG, adapted to PPG's
    // normalized-error metric. Falls back to keeping everything if the
    // filter would otherwise reject the whole set (degenerate reference).
    std::vector<std::vector<double>> filteredBeats;
    {
        const int w = static_cast<int>(aligned.beats.front().size());
        // (a) reference = column-wise median across ALL candidates.
        std::vector<double> reference(w, NaN);
        for (int c = 0; c < w; ++c) {
            std::vector<double> col;
            col.reserve(aligned.beats.size());
            for (const auto& sl : aligned.beats)
                if (c < (int)sl.size() && !std::isnan(sl[c])) col.push_back(sl[c]);
            if (col.empty()) continue;
            const size_t nc = col.size();
            const size_t mid = nc / 2;
            std::nth_element(col.begin(), col.begin() + mid, col.end());
            reference[c] = (nc % 2)
                ? col[mid]
                : 0.5 * (*std::max_element(col.begin(), col.begin() + mid) + col[mid]);
        }
        // (b) per-pulse accept/reject against the reference.
        filteredBeats.reserve(aligned.beats.size());
        for (const auto& bt : aligned.beats) {
            const double err = pulse_matched_filter::normalizedError(bt, reference);
            if (err < 0.05) filteredBeats.push_back(bt);
        }
        if (filteredBeats.empty()) filteredBeats = aligned.beats;   // degenerate ref -> keep all
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
 *           (1) self-detect systolic peaks in TWO steps:
 *               1a) derivative-max census (ppg_matched_filter): finds each
 *                   pulse's steepest upstroke, the sharpest and least-
 *                   variable landmark in an arterial waveform;
 *               1b) apex-walk: from each upstroke, walk forward to the
 *                   local maximum (the true systolic peak), bounded by
 *                   the next upstroke so we can't cross into the next beat.
 *           (2) slice [prevPeak, thisPeak] segments and batch them through
 *               find_foot_pulseox to locate each beat's true foot;
 *           (3) re-slice on a shared axis so every beat's OWN foot lands at
 *               a fixed column (padSamples), one foot-to-foot interval
 *               (median-length fallback for the last beat) plus trailing
 *               pad wide, NaN-padding short beats;
 *           (4) two-pass matched-filter QC:
 *               4a) reference template = column-wise median of ALL
 *                   candidate beats (robust to outliers, no need to pick
 *                   a fixed "seed" count);
 *               4b) score every candidate by normalized error against the
 *                   reference; keep beats whose error is below 5%.
 *           (5) column-wise NaN-skipping median across survivors ->
 *               final template. (This is the "re-median from survivors"
 *               half of the two-pass approach.)
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

    // ---- (1) self-detect systolic peaks in TWO steps, per spec:
    //   1a) derivative-max census: run ppg_matched_filter's derivative-max
    //       detector to get a rough list of where each pulse's steepest
    //       upstroke sits. The upstroke is the sharpest, least-variable
    //       part of an arterial pulse -- more reliable to detect than the
    //       apex or dicrotic notch, both of which vary beat-to-beat.
    //   1b) apex-walk: from each detected upstroke, walk forward through
    //       the signal to the local maximum -- the actual systolic peak.
    //       Bounded by the next upstroke's location so we can't overshoot
    //       into the following beat.
    // -------------------------------------------------------------------
    const int minSep = std::max(1, static_cast<int>(std::llround(0.25 * channelRate)));
    const std::vector<int> upstrokes =
        ppg_matched_filter::derivativePulseLocations(signal, minSep);
    if (upstrokes.size() < 2) return;

    std::vector<int> peaks;
    peaks.reserve(upstrokes.size());
    for (size_t k = 0; k < upstrokes.size(); ++k) {
        const int start = upstrokes[k];
        // Search up to the next upstroke (exclusive), or to the end of the
        // signal for the last one; cap at start + minSep as a safety belt
        // in case an upstroke got dropped and the "next" one is far away.
        const int hardEnd = (k + 1 < upstrokes.size())
            ? upstrokes[k + 1]
            : n;
        const int end = std::min(hardEnd, start + minSep);
        int pk = start;
        double pkVal = -Inf;
        for (int i = start; i < end && i < n; ++i) {
            const double v = signal[i];
            if (std::isnan(v)) continue;
            if (v > pkVal) { pkVal = v; pk = i; }
        }
        if (std::isfinite(pkVal)) peaks.push_back(pk);
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

    // ---- (4) Matched-filter QC, two-pass, per spec:
    //   4a) build a REFERENCE template as the column-wise NaN-skipping
    //       median across ALL candidate beats. The median is robust to
    //       outliers without needing to pick a fixed "seed" count.
    //   4b) score every candidate against the reference by normalized
    //       error ||beat - ref|| / ||ref||; keep beats whose error is
    //       below 5%.
    // The final template (step 5 below) is rebuilt from the survivors,
    // giving the two-pass: median-of-all -> reject high-error ->
    // re-median. Same wave-score pruning logic as ECG, adapted to PPG's
    // normalized-error metric. Falls back to keeping everything if the
    // filter would otherwise reject the whole set (degenerate reference).
    // ---------------------------------------------------------------------
    std::vector<std::vector<double>> filteredBeats;
    {
        // 4a) reference = column-wise median of ALL candidates.
        std::vector<double> reference(width, NaN);
        for (int c = 0; c < width; ++c) {
            std::vector<double> col;
            col.reserve(beats.size());
            for (const auto& sl : beats)
                if (!std::isnan(sl[c])) col.push_back(sl[c]);
            if (col.empty()) continue;
            const size_t nc = col.size();
            const size_t mid = nc / 2;
            std::nth_element(col.begin(), col.begin() + mid, col.end());
            reference[c] = (nc % 2)
                ? col[mid]
                : 0.5 * (*std::max_element(col.begin(), col.begin() + mid) + col[mid]);
        }

        // 4b) per-pulse accept/reject on normalized error against reference.
        filteredBeats.reserve(beats.size());
        for (const auto& bt : beats) {
            const double err = pulse_matched_filter::normalizedError(bt, reference);
            if (err < 0.05) filteredBeats.push_back(bt);
        }
        if (filteredBeats.empty()) filteredBeats = beats;   // degenerate ref -> keep all
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