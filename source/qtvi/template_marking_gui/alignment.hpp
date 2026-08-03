#pragma once
//
// alignment.hpp
//
// Per-bin beat alignment for ECG and PPG. Runs BEFORE any normalization
// or template averaging.
//
//   1) Slice one beat per RR window as 0.4*RR before R to 1.3*RR after R (integer-truncated).
//   2) Tukey outlier rejection: ECG is based on RR length and distance from max to min, PPG is based on 50% upslope location and distance from max to min.
//   3) Horizontal align the ECG such that R peaks are aligned, if it is the Q peak screening then after align by Q peak
//   4) Vertical DC shift:
//        - ECG: match each beat's PR-baseline mean to the reference beat's.
//        - PPG: match each beat's foot-baseline mean (a small window around
//               its own foot column) to the reference beat's.
//

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>
#include <functional>
#include "feature_marks.hpp"


namespace alignment {

    // Which isoelectric segment supplied a beat's DC baseline. TP (T-end to
    // next P-onset) is preferred over PQ; NONE means neither was usable for
    // this beat -- it is left un-shifted, and CreateEcgTemplates.hpp excludes
    // it from the per-sample amplitude aggregation (median/std) rather than
    // silently contribute an unreliable, unadjusted amplitude.
    enum class BaselineSource { TP, PQ, NONE };
    constexpr double percent_interval_preceeding_rpeak = 0.3; //how far before the R peak the snip goes, in terms of percent of the RR interval length
    constexpr double percent_interval_following_rpeak = 1.5;   //how far after the R peak the snip goes, in terms of percent of the RR interval length

    enum class AnchorType { P_ONSET, P_PEAK, Q_ONSET, R_PEAK, J_POINT, T_PEAK };

    // Sample counts for a given RR (integer-truncated).
    inline int64_t rr_before_samples(int64_t rr) {
        return static_cast<int64_t>(percent_interval_preceeding_rpeak * rr);
    }
    inline int64_t rr_after_samples(int64_t rr) {
        return static_cast<int64_t>(percent_interval_following_rpeak * rr);
    }

    // Tukey outlier mask. keep[i] = false iff values[i] falls outside
    // [Q1 - k*IQR, Q3 + k*IQR]. NaN entries are left keep=true; callers
    // must skip them explicitly when needed.
    inline std::vector<bool> keep_within_tukey(
        const std::vector<double>& values, double k)
    {
        std::vector<bool> keep(values.size(), true);
        if (values.size() < 4) return keep;

        std::vector<double> sorted;
        sorted.reserve(values.size());
        for (double v : values) if (!std::isnan(v)) sorted.push_back(v);
        if (sorted.size() < 4) return keep;
        std::sort(sorted.begin(), sorted.end());

        auto quantile = [&](double q) {
            const double h = q * (sorted.size() - 1);
            const size_t lo = static_cast<size_t>(std::floor(h));
            const size_t hi = static_cast<size_t>(std::ceil(h));
            const double frac = h - lo;
            return sorted[lo] * (1.0 - frac) + sorted[hi] * frac;
            };
        const double q1 = quantile(0.25);
        const double q3 = quantile(0.75);
        const double iqr = q3 - q1;
        if (iqr <= 0.0) return keep;

        const double lo_b = q1 - k * iqr;
        const double hi_b = q3 + k * iqr;
        for (size_t i = 0; i < values.size(); ++i) {
            if (std::isnan(values[i])) continue;
            if (values[i] < lo_b || values[i] > hi_b) keep[i] = false;
        }
        return keep;
    }
    struct ecg_beat_set {
        //a set of individually-aligned beats that will become a template.
        std::vector<std::vector<double>> beats;
        std::vector<size_t> r_indices;
        std::vector<int>    rr_lens;
        int    median_length = -1;
        size_t total_beats = 0;
        int r_aligned_col = -1;
        int q_aligned_col = -1;
        int ref_beat_index = -1;
        // Parallel to `beats`: which segment supplied each beat's DC
        // baseline (TP/PQ/NONE). Kept in sync with `beats` through every
        // apply_mask() compaction, including the wave-score pruning pass.
        std::vector<BaselineSource> baseline_source;
        // Parallel to `beats`: PQ_level - TP_level whenever BOTH stage
        // estimates succeeded for that beat (NaN otherwise). QC metric per
        // spec I-1; not used to drive any decision, just recorded.
        std::vector<double> tp_pq_delta;
    };

    inline ecg_beat_set extract_beats_and_align(const std::vector<double>& signal, const std::vector<size_t>& rPeaks, double fs) {
        ecg_beat_set out;
        const int64_t N = static_cast<int64_t>(signal.size());
        if (N == 0 || rPeaks.size() < 2) {
            fprintf(stderr, "[align] no beats: signal=%lld rPeaks=%zu\n",
                static_cast<long long>(N), rPeaks.size());
            return out;
        }

        // Compact the parallel (beats, r_indices, rr_lens) vectors, keeping
        // only entries where keep[i] is true.
        auto apply_mask = [&](const std::vector<bool>& keep) {
            std::vector<std::vector<double>> kb;
            std::vector<size_t> kr;
            std::vector<int>    km;
            std::vector<BaselineSource> ks;
            std::vector<double> kd;
            const bool haveSrc = out.baseline_source.size() == out.beats.size();
            const bool haveDelta = out.tp_pq_delta.size() == out.beats.size();
            for (size_t i = 0; i < keep.size(); ++i) {
                if (!keep[i]) continue;
                kb.push_back(std::move(out.beats[i]));
                kr.push_back(out.r_indices[i]);
                km.push_back(out.rr_lens[i]);
                if (haveSrc) ks.push_back(out.baseline_source[i]);
                if (haveDelta) kd.push_back(out.tp_pq_delta[i]);
            }
            out.beats = std::move(kb);
            out.r_indices = std::move(kr);
            out.rr_lens = std::move(km);
            if (haveSrc) out.baseline_source = std::move(ks);
            if (haveDelta) out.tp_pq_delta = std::move(kd);
            };

        // ---- slice every beat ------------------------------------------
        for (size_t i = 0; i + 1 < rPeaks.size(); ++i) {
            const int64_t r0 = static_cast<int64_t>(rPeaks[i]);
            const int64_t rr = static_cast<int64_t>(rPeaks[i + 1]) - r0;
            if (rr <= 3) continue;

            const int64_t before = rr_before_samples(rr);
            const int64_t after = rr_after_samples(rr);
            const int64_t len = before + after;
            const int64_t start = r0 - before;
            const int64_t end = r0 + after;

            std::vector<double> beat(static_cast<size_t>(len),
                std::numeric_limits<double>::quiet_NaN());
            const int64_t cs = std::max<int64_t>(0, start);
            const int64_t ce = std::min<int64_t>(N, end);
            for (int64_t k = cs; k < ce; ++k)
                beat[static_cast<size_t>(k - start)] = signal[static_cast<size_t>(k)];

            out.beats.push_back(std::move(beat));
            out.r_indices.push_back(rPeaks[i]);
            out.rr_lens.push_back(static_cast<int>(rr));
        }
        if (out.beats.empty()) return out;
        out.baseline_source.assign(out.beats.size(), BaselineSource::NONE);

        //Tukey rejection: R-R interval (1.5*IQR)
        {
            std::vector<double> lens_d(out.rr_lens.begin(), out.rr_lens.end());
            apply_mask(keep_within_tukey(lens_d, 1.5));
        }
        if (out.beats.empty()) return out;
        //Tukey rejection: R peak from template min distance (1.5*IQR)
        {
            std::vector<double> amps;
            amps.reserve(out.beats.size());
            for (size_t i = 0; i < out.beats.size(); ++i) {
                const int r_col = static_cast<int>(rr_before_samples(out.rr_lens[i]));
                const auto& beat = out.beats[i];

                if (r_col >= (int)beat.size() || std::isnan(beat[r_col])) {
                    amps.push_back(std::numeric_limits<double>::quiet_NaN());
                    continue;
                }
                double min_val = std::numeric_limits<double>::infinity();
                for (double v : beat)
                    if (!std::isnan(v) && v < min_val) min_val = v;

                amps.push_back(std::isfinite(min_val)
                    ? beat[r_col] - min_val
                    : std::numeric_limits<double>::quiet_NaN());
            }
            auto keepA = keep_within_tukey(amps, 1.5);
            for (size_t i = 0; i < keepA.size(); ++i)
                if (std::isnan(amps[i])) keepA[i] = false;
            apply_mask(keepA);
        }
        out.total_beats = out.beats.size();
        if (out.beats.empty()) return out;

        // the Sangala document says use mode, in this case we use median length
        {
            std::vector<int> lens = out.rr_lens;
            std::sort(lens.begin(), lens.end());
            out.median_length = lens[lens.size() / 2];
        }

        // ---- reference beat = first median-length beat -----------------
        // Used only as the PR-baseline DC-alignment reference. No RMS shape
        // clustering: the final template is the column-wise median.
        for (size_t i = 0; i < out.beats.size(); ++i) {
            if (out.rr_lens[i] == out.median_length) {
                out.ref_beat_index = static_cast<int>(i);
                break;
            }
        }
        // ---- Pass 1: R-align on shared axis ----------------------------
        int max_rr_len = 0;
        for (int L : out.rr_lens)
            if (L > max_rr_len) max_rr_len = L;
        if (max_rr_len <= 0) return out;

        const int R_anchor = static_cast<int>(rr_before_samples(max_rr_len));
        const int shared_w = R_anchor + static_cast<int>(rr_after_samples(max_rr_len));

        const double NaND = std::numeric_limits<double>::quiet_NaN();
        std::vector<std::vector<double>> aligned;
        aligned.reserve(out.beats.size());
        for (size_t i = 0; i < out.beats.size(); ++i) {
            const auto& b = out.beats[i];
            const int L = static_cast<int>(b.size());
            const int r_in_beat = static_cast<int>(rr_before_samples(out.rr_lens[i]));
            const int prepend = R_anchor - r_in_beat;
            std::vector<double> a(shared_w, NaND);
            for (int k = 0; k < L; ++k) {
                const int dst = prepend + k;
                if (dst >= 0 && dst < shared_w) a[dst] = b[k];
            }
            aligned.push_back(std::move(a));
        }
        out.beats = std::move(aligned);
        out.r_aligned_col = R_anchor;

        // ---- Pass 3: two-stage TP-then-PQ vertical DC alignment --------
        // Stage 1 (TP): this beat's own T-end -> the NEXT beat's P-onset.
        // Stage 2 (PQ): this beat's own P-end -> this beat's own Q-onset.
        // PQ is the higher-priority reference (spec): when BOTH estimates
        // are available, the beat is shifted using PQ's level, not TP's --
        // TP is used only when PQ's landmarks aren't usable. Both are always
        // attempted whenever possible (not short-circuited on TP success)
        // so the TP-to-PQ delta can be recorded as a per-beat QC metric even
        // on beats where PQ ends up winning.
        //
        // Both estimators report level as the MEDIAN over the flattest
        // (lowest-variance) min_w-wide sub-window of their respective
        // landmark-bounded segment -- flatness search still uses variance
        // (cheapest way to find "flat"), only the reported level is median
        // rather than mean.
        //
        // The next beat's own R, in this beat's own SHARED aligned column
        // space, sits at R_anchor + rr_lens[i] (this beat's own R was moved
        // to R_anchor by Pass 1; the next R was the same distance further in
        // the original slice, so the same shift lands it there) -- pure
        // arithmetic, no extra detection needed to locate it.
        //
        // Degrades gracefully: baseline_source[i] records which segment (if
        // any) supplied this beat's baseline; tp_pq_delta[i] records
        // PQ_level - TP_level whenever both succeeded (NaN otherwise). NONE
        // means neither was usable; that beat is left un-shifted, and
        // CreateEcgTemplates.hpp excludes it from the per-sample amplitude
        // aggregation (median/std) rather than silently contribute an
        // unreliable, unadjusted amplitude.
        if (out.ref_beat_index >= 0 && out.ref_beat_index < static_cast<int>(out.beats.size()))
        {
            const int min_w = std::max(1, static_cast<int>(std::lround(0.010 * fs)));  // 10 ms window

            // Median over the lowest-variance min_w-wide sub-window of
            // [lo, hi] (inclusive, in this beat's own aligned column space).
            auto flattest_median = [&](const std::vector<double>& beat, int lo, int hi)
                -> std::pair<double, bool>
                {
                    const int hi_bound = std::min(hi, static_cast<int>(beat.size()) - min_w);
                    if (hi_bound < lo)
                        return { std::numeric_limits<double>::quiet_NaN(), false };

                    double best_var = std::numeric_limits<double>::infinity();
                    int best_start = -1;
                    for (int s = lo; s <= hi_bound; ++s) {
                        double sum = 0.0, sumsq = 0.0; int n = 0;
                        for (int k = s; k < s + min_w; ++k) {
                            const double v = beat[k];
                            if (std::isnan(v)) continue;
                            sum += v; sumsq += v * v; ++n;
                        }
                        if (n < static_cast<int>(min_w * 0.7)) continue;
                        const double mean = sum / n;
                        const double var = sumsq / n - mean * mean;
                        if (var < best_var) { best_var = var; best_start = s; }
                    }
                    if (best_start < 0)
                        return { std::numeric_limits<double>::quiet_NaN(), false };

                    std::vector<double> win;
                    win.reserve(min_w);
                    for (int k = best_start; k < best_start + min_w; ++k)
                        if (!std::isnan(beat[k])) win.push_back(beat[k]);
                    if (win.empty())
                        return { std::numeric_limits<double>::quiet_NaN(), false };
                    std::sort(win.begin(), win.end());
                    const size_t m = win.size() / 2;
                    const double med = (win.size() % 2 == 0)
                        ? 0.5 * (win[m - 1] + win[m]) : win[m];
                    return { med, true };
                };

            // Stage 1, TP window: [this beat's own T-end, next beat's
            // P-onset]. {-1,-1} if the next R isn't inside this beat's slice
            // or the landmarks come back invalid/out of order (e.g. HR fast
            // enough that T runs into the next P).
            //
            // detect_p_peak/detect_q_begin (and likely other single-beat
            // detectors) search from absolute column 0 up to r_idx -- safe
            // only when r_idx is the array's OWN first/only R. Calling them
            // directly on this beat's SHARED array with next_r reaches back
            // across THIS beat's own, much-larger QRS spike, which wins any
            // min/max search meant to find the next beat's P or Q. Rather
            // than modify those detectors (used correctly elsewhere in the
            // true single-beat context), extract a small, genuinely-
            // isolated local window around next_r first, so every detector
            // sees exactly the kind of array it was designed for.
            auto tp_window = [&](size_t i) -> std::pair<int, int> {
                const auto& beat = out.beats[i];
                const int N = static_cast<int>(beat.size());
                const int t_end = FeatureMarks::detect_t_end(beat, R_anchor, fs);
                const int next_r = R_anchor + out.rr_lens[i];
                if (next_r <= 0 || next_r >= N) return { -1, -1 };

                // Local window: comfortably covers a PR interval before
                // next_r (400 samples ~ 400ms at 1000 Hz, generous margin)
                // with a little slack past it, clamped to this beat's own
                // array bounds.
                const int localLo = std::max(0, next_r - 400);
                const int localHi = std::min(N, next_r + 50);
                if (localHi - localLo < 10) return { -1, -1 };
                std::vector<double> local(beat.begin() + localLo, beat.begin() + localHi);
                const int local_r = next_r - localLo;

                const int p_pk_local = FeatureMarks::detect_p_peak(local, local_r);
                const int p_on_local = FeatureMarks::compute_p_begin(local, p_pk_local, fs, local_r);
                if (p_on_local < 0) return { -1, -1 };
                const int p_on = p_on_local + localLo;   // back to shared column space

                if (t_end < 0 || p_on < 0 || t_end >= N || p_on >= N || p_on <= t_end)
                    return { -1, -1 };
                return { t_end, p_on };
                };

            // Stage 2, PQ window: [this beat's own P-end, this beat's own
            // Q-onset]. Both landmarks are real, fitted anchors (Phase A/B),
            // not a fixed lookback -- matches the spec's literal PQ
            // definition ("end of P to immediately before Q onset").
            // {-1,-1} if the landmarks come back invalid/out of order.
            auto pq_window = [&](size_t i) -> std::pair<int, int> {
                const auto& beat = out.beats[i];
                const int N = static_cast<int>(beat.size());
                const int p_end_seed = FeatureMarks::detect_p_end(beat, R_anchor);
                const int p_end = FeatureMarks::compute_p_end(beat, p_end_seed, fs, R_anchor);
                const int q_seed = FeatureMarks::detect_q_begin(beat, R_anchor);
                const int q_on = FeatureMarks::compute_q_onset(beat, q_seed, fs, R_anchor);
                if (p_end < 0 || q_on < 0 || p_end >= N || q_on >= N || q_on <= p_end)
                    return { -1, -1 };
                return { p_end, q_on };
                };

            // Compute both stage estimates for beat i (never short-circuited,
            // so the delta can be recorded even when PQ will be used).
            auto both_levels = [&](size_t i) -> std::tuple<double, bool, double, bool> {
                double tpLvl = std::numeric_limits<double>::quiet_NaN(), pqLvl = tpLvl;
                bool tpOk = false, pqOk = false;
                const auto [tlo, thi] = tp_window(i);
                if (tlo >= 0) { const auto [lvl, ok] = flattest_median(out.beats[i], tlo, thi); tpLvl = lvl; tpOk = ok; }
                const auto [plo, phi] = pq_window(i);
                if (plo >= 0) { const auto [lvl, ok] = flattest_median(out.beats[i], plo, phi); pqLvl = lvl; pqOk = ok; }
                return { tpLvl, tpOk, pqLvl, pqOk };
                };

            // PQ is the higher-priority reference: use it whenever available,
            // falling back to TP only when PQ's landmarks aren't usable.
            auto pick = [&](double tpLvl, bool tpOk, double pqLvl, bool pqOk,
                BaselineSource& src) -> double {
                    if (pqOk) { src = BaselineSource::PQ; return pqLvl; }
                    if (tpOk) { src = BaselineSource::TP; return tpLvl; }
                    src = BaselineSource::NONE;
                    return std::numeric_limits<double>::quiet_NaN();
                };

            const auto [refTp, refTpOk, refPq, refPqOk] = both_levels(static_cast<size_t>(out.ref_beat_index));
            BaselineSource refSrc = BaselineSource::NONE;
            const double target = pick(refTp, refTpOk, refPq, refPqOk, refSrc);

            // If the default reference beat has no usable PQ/TP, don't abandon
            // the whole bin: find any beat with a usable PQ (else TP) and use
            // it as the reference so beats that DO resolve still align.
            double refTarget = target;
            if (refSrc == BaselineSource::NONE) {
                for (size_t i = 0; i < out.beats.size(); ++i) {
                    const auto [tTp, tTpOk, tPq, tPqOk] = both_levels(i);
                    if (tPqOk) { refTarget = tPq; refSrc = BaselineSource::PQ; break; }
                    if (tTpOk && refSrc == BaselineSource::NONE) {
                        refTarget = tTp; refSrc = BaselineSource::TP;
                    }
                }
            }

            out.tp_pq_delta.assign(out.beats.size(), std::numeric_limits<double>::quiet_NaN());

            // Diagnostic accumulators (spec-neutral: recorded only, doesn't
            // change any beat).
            int diag_pq = 0, diag_tp = 0, diag_none = 0;
            std::vector<double> diag_shifts;
            diag_shifts.reserve(out.beats.size());

            if (refSrc != BaselineSource::NONE && !std::isnan(refTarget)) {
                for (size_t i = 0; i < out.beats.size(); ++i) {
                    const auto [tpLvl, tpOk, pqLvl, pqOk] = both_levels(i);
                    if (tpOk && pqOk) out.tp_pq_delta[i] = pqLvl - tpLvl;   // QC metric, whenever both succeed

                    BaselineSource src = BaselineSource::NONE;
                    const double lvl = pick(tpLvl, tpOk, pqLvl, pqOk, src);
                    out.baseline_source[i] = src;
                    if (src == BaselineSource::PQ) ++diag_pq;
                    else if (src == BaselineSource::TP) ++diag_tp;
                    else ++diag_none;
                    if (src == BaselineSource::NONE || std::isnan(lvl))
                        continue;   // graceful degrade: leave un-shifted, flagged NONE
                    const double d = lvl - refTarget;
                    diag_shifts.push_back(d);
                    if (d == 0.0) continue;
                    for (double& v : out.beats[i])
                        if (!std::isnan(v)) v -= d;
                }
            }
            // One line per aligned matrix. Shows whether shifts are sane
            // (small std, tight range) or wild (large std, big range).
            {
                const char* refTag = (refSrc == BaselineSource::PQ) ? "PQ"
                    : (refSrc == BaselineSource::TP) ? "TP" : "NONE";
                double smin = std::numeric_limits<double>::infinity();
                double smax = -std::numeric_limits<double>::infinity();
                double ssum = 0.0, ssumsq = 0.0;
                for (double s : diag_shifts) {
                    if (s < smin) smin = s;
                    if (s > smax) smax = s;
                    ssum += s; ssumsq += s * s;
                }
                const size_t nS = diag_shifts.size();
                const double smean = (nS > 0) ? ssum / nS : 0.0;
                const double sstd = (nS > 1)
                    ? std::sqrt(std::max(0.0, ssumsq / nS - smean * smean)) : 0.0;
                double smed = 0.0;
                if (!diag_shifts.empty()) {
                    auto tmp = diag_shifts;
                    std::sort(tmp.begin(), tmp.end());
                    smed = tmp[tmp.size() / 2];
                }
            }
        }
        // ---- Wave-score pruning ----------------------------------------
        // Iterative template-matching QC, run AFTER R-alignment and PQ-
        // baseline DC-shift (so beats are already on a common axis and
        // vertical offset). Three tightening passes (4.0 -> 3.0 -> 2.5 SD);
        // each pass rebuilds a column-wise NaN-skipping median template from
        // the CURRENT survivors, then drops any beat whose:
        //   (a) Pearson correlation with that template (over columns where
        //       both are non-NaN) is below 0.30, OR
        //   (b) worst single-sample deviation from the template exceeds
        //       (pass threshold) x the template's own per-sample residual SD
        //       (one scalar per pass, not per-column, to avoid overfitting a
        //       591-column threshold vector to a handful of beats).
        // ASSUMPTION: the spec's three SD numbers are read as three
        // successive tightening passes, not three simultaneous per-metric
        // thresholds -- flag if that's not the intended reading.
        {
            auto column_median = [&](const std::vector<std::vector<double>>& beats, int width) {
                std::vector<double> tmpl(width, NaND);
                std::vector<double> col;
                col.reserve(beats.size());
                for (int c = 0; c < width; ++c) {
                    col.clear();
                    for (const auto& b : beats)
                        if (!std::isnan(b[c])) col.push_back(b[c]);
                    if (col.empty()) continue;
                    std::sort(col.begin(), col.end());
                    const size_t nc = col.size();
                    tmpl[c] = (nc % 2 == 0)
                        ? 0.5 * (col[nc / 2 - 1] + col[nc / 2])
                        : col[nc / 2];
                }
                return tmpl;
                };

            auto pearson = [&](const std::vector<double>& a, const std::vector<double>& b) -> double {
                double sa = 0, sb = 0, saa = 0, sbb = 0, sab = 0; int n = 0;
                for (size_t k = 0; k < a.size(); ++k) {
                    if (std::isnan(a[k]) || std::isnan(b[k])) continue;
                    sa += a[k]; sb += b[k]; saa += a[k] * a[k];
                    sbb += b[k] * b[k]; sab += a[k] * b[k]; ++n;
                }
                if (n < 4) return 0.0;   // too few overlapping samples to trust
                const double ma = sa / n, mb = sb / n;
                const double cov = sab / n - ma * mb;
                const double va = saa / n - ma * ma, vb = sbb / n - mb * mb;
                if (va <= 0.0 || vb <= 0.0) return 0.0;
                return cov / std::sqrt(va * vb);
                };

            const double corr_min = 0.30;
            const double sd_thresholds[3] = { 4.0, 3.0, 2.5 };

            // Freeze the residual-SD SCALE once, from the population as it
            // stands right before wave-score pruning starts (post length/
            // amplitude Tukey). Recomputing this scale fresh each pass from
            // an ever-shrinking, ever-more-homogeneous survivor set makes it
            // collapse toward zero -- the absolute threshold (sdThresh x SD)
            // then shrinks faster than intended and the last pass rejects
            // everything, including clean beats. Only the TEMPLATE (shape)
            // and the surviving set are allowed to refine across passes;
            // the yardstick they're measured against does not shrink with
            // them.
            double frozenSD = -1.0;
            {
                const std::vector<double> tmpl0 = column_median(out.beats, shared_w);
                double sumsq = 0.0; long n = 0;
                for (const auto& b : out.beats)
                    for (int c = 0; c < shared_w; ++c) {
                        const double v = b[c], t = tmpl0[c];
                        if (std::isnan(v) || std::isnan(t)) continue;
                        sumsq += (v - t) * (v - t); ++n;
                    }
                if (n >= 2) frozenSD = std::sqrt(sumsq / (n - 1));
            }

            if (frozenSD > 0.0) {
                for (double sdThresh : sd_thresholds) {
                    if (out.beats.size() < 4) break;   // not enough left to keep pruning meaningfully

                    const std::vector<double> tmpl = column_median(out.beats, shared_w);

                    std::vector<bool> keep(out.beats.size(), true);
                    for (size_t i = 0; i < out.beats.size(); ++i) {
                        const double r = pearson(out.beats[i], tmpl);
                        // RMS deviation over this beat's own overlapping columns --
                        // NOT the single worst sample. With hundreds of columns the
                        // max single-sample deviation is always several SDs out
                        // (extreme-value statistics), regardless of whether the
                        // beat is a good match; RMS is the well-behaved, standard
                        // "how many SDs away is this beat" measure.
                        double sumsq = 0.0; int n = 0;
                        for (int c = 0; c < shared_w; ++c) {
                            const double v = out.beats[i][c], t = tmpl[c];
                            if (std::isnan(v) || std::isnan(t)) continue;
                            sumsq += (v - t) * (v - t); ++n;
                        }
                        const double rms = (n > 0) ? std::sqrt(sumsq / n) : 0.0;
                        keep[i] = (r >= corr_min) && (rms <= sdThresh * frozenSD);
                    }
                    apply_mask(keep);
                }
            }
        }

        out.q_aligned_col = out.r_aligned_col;
        return out;
    }

    inline int find_q_column(const std::vector<double>& b, int R_anchor, double fs) {
        /*Finds Q column by negative 1st derivative followed by positive first dirivative and the whole thing
        has positive second derivative. Closets to R peak wins. If none found, */
        const int Wsh = static_cast<int>(b.size());
        if (Wsh == 0 || R_anchor < 1 || fs <= 0.0) return -1;

        const int q_window = std::max(2, static_cast<int>(std::lround(0.100 * fs)));
        const int lo = std::max(1, R_anchor - q_window);

        for (int k = R_anchor - 1; k > lo; --k) {
            if (k + 1 >= Wsh) continue;
            if (std::isnan(b[k - 1]) || std::isnan(b[k]) || std::isnan(b[k + 1])) continue;

            const double d1_prev = b[k] - b[k - 1];              // slope into k
            const double d1_next = b[k + 1] - b[k];              // slope out of k
            const double d2 = b[k + 1] - 2.0 * b[k] + b[k - 1];  // curvature at k

            if (d1_prev <= 0.0 && d1_next >= 0.0 && d2 > 0.0) return k;
        }
        std::cout << "No Q peak found in QAlign, falling back to steepest slope in R-upstroke region\n";
        int fallback = -1;
        double maxSlope = -std::numeric_limits<double>::infinity();
        for (int k = lo; k < R_anchor && k + 1 < Wsh; ++k) {
            if (std::isnan(b[k]) || std::isnan(b[k + 1])) continue;
            const double s = b[k + 1] - b[k];
            if (s > maxSlope) { maxSlope = s; fallback = k; }
        }
        return fallback;
    }

    // =========================================================================
    // Q-align a bin's cached (R-aligned) snippets and re-median.
    //
    // Mirrors R-alignment's structure: R-align picks a reference "median
    // snippet" (the first beat whose finite length equals the median) and
    // aligns to it. Here we do the same, but on Q. The alignment reference is
    // the location of the R peak on the beat of median length.  
    // find_q only looks in the R-upstroke region. We find Q on that median, then
    // shift every snippet so its own Q lands on that marker, and re-median the
    // shifted snippets into the Q template. The snippets stay full length
    // (never cropped); the window only bounds where the marker is found.
    //
    // Because the snippets are re-aligned on Q, each beat's R moves by its own
    // shift, so the Q template's R spike sits at R_anchor + median(shift). That
    // column is returned in r_col -- computed from the applied shifts, not a
    // window search (which would risk landing on Q or S).
    //
    // R_anchor is the R-pass R column; fs is the ECG rate; refMedian is the
    // R-aligned template used to locate the Q marker. compute_iqr is retained
    // for call-site compatibility but no longer gates anything -- the gray-band
    // IQR (a robust spread, not a standard deviation) is always filled.
    // q_aligned_col is the marker Q landed on; r_col is the R fiducial (R_anchor).
    struct QAlignResult {
        std::vector<double> tmpl;
        std::vector<double> iqr;
        int q_aligned_col = -1;
        int r_col = -1;
    };

    inline QAlignResult align_beat_matrix(
        const std::vector<std::vector<double>>& beatsIn,
        int R_anchor, double fs, bool compute_iqr,
        const std::vector<double>& ref_beat_of_median_length,
        const std::function<int(const std::vector<double>&)>& locate)
    {
        QAlignResult res;
        if (beatsIn.empty()) return res;

        std::vector<std::vector<double>> beats = beatsIn;
        const int Wsh = static_cast<int>(beats.front().size());

        // Q marker = Q of the column-wise median of all beats (== the R
        // template). Every snippet is then shifted so its own Q lands here.
        const int marker = locate(ref_beat_of_median_length);
        res.q_aligned_col = (marker >= 0) ? marker : R_anchor;

        std::vector<int> shifts;   // per-beat Q-align shift (for R's new column)
        if (marker >= 0) {
            const double NaNv = std::numeric_limits<double>::quiet_NaN();
            shifts.reserve(beats.size());
            for (size_t i = 0; i < beats.size(); ++i) {
                const int mi = locate(beats[i]);
                if (mi < 0) continue;
                const int shift = marker - mi;   // may be negative
                shifts.push_back(shift);
                if (shift == 0) continue;
                std::vector<double> a(Wsh, NaNv);
                for (int k = 0; k < Wsh; ++k) {
                    const int dst = k + shift;
                    if (dst >= 0 && dst < Wsh) a[dst] = beats[i][k];
                }
                beats[i] = std::move(a);
            }
        }

        // Column-wise median (+ optional IQR spread) over the Q-aligned beats.
        // One gather per column, reused scratch, nth_element for the median.
        res.tmpl.assign(Wsh, std::numeric_limits<double>::quiet_NaN());
        res.iqr.assign(Wsh, 0.0);

        std::vector<double> col;
        col.reserve(beats.size());
        for (int c = 0; c < Wsh; ++c) {
            col.clear();
            for (const auto& b : beats) {
                const double v = b[c];
                if (!std::isnan(v)) col.push_back(v);
            }
            const size_t nc = col.size();
            if (nc == 0) continue;

            const size_t mid = nc / 2;
            std::nth_element(col.begin(), col.begin() + mid, col.end());
            const double hi = col[mid];
            res.tmpl[c] = (nc % 2 == 0)
                ? 0.5 * (*std::max_element(col.begin(), col.begin() + mid) + hi)
                : hi;

            if (nc >= 2) {
                // Robust symmetric spread: the IQR (25th-75th percentile)
                // . Not inflated by the
                // few outlier / 1-sample-off beats at the R spike, so the band
                // tracks the median line instead of ballooning at R.
                std::vector<double> s(col);
                const size_t q1i = nc / 4;
                const size_t q3i = (3 * nc) / 4;
                std::nth_element(s.begin(), s.begin() + q1i, s.end());
                const double q1 = s[q1i];
                std::nth_element(s.begin() + q1i, s.begin() + q3i, s.end());
                const double q3 = s[q3i];
                res.iqr[c] = q3 - q1;
            }
        }

        // Passed-in R, tracked through the Q-align shift (median of the applied shifts).
        int med_shift = 0;
        if (!shifts.empty()) {
            std::sort(shifts.begin(), shifts.end());
            med_shift = shifts[shifts.size() / 2];
        }
        int rc = std::clamp(R_anchor + med_shift, 0, std::max(0, Wsh - 1));

        // Snap to the local peak within +/-5 ms of that passed-in position.
        // The window is far too tight to reach Q or S, so it only cleans up
        // sub-window drift of the R spike.
        const int w5 = std::max(1, static_cast<int>(std::lround(0.005 * fs)));
        const int rlo = std::max(0, rc - w5);
        const int rhi = std::min(Wsh - 1, rc + w5);
        double rbest = -std::numeric_limits<double>::infinity();
        for (int i = rlo; i <= rhi; ++i) {
            if (!std::isnan(res.tmpl[i]) && res.tmpl[i] > rbest) { rbest = res.tmpl[i]; rc = i; }
        }
        res.r_col = rc;
        return res;
    }

    // =========================================================================
    struct PpgBeatSet {
        std::vector<std::vector<double>> beats;   // NaN-padded, 50%-upslope-aligned
        std::vector<int> peak_cols;               // per-beat systolic peak column (varies)
        std::vector<int> foot_cols;               // per-beat foot column (varies)
        int    median_length = -1;
        size_t total_beats = 0;
        int    up50_aligned_col = -1;             // shared column all half-height points land on
        int    foot_aligned_col = -1;             // feet scatter (per-beat); not a shared column
        int    peak_aligned_col = -1;             // peaks scatter (per-beat); not a shared column
        int    ref_beat_index = -1;
    };

    inline PpgBeatSet extract_ppg_beats_and_align(const std::vector<double>& signal, const std::vector<size_t>& rPeaks, double fs)
    {
        PpgBeatSet out;
        const int64_t N = static_cast<int64_t>(signal.size());
        if (N == 0 || rPeaks.size() < 2) return out;

        struct Raw { std::vector<double> data; int peak; int foot; int up50; };
        std::vector<Raw> raw;
        std::vector<int> rr_lens;
        raw.reserve(rPeaks.size());

        // Compact the parallel (raw, rr_lens) vectors, keeping only entries
        // where keep[i] is true.
        auto apply_mask = [&](const std::vector<bool>& keep) {
            std::vector<Raw> filt;
            std::vector<int> filt_lens;
            filt.reserve(raw.size());
            filt_lens.reserve(raw.size());
            for (size_t i = 0; i < raw.size(); ++i)
                if (keep[i]) {
                    filt.push_back(std::move(raw[i]));
                    filt_lens.push_back(rr_lens[i]);
                }
            raw.swap(filt);
            rr_lens.swap(filt_lens);
            };

        // Physiological cap on the RR used to SIZE the beat window: 3.0 s
        // (20 bpm). On long signal dropouts/artifacts the peak-finder can
        // report an enormous RR; a single such beat sizes the shared window
        // (Pass 1: shared_w = up50_anchor + max_tail) for the WHOLE bin,
        // NaN-padding every beat out to tens of thousands of columns and
        // hanging the template build. Clamp the window-sizing RR (rr_w).
        // The beat is still sliced and kept -- only its width is clamped --
        // and the raw rr is still recorded in rr_lens so median_length and
        // the Tukey RR rejection see the true interval.
        const int64_t rr_cap = (fs > 0.0) ? static_cast<int64_t>(3.0 * fs) : 0;

        // ---- slice + per-beat peak/foot --------------------------------
        for (size_t i = 0; i + 1 < rPeaks.size(); ++i) {
            const int64_t r0 = static_cast<int64_t>(rPeaks[i]);
            const int64_t rr = static_cast<int64_t>(rPeaks[i + 1]) - r0;
            if (rr <= 3) continue;

            const int64_t rr_w = (rr_cap > 0) ? std::min(rr, rr_cap) : rr;
            const int64_t before = rr_before_samples(rr_w);
            const int64_t after = rr_after_samples(rr_w);
            const int64_t len = before + after;
            const int64_t start = r0 - before;
            const int64_t end = r0 + after;

            std::vector<double> beat(static_cast<size_t>(len),
                std::numeric_limits<double>::quiet_NaN());
            const int64_t cs = std::max<int64_t>(0, start);
            const int64_t ce = std::min<int64_t>(N, end);
            for (int64_t k = cs; k < ce; ++k)
                beat[static_cast<size_t>(k - start)] = signal[static_cast<size_t>(k)];

            const int r_col = static_cast<int>(before);
            // Peak search is bounded to THIS beat's own R-R window: from
            // r_col+1 to the next R (which sits at r_col + rr). The full
            // beat slice extends past next R (beat length = 1.8*rr, next R
            // at 1.3*rr), so argmax-over-whole-beat would easily land on
            // the NEXT beat's peak whenever it's taller. That mislocated
            // "peak" then anchors up50 detection on the next beat's
            // upstroke, and up50-alignment then shifts every beat such
            // that individual beats' data effectively ends at their own
            // peak in the shared frame -- producing the peak-cutoff
            // plummet in the displayed template.
            const int peakSearchEnd = std::min(
                static_cast<int>(r_col + rr_w),
                static_cast<int>(beat.size()));
            int peak = -1;
            double pv = -std::numeric_limits<double>::infinity();
            for (int k = r_col + 1; k < peakSearchEnd; ++k)
                if (!std::isnan(beat[k]) && beat[k] > pv) { pv = beat[k]; peak = k; }
            if (peak < 0) continue;

            int foot = 0;
            double fv = std::numeric_limits<double>::infinity();
            for (int k = 0; k <= peak; ++k)
                if (!std::isnan(beat[k]) && beat[k] < fv) { fv = beat[k]; foot = k; }

            // 50%-upslope (half-height) crossing on [foot, peak]. This is the
            // horizontal alignment fiducial: it sits on the steep upstroke,
            // so its column is well-localized (unlike the flat apex or the
            // shallow foot). Linear-interpolate the crossing, round to a
            // column for the integer shift.
            const double half = fv + 0.5 * (pv - fv);
            int up50 = -1;
            for (int k = foot + 1; k <= peak; ++k) {
                if (std::isnan(beat[k]) || std::isnan(beat[k - 1])) continue;
                if (beat[k] >= half && beat[k - 1] < half) {
                    const double denom = beat[k] - beat[k - 1];
                    const double frac = (denom != 0.0)
                        ? (half - beat[k - 1]) / denom : 0.0;
                    up50 = static_cast<int>(std::lround((k - 1) + frac));
                    break;
                }
            }
            if (up50 < 0) continue;   // no clean upslope crossing; drop beat

            // NOTE: peak/foot are stored raw (no subsample_refine call) --
            // they only feed the peak-position rejection check below and
            // the (unread by any caller) peak_cols/foot_cols output,
            // neither of which needs sub-sample precision. The expensive
            // 4x-upsample fit-and-select refinement was pure overhead here;
            // the real fiducials used downstream (outPeakCol/outFootCol)
            // are recomputed independently from the final median template.
            raw.push_back({ std::move(beat), peak, foot, up50 });
            rr_lens.push_back(static_cast<int>(rr));
        }
        if (raw.empty()) return out;

        // ---- representative (median) length ----------------------------
        // Middle element of the sorted lengths (a length some beat has).
        {
            std::vector<int> lens = rr_lens;
            std::sort(lens.begin(), lens.end());
            out.median_length = lens[lens.size() / 2];
        }
        for (size_t i = 0; i < raw.size(); ++i) {
            if (rr_lens[i] == out.median_length) {
                out.ref_beat_index = static_cast<int>(i);
                break;
            }
        }

        // ---- Pass 1: 50%-upslope align on a shared axis (shift + NaN) ---
        // Anchor every beat's half-height point to the MEDIAN up50 column
        // among survivors. Beats with up50 > anchor get their leading samples
        // clipped (dst < 0 dropped by the guard below); this is accepted.
        std::vector<int> up50s;
        up50s.reserve(raw.size());
        for (const auto& r : raw) up50s.push_back(r.up50);
        std::sort(up50s.begin(), up50s.end());
        const int up50_anchor = up50s[up50s.size() / 2];   // median

        int max_tail = 0;         // max (beat_len - up50) over survivors
        for (const auto& r : raw) {
            const int tail = static_cast<int>(r.data.size()) - r.up50;
            if (tail > max_tail) max_tail = tail;
        }
        const int shared_w = up50_anchor + max_tail;
        if (shared_w <= 0) return out;
        const double NaND = std::numeric_limits<double>::quiet_NaN();

        out.beats.reserve(raw.size());
        out.peak_cols.reserve(raw.size());
        out.foot_cols.reserve(raw.size());
        for (const auto& b : raw) {
            const int prepend = up50_anchor - b.up50;   // may be < 0 now
            std::vector<double> a(shared_w, NaND);
            for (int k = 0; k < (int)b.data.size(); ++k) {
                const int dst = prepend + k;
                if (dst >= 0 && dst < shared_w) a[dst] = b.data[k];   // clips left overflow
            }
            out.beats.push_back(std::move(a));
            out.peak_cols.push_back(prepend + b.peak);   // may be < 0 for clipped beats
            out.foot_cols.push_back(prepend + b.foot);   // may be < 0 for clipped beats
        }
        out.up50_aligned_col = up50_anchor;

        // ---- Rejection: peak-column distance from the median peak column,
        // measured HERE (after alignment) rather than in each beat's own
        // local frame, because peak/foot/r_col all scale with that beat's
        // own RR before alignment -- raw column values aren't comparable
        // across beats until they share this common up50-anchored axis.
        //
        // The first 100 beats (or all, if fewer are available) seed the
        // reference median peak column and are kept unconditionally --
        // there's no reference to reject them against yet. Every beat after
        // that is rejected once its peak sits 5% of its own RR or further
        // from that median column. This guarantees at least
        // min(100, total_beats) beats always survive.
        {
            const size_t nSeed = std::min<size_t>(100, out.peak_cols.size());
            std::vector<int> seedPeaks(out.peak_cols.begin(), out.peak_cols.begin() + nSeed);
            std::sort(seedPeaks.begin(), seedPeaks.end());
            const int medianPeakCol = seedPeaks[seedPeaks.size() / 2];

            std::vector<bool> keep(out.peak_cols.size(), true);
            for (size_t i = nSeed; i < out.peak_cols.size(); ++i) {
                const double dist = std::abs(out.peak_cols[i] - medianPeakCol);
                keep[i] = (rr_lens[i] > 0) && (dist < 0.05 * rr_lens[i]);
            }

            std::vector<std::vector<double>> fBeats;
            std::vector<int> fPeak, fFoot, fRr;
            fBeats.reserve(out.beats.size());
            fPeak.reserve(out.beats.size());
            fFoot.reserve(out.beats.size());
            fRr.reserve(out.beats.size());
            for (size_t i = 0; i < out.beats.size(); ++i) {
                if (!keep[i]) continue;
                fBeats.push_back(std::move(out.beats[i]));
                fPeak.push_back(out.peak_cols[i]);
                fFoot.push_back(out.foot_cols[i]);
                fRr.push_back(rr_lens[i]);
            }
            out.beats.swap(fBeats);
            out.peak_cols.swap(fPeak);
            out.foot_cols.swap(fFoot);
            rr_lens.swap(fRr);

            // Re-anchor ref_beat_index to the filtered set (same
            // median_length target chosen above); fall back to the closest
            // surviving RR if the original reference beat itself got cut.
            out.ref_beat_index = -1;
            for (size_t i = 0; i < rr_lens.size(); ++i)
                if (rr_lens[i] == out.median_length) { out.ref_beat_index = static_cast<int>(i); break; }
            if (out.ref_beat_index < 0 && !rr_lens.empty()) {
                size_t best = 0;
                int bestDiff = std::abs(rr_lens[0] - out.median_length);
                for (size_t i = 1; i < rr_lens.size(); ++i) {
                    const int d = std::abs(rr_lens[i] - out.median_length);
                    if (d < bestDiff) { bestDiff = d; best = i; }
                }
                out.ref_beat_index = static_cast<int>(best);
            }
        }
        out.total_beats = out.beats.size();
        if (out.beats.empty()) return out;

        // ---- Pass 2: min-baseline vertical DC match ---------------------
        // Match each beat's baseline (windowed mean around its OWN minimum
        // sample) to the reference beat's. We use the per-beat argmin rather
        // than the stored foot column: after median-anchoring, a beat's foot
        // can be clipped off the left edge (foot_cols[i] < 0), but the argmin
        // over surviving samples is always a valid in-range column. For
        // unclipped beats the min IS the foot; for clipped beats it's the
        // lowest surviving point, a good-enough baseline proxy. Windowed mean
        // (not the single argmin sample) avoids order-statistic bias. A
        // constant vertical shift can't disturb horizontal alignment.
        if (out.ref_beat_index >= 0
            && out.ref_beat_index < static_cast<int>(out.beats.size()))
        {
            const int fb_w = std::max(1, out.median_length / 50);
            auto min_baseline = [&](size_t i) -> double {
                const auto& beat = out.beats[i];
                int mc = -1;
                double mv = std::numeric_limits<double>::infinity();
                for (int k = 0; k < (int)beat.size(); ++k) {
                    const double v = beat[k];
                    if (!std::isnan(v) && v < mv) { mv = v; mc = k; }
                }
                if (mc < 0) return std::numeric_limits<double>::quiet_NaN();
                const int lo = std::max(0, mc - fb_w);
                const int hi = std::min(shared_w, mc + fb_w + 1);
                double sum = 0.0; int n = 0;
                for (int k = lo; k < hi; ++k) {
                    const double v = beat[k];
                    if (!std::isnan(v)) { sum += v; ++n; }
                }
                return n >= 1 ? sum / n
                    : std::numeric_limits<double>::quiet_NaN();
                };

            const double target = min_baseline(static_cast<size_t>(out.ref_beat_index));
            if (!std::isnan(target)) {
                for (size_t i = 0; i < out.beats.size(); ++i) {
                    const double b_base = min_baseline(i);
                    if (std::isnan(b_base)) continue;
                    const double d = b_base - target;
                    if (d == 0.0) continue;
                    for (double& v : out.beats[i])
                        if (!std::isnan(v)) v -= d;
                }
            }
        }

        return out;
    }

}   // namespace alignment