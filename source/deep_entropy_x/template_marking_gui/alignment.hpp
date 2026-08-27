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
        // Parallel to `beats`: the Section 4.6 rhythm verdict, assigned
        // AFTER SLICING AND BEFORE ANY PRUNING, and compacted through every
        // apply_mask() alongside `beats`.
        //
        // Why here and not downstream. The first Tukey pass below rejects on
        // RR LENGTH, and a premature beat is short by definition -- so
        // alignment itself throws ectopics out as length outliers, for outlier
        // reasons, with no flag and no record. Anything that asks "was this
        // beat premature" after that pass is asking about the beats that
        // survived not being premature. Assigned here, the verdict exists for
        // every sliced beat and survives whatever the pruning does next.
        //
        // premature: RR(i) < 0.80 * median of the trailing ten.
        // voted:     the 5-of-8 rule over the raw premature flags -- the
        //            middle of a run, where the trailing median has itself
        //            gone short so the beat stops reading as premature alone.
        std::vector<char> premature;
        std::vector<char> voted;
        // Counts as they stood BEFORE pruning, so the flag survivors can be
        // compared against what the record actually contained.
        int n_premature_presliced = 0;
        int n_voted_presliced = 0;

        // Parallel to `beats`: PQ_level - TP_level whenever BOTH stage
        // estimates succeeded for that beat (NaN otherwise). QC metric per
        // spec I-1; not used to drive any decision, just recorded.
        std::vector<double> tp_pq_delta;
        // Parallel to `beats`: the per-stage vertical DC shifts applied by the
        // two-stage leveling. tp_shift = Stage 1 (TP) amount; pq_shift =
        // Stage 2 (PQ finalize) amount. NaN where that stage didn't apply.
        // Returned for the per-beat move log.
        std::vector<double> tp_shift;
        std::vector<double> pq_shift;
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
            std::vector<char> kp, kv;
            const bool haveSrc = out.baseline_source.size() == out.beats.size();
            const bool haveDelta = out.tp_pq_delta.size() == out.beats.size();
            const bool haveFlags = out.premature.size() == out.beats.size()
                && out.voted.size() == out.beats.size();
            for (size_t i = 0; i < keep.size(); ++i) {
                // A RHYTHM-FLAGGED BEAT IS NEVER PRUNED.
                //
                // The first Tukey pass below rejects on RR LENGTH at 1.5*IQR,
                // and a premature beat is short by definition -- so without
                // this exemption alignment discards the ectopy as a length
                // outlier, for outlier reasons, and every beat that survives
                // is one that was not premature. Measured on a record with 9
                // scripted PVCs: all 9 were pruned, the flag vector came out
                // empty of positives, and the per-beat output read NORMAL
                // throughout while being entirely correct about the beats it
                // still had.
                //
                // 4.6 requires these beats "excluded from the reference
                // template but RETAINED with flags". Excluding them is the
                // ectopic mask's job (create_ecg_templates.hpp), and it keeps
                // a record; dropping them here leaves nothing to retain.
                // Pruning is for detector errors, not for real ectopy.
                const bool flagged = haveFlags
                    && (out.premature[i] || out.voted[i]);
                if (!keep[i] && !flagged) continue;   // rejected, discarded
                kb.push_back(std::move(out.beats[i]));
                kr.push_back(out.r_indices[i]);
                km.push_back(out.rr_lens[i]);
                if (haveSrc) ks.push_back(out.baseline_source[i]);
                if (haveDelta) kd.push_back(out.tp_pq_delta[i]);
                if (haveFlags) { kp.push_back(out.premature[i]); kv.push_back(out.voted[i]); }
            }
            out.beats = std::move(kb);
            out.r_indices = std::move(kr);
            out.rr_lens = std::move(km);
            if (haveSrc) out.baseline_source = std::move(ks);
            if (haveDelta) out.tp_pq_delta = std::move(kd);
            if (haveFlags) { out.premature = std::move(kp); out.voted = std::move(kv); }
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

        // ---- Section 4.6 rhythm verdict: after the slice, before the -----
        // ---- pruning. Indexed with `beats`; compacted with them.      -----
        // out.rr_lens[i] is beat i's own interval, so this is exact -- one
        // interval per beat, none missing, which is precisely what stops being
        // true the moment the Tukey passes below run.
        {
            const size_t nb = out.beats.size();
            out.premature.assign(nb, 0);
            out.voted.assign(nb, 0);
            if (nb >= 12) {
                // Beat t is premature when the interval BEFORE it is short.
                // rr_lens[i] holds the interval AFTER beat i (R[i+1] - R[i]),
                // because that is the span the slice covers -- so beat t's
                // PRECEDING interval is rr_lens[t-1]. Testing rr_lens[t]
                // directly flags the beat that precedes each PVC instead of
                // the PVC, which is an off-by-one that looks exactly like a
                // detector that "nearly works".
                for (size_t t = 11; t < nb; ++t) {
                    std::vector<double> w(out.rr_lens.begin() + (t - 11),
                        out.rr_lens.begin() + (t - 1));   // the ten before it
                    std::sort(w.begin(), w.end());
                    const double med = w[w.size() / 2];
                    if (med > 0.0 && out.rr_lens[t - 1] < 0.80 * med)
                        out.premature[t] = 1;
                }
                // The vote reads the RAW flags, never its own output: feeding
                // it back would let one run grow along the whole record.
                for (size_t t = 0; t < nb; ++t) {
                    const size_t lo = (t > 4) ? t - 4 : 0;
                    const size_t hi = std::min(nb, t + 4);
                    int c = 0;
                    for (size_t i = lo; i < hi; ++i) c += out.premature[i];
                    if (c >= 5) out.voted[t] = 1;
                }
                for (size_t i = 0; i < nb; ++i) {
                    out.n_premature_presliced += out.premature[i];
                    if (out.voted[i] && !out.premature[i]) ++out.n_voted_presliced;
                }
            }
        }

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
            // RR-fraction bounds for the TP baseline window, and ms-before-R
            // bounds for the PQ baseline window. Conservative defaults; tune
            // against real recordings if the leveling looks off.
            const double kTpLoFrac = 0.55;   // start of TP flat, fraction of RR after R
            const double kTpHiFrac = 0.85;   // end of TP flat, before next P
            const double kPqPreRMs = 80.0;   // PQ window opens this many ms before R
            const double kPqGuardMs = 20.0;  // and closes this many ms before R (guard vs Q)

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

            /* Stage 1: TP Window. The baseline segment is the TP window ESTIMATED
            * by range:
            * lo = R_anchor + 0.55 · RR
            * hi = R_anchor + 0.85 · RR
            * And then the flattest 10ms wide segment is located via a minimum variance scan
            */
            auto tp_window = [&](size_t i) -> std::pair<int, int> {
                const int N = static_cast<int>(out.beats[i].size());
                const int rr = out.rr_lens[i];
                if (rr <= 0) return { -1, -1 };
                const int lo = R_anchor + static_cast<int>(std::lround(kTpLoFrac * rr));
                const int hi = R_anchor + static_cast<int>(std::lround(kTpHiFrac * rr));
                const int clo = std::max(0, lo);
                const int chi = std::min(N, hi);
                if (chi - clo < 3) return { -1, -1 };
                return { clo, chi };
                };

            // Stage 2, PQ window:Segment is always 80 ms before R to 20 ms before R
            auto pq_window = [&](size_t i) -> std::pair<int, int> {
                const int N = static_cast<int>(out.beats[i].size());
                const int lo = R_anchor - static_cast<int>(std::lround(kPqPreRMs * 0.001 * fs));
                const int hi = R_anchor - static_cast<int>(std::lround(kPqGuardMs * 0.001 * fs));
                const int clo = std::max(0, lo);
                const int chi = std::min(N, hi);
                if (chi - clo < 3) return { -1, -1 };
                return { clo, chi };
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

            const auto [refTp, refTpOk, refPq, refPqOk] = both_levels(static_cast<size_t>(out.ref_beat_index));
            BaselineSource refSrc = BaselineSource::NONE;

            /*Two-stage vertical alignment. - Stage 1: level each beat on the TP segment (end of T to just before onset of
            next P). Estimate level with median over the flattest sub-window bounded by fitted landmarks. - Stage
            2: finalize the zero on the PQ segment (end of P to immediately before Q onset). PQ is the higher-priority reference.*/
            double refTpTarget = refTpOk ? refTp : std::numeric_limits<double>::quiet_NaN();
            double refPqTarget = refPqOk ? refPq : std::numeric_limits<double>::quiet_NaN();
            if (std::isnan(refTpTarget) || std::isnan(refPqTarget)) {
                for (size_t i = 0; i < out.beats.size(); ++i) {
                    if (!std::isnan(refTpTarget) && !std::isnan(refPqTarget)) break;
                    const auto [tTp, tTpOk, tPq, tPqOk] = both_levels(i);
                    if (tTpOk && std::isnan(refTpTarget)) refTpTarget = tTp;
                    if (tPqOk && std::isnan(refPqTarget)) refPqTarget = tPq;
                }
            }
            const bool haveAnyRef = !std::isnan(refTpTarget) || !std::isnan(refPqTarget);

            out.tp_pq_delta.assign(out.beats.size(), std::numeric_limits<double>::quiet_NaN());

            // Diagnostic accumulators (spec-neutral: recorded only, doesn't
            // change any beat).
            int diag_pq = 0, diag_tp = 0, diag_none = 0;
            std::vector<double> diag_shifts;
            diag_shifts.reserve(out.beats.size());

            out.tp_shift.assign(out.beats.size(), std::numeric_limits<double>::quiet_NaN());
            out.pq_shift.assign(out.beats.size(), std::numeric_limits<double>::quiet_NaN());
            if (haveAnyRef) {
                for (size_t i = 0; i < out.beats.size(); ++i) {
                    const auto [tpLvl, tpOk, pqLvl, pqOk] = both_levels(i);
                    if (tpOk && pqOk) out.tp_pq_delta[i] = pqLvl - tpLvl;

                    // Stage 1 (TP): coarse level, applied only if BOTH this
                    // beat and the reference have a usable TP.
                    double applied = 0.0;
                    bool leveled = false;
                    if (tpOk && !std::isnan(refTpTarget)) {
                        const double d = tpLvl - refTpTarget;
                        for (double& v : out.beats[i]) if (!std::isnan(v)) v -= d;
                        applied += d; leveled = true;
                        out.tp_shift[i] = d;
                    }
                    // Stage 2 (PQ finalize): authoritative zero. Recompute PQ
                    // level AFTER stage 1 shift, match it to the reference PQ.
                    if (!std::isnan(refPqTarget)) {
                        const auto [plo, phi] = pq_window(i);
                        if (plo >= 0) {
                            const auto [lvl2, ok2] = flattest_median(out.beats[i], plo, phi);
                            if (ok2) {
                                const double d = lvl2 - refPqTarget;
                                for (double& v : out.beats[i]) if (!std::isnan(v)) v -= d;
                                applied += d; leveled = true;
                                out.pq_shift[i] = d;
                                out.baseline_source[i] = BaselineSource::PQ;
                            }
                        }
                    }
                    if (!leveled) { out.baseline_source[i] = BaselineSource::NONE; continue; }
                    if (out.baseline_source[i] != BaselineSource::PQ)
                        out.baseline_source[i] = BaselineSource::TP;   // TP-only
                    diag_shifts.push_back(applied);
                    if (out.baseline_source[i] == BaselineSource::PQ) ++diag_pq;
                    else ++diag_tp;
                }
                diag_none = (int)out.beats.size() - diag_pq - diag_tp;
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
    struct aligned_beats {
        std::vector<double> tmpl;
        std::vector<double> iqr;
        std::vector<std::vector<double>> beats;
        int q_aligned_col = -1;
        int r_col = -1;
    };

    inline aligned_beats align_beat_matrix(
        const std::vector<std::vector<double>>& beatsIn,
        int R_anchor, double fs, bool compute_iqr,
        const std::vector<double>& ref_beat_of_median_length,
        const std::function<double(const std::vector<double>&)>& locate)
    {
        aligned_beats res;
        if (beatsIn.empty()) return res;

        std::vector<std::vector<double>> beats = beatsIn;
        const int Wsh = static_cast<int>(beats.front().size());

        // Q marker = Q of the column-wise median of all beats (== the R
        // template). Every snippet is then shifted so its own Q lands here.
        // AnchorLocator returns a sub-sample position now; this path shifts by
        // whole samples, so it rounds. Numerically identical to the previous
        // int-returning locator, which rounded inside the finder.
        const int marker = static_cast<int>(std::lround(
            locate(ref_beat_of_median_length)));
        res.q_aligned_col = (marker >= 0) ? marker : R_anchor;

        std::vector<int> shifts;   // per-beat Q-align shift (for R's new column)
        if (marker >= 0) {
            const double NaNv = std::numeric_limits<double>::quiet_NaN();
            shifts.reserve(beats.size());
            for (size_t i = 0; i < beats.size(); ++i) {
                const int mi = static_cast<int>(std::lround(locate(beats[i])));
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
                // Per-column standard deviation (matches the QC noise metric,
                // which also uses std). Computed over the same non-NaN values
                // used for the median. Field is still named 'iqr' downstream
                // (serializer/viewer) but now carries std.
                double mean = 0.0;
                for (double v : col) mean += v;
                mean /= static_cast<double>(nc);
                double ss = 0.0;
                for (double v : col) { const double d = v - mean; ss += d * d; }
                res.iqr[c] = std::sqrt(ss / static_cast<double>(nc - 1));  // sample std
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
        res.beats = std::move(beats);
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
            // Upstroke-located FIRST peak, not the tallest sample in the window.
            // This is the site that matters most: the peak found here brackets
            // the foot and the up50 half-height crossing below, and a beat whose
            // up50 cannot be found is DISCARDED. With argmax, a pulse whose
            // reflected wave exceeds systole had its up50 searched on the
            // notch-to-P2 rise, which frequently has no clean single crossing --
            // so those beats were dropped, and the surviving count collapsed.
            const int peak = FeatureMarks::detect_ppg_upstroke_peak(beat, r_col + 1,
                peakSearchEnd);
            if (peak < 0) continue;

            // Foot and up50 through the SAME primitives the display fiducials
            // use (FeatureMarks::trough_in / amplitude_crossing), so the
            // alignment axis and the markers drawn on it are defined
            // identically. Both were hand-rolled here, which is why a fix to
            // one never reached the other.
            const int foot = FeatureMarks::trough_in(beat, 0, peak);
            if (foot < 0) continue;

            // 50%-upslope (half-height) crossing on [foot, peak]. This is the
            // horizontal alignment fiducial: it sits on the steep upstroke, so
            // its column is well-localized (unlike the flat apex or the shallow
            // foot). Interpolated first-upward-crossing, and it can FAIL --
            // both properties matter here and neither is provided by
            // amplitude_crossing, which is for display markers.
            // first_crossing returns a sub-sample crossing now; the up50 axis
            // is integer-column, so it rounds here.
            const double up50D = FeatureMarks::first_crossing(beat, foot, peak, 0.50);
            const int up50 = (up50D >= 0.0)
                ? static_cast<int>(std::lround(up50D)) : -1;
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