#pragma once
//
// alignment.hpp
//
// Per-bin beat alignment for ECG and PPG. Runs BEFORE any normalization
// or template averaging.
//
//   1) Slice one beat per RR window as
//        [R_i - percent_interval_preceeding_rpeak*RR,
//         R_i + percent_interval_following_rpeak*RR].
//   2) Tukey outlier rejection: length (1.5), systolic amplitude (1.5),
//      and fiducial-position (3.0, PPG only) using [Q1-k*IQR, Q3+k*IQR].
//   3) Horizontal align onto a shared axis (integer shift + NaN pad, no
//      resampling):
//        - ECG: anchor R at percent_interval_preceeding_rpeak * max_rr_len.
//        - PPG: anchor the 50%-upslope (half-height) point at the largest
//               such column among survivors. The half-height point sits on
//               the steep systolic upstroke, so it is well-localized (like
//               ECG's R) and pinning it does NOT create an apex cusp; peaks
//               scatter and the apex averages honestly.
//   4) Vertical DC shift:
//        - ECG: match each beat's PR-baseline mean to the reference beat's.
//        - PPG: match each beat's foot-baseline mean (a small window around
//               its own foot column) to the reference beat's.
//

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <utility>
#include <vector>
namespace alignment {
    constexpr double percent_interval_preceeding_rpeak = 0.4; //how far before the R peak the snip goes, in terms of percent of the RR interval length
    constexpr double percent_interval_following_rpeak = 1.3;   //how far after the R peak the snip goes, in terms of percent of the RR interval length

    // NOTE: the second (Q-aligned) template pass no longer runs through this
    // function. It reuses the R-aligned + DC-aligned beats cached during the
    // R pass; the Q pass reuses these templates and shifts ECG to a common Q
    // column (see find_q_column below and qAlignTemplatesFromCache).
    // The old global g_q_align switch and the Pass-4 block are gone.

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
        size_t median_count = 0;
        size_t total_beats = 0;
        int r_aligned_col = -1;
        int q_aligned_col = -1;
        int ref_beat_index = -1;
    };

    inline ecg_beat_set extract_beats_and_align(const std::vector<double>& signal,
        const std::vector<size_t>& rPeaks)
    {
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
            for (size_t i = 0; i < keep.size(); ++i) {
                if (!keep[i]) continue;
                kb.push_back(std::move(out.beats[i]));
                kr.push_back(out.r_indices[i]);
                km.push_back(out.rr_lens[i]);
            }
            out.beats = std::move(kb);
            out.r_indices = std::move(kr);
            out.rr_lens = std::move(km);
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

        // ---- representative (median) length ----------------------------
        // Use the MEDIAN RR length as the representative beat length (not the
        // modal one). Take the middle element of the sorted lengths so the
        // value is guaranteed to be a length that some beat actually has,
        // which keeps the exact-match shape clustering below non-empty.
        {
            std::vector<int> lens = out.rr_lens;
            std::sort(lens.begin(), lens.end());
            out.median_length = lens[lens.size() / 2];
            out.median_count = 0;
            for (int L : out.rr_lens)
                if (L == out.median_length) ++out.median_count;
        }

        // ---- reference beat = first median-length beat -----------------
        // Used only as the PR-baseline DC-alignment reference. No RMS shape
        // clustering: the final template is the column-wise median.
        for (size_t i = 0; i < out.beats.size(); ++i)
            if (out.rr_lens[i] == out.median_length) {
                out.ref_beat_index = static_cast<int>(i);
                break;
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

        // ---- Pass 3: PR-baseline vertical DC alignment -----------------
        if (out.ref_beat_index >= 0
            && out.ref_beat_index < static_cast<int>(out.beats.size()))
        {
            const int pr_w = std::max(3, out.median_length / 20);
            const int pr_gap = std::max(1, out.median_length / 50);
            const int pr_lo = std::max(0, R_anchor - pr_gap - pr_w);
            const int pr_hi = std::max(pr_lo, R_anchor - pr_gap);

            auto pr_baseline = [&](const std::vector<double>& beat)
                -> std::pair<double, int>
                {
                    double sum = 0.0; int n = 0;
                    for (int k = pr_lo; k < pr_hi
                        && k < static_cast<int>(beat.size()); ++k) {
                        const double v = beat[k];
                        if (!std::isnan(v)) { sum += v; ++n; }
                    }
                    if (n < 3)
                        return { std::numeric_limits<double>::quiet_NaN(), n };
                    return { sum / n, n };
                };

            const auto [target, target_n] =
                pr_baseline(out.beats[out.ref_beat_index]);

            if (!std::isnan(target)) {
                for (auto& beat : out.beats) {
                    const auto [b_base, b_n] = pr_baseline(beat);
                    if (std::isnan(b_base)) continue;
                    const double d = b_base - target;
                    if (d == 0.0) continue;
                    for (double& v : beat)
                        if (!std::isnan(v)) v -= d;
                }
            }
        }

        // Q-align is no longer done here. extract_beats_and_align now only
        // R-aligns (+ DC-aligns); the Q pass reuses these beats via
        // qAlignTemplatesFromCache in build_templates.hpp (reuse + Q-shift).
        out.q_aligned_col = out.r_aligned_col;
        return out;
    }

    // =========================================================================
    // Locate the Q trough of an ECG beat/template: steepest R-upstroke, walk
    // back to its onset, then the first local minimum before it. R_anchor is
    // the R column; fs is the ECG rate. Returns -1 if no Q can be located.
    inline int find_q_column(const std::vector<double>& b, int R_anchor, double fs) {
        const int Wsh = static_cast<int>(b.size());
        if (Wsh == 0 || R_anchor < 1 || fs <= 0.0) return -1;

        int finite = 0;
        for (double v : b) if (!std::isnan(v)) ++finite;
        const int w20 = std::max(2, static_cast<int>(std::lround(0.020 * fs)));
        const int fb = std::max(1, static_cast<int>(
            std::lround(0.10 * (finite > 0 ? finite : Wsh))));

        double maxSlope = 0.0; int steep = R_anchor;
        for (int k = std::max(1, R_anchor - 3 * w20); k < R_anchor && k + 1 < Wsh; ++k) {
            if (std::isnan(b[k]) || std::isnan(b[k + 1])) continue;
            const double s = b[k + 1] - b[k];
            if (s > maxSlope) { maxSlope = s; steep = k; }
        }
        int onset = steep;
        if (maxSlope > 0.0) {
            for (int k = steep; k > std::max(1, R_anchor - 4 * w20); --k) {
                if (std::isnan(b[k]) || std::isnan(b[k - 1])) break;
                onset = k;
                if ((b[k] - b[k - 1]) < 0.2 * maxSlope) break;
            }
        }
        const int lo = std::max(1, onset - 2 * w20);
        for (int k = onset - 1; k > lo; --k) {
            if (std::isnan(b[k - 1]) || std::isnan(b[k]) || std::isnan(b[k + 1])) continue;
            if (b[k] <= b[k - 1] && b[k] <= b[k + 1]) return k;
        }
        const int q = onset - fb;
        return (q >= 0 && q < Wsh) ? q : -1;
    }

    // =========================================================================
    // Q-align a bin's cached (R-aligned) snippets and re-median.
    //
    // Mirrors R-alignment's structure: R-align picks a reference "median
    // snippet" (the first beat whose finite length equals the median) and
    // aligns to it. Here we do the same, but on Q. The alignment reference is
    // the column-wise MEDIAN of all the beats over the -0.25..+0.75 RR window
    // around R -- which is exactly the R-aligned template (refMedian), since
    // that template is that column median and find_q only looks in the
    // R-upstroke region (inside that window). We find Q on that median, then
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

    inline QAlignResult q_align_beat_matrix(
        const std::vector<std::vector<double>>& beatsIn,
        int R_anchor, double fs, bool compute_iqr,
        const std::vector<double>& refMedian)
    {
        QAlignResult res;
        if (beatsIn.empty()) return res;

        std::vector<std::vector<double>> beats = beatsIn;
        const int Wsh = static_cast<int>(beats.front().size());

        // Q marker = Q of the column-wise median of all beats (== the R
        // template). Every snippet is then shifted so its own Q lands here.
        const int q_marker = find_q_column(refMedian, R_anchor, fs);
        res.q_aligned_col = (q_marker >= 0) ? q_marker : R_anchor;

        std::vector<int> shifts;   // per-beat Q-align shift (for R's new column)
        if (q_marker >= 0 && fs > 0.0) {
            const double NaNv = std::numeric_limits<double>::quiet_NaN();
            shifts.reserve(beats.size());
            for (size_t i = 0; i < beats.size(); ++i) {
                const int qi = find_q_column(beats[i], R_anchor, fs);
                if (qi < 0) continue;
                const int shift = q_marker - qi;   // may be negative
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
        size_t median_count = 0;
        size_t total_beats = 0;
        int    up50_aligned_col = -1;             // shared column all half-height points land on
        int    foot_aligned_col = -1;             // feet scatter (per-beat); not a shared column
        int    peak_aligned_col = -1;             // peaks scatter (per-beat); not a shared column
        int    ref_beat_index = -1;
    };

    inline PpgBeatSet extract_ppg_beats_and_align(
        const std::vector<double>& signal,
        const std::vector<size_t>& rPeaks)
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

        // ---- slice + per-beat peak/foot --------------------------------
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

            const int r_col = static_cast<int>(before);
            int peak = -1;
            double pv = -std::numeric_limits<double>::infinity();
            for (int k = r_col + 1; k < (int)beat.size(); ++k)
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

            raw.push_back({ std::move(beat), peak, foot, up50 });
            rr_lens.push_back(static_cast<int>(rr));
        }
        if (raw.empty()) return out;

        // ---- Tukey rejection: length (1.5*IQR) --------------------------
        {
            std::vector<double> lens_d(rr_lens.begin(), rr_lens.end());
            apply_mask(keep_within_tukey(lens_d, 1.5));
        }
        if (raw.empty()) return out;

        // ---- Tukey rejection: systolic amplitude peak-foot (1.5*IQR) ----
        {
            std::vector<double> amps;
            amps.reserve(raw.size());
            for (const auto& r : raw) amps.push_back(r.data[r.peak] - r.data[r.foot]);
            apply_mask(keep_within_tukey(amps, 1.5));
        }
        if (raw.empty()) return out;

        // ---- Tukey rejection: 50%-upslope position (1.5*IQR) ------------
        // up50 is the horizontal alignment fiducial and is searched per beat
        // (unlike ECG's analytic R column), so we guard its position here.
        {
            std::vector<double> poss;
            poss.reserve(raw.size());
            for (const auto& r : raw) poss.push_back(static_cast<double>(r.up50));
            apply_mask(keep_within_tukey(poss, 1.5));
        }
        out.total_beats = raw.size();
        if (raw.empty()) return out;

        // ---- representative (median) length ----------------------------
        // Middle element of the sorted lengths (a length some beat has).
        {
            std::vector<int> lens = rr_lens;
            std::sort(lens.begin(), lens.end());
            out.median_length = lens[lens.size() / 2];
            out.median_count = 0;
            for (int L : rr_lens) if (L == out.median_length) ++out.median_count;
        }
        for (size_t i = 0; i < raw.size(); ++i) {
            if (rr_lens[i] == out.median_length) {
                out.ref_beat_index = static_cast<int>(i);
                break;
            }
        }

        // ---- Pass 1: 50%-upslope align on a shared axis (shift + NaN) ---
        // Anchor every beat's half-height point to the largest up50 column
        // among survivors. Anchoring on the max fiducial column (as ECG does
        // with 0.3*max) guarantees prepend >= 0, so beats are only ever
        // NaN-prepended, never left-clipped. The half-height point sits on the
        // steep upstroke, so its column is well-localized and pinning it does
        // not manufacture an apex cusp; peaks and feet both scatter.
        int up50_anchor = 0;      // max up50 column over survivors
        int max_tail = 0;         // max (beat_len - up50) over survivors
        for (const auto& r : raw) {
            if (r.up50 > up50_anchor) up50_anchor = r.up50;
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
            const int prepend = up50_anchor - b.up50;   // >= 0 by construction
            std::vector<double> a(shared_w, NaND);
            for (int k = 0; k < (int)b.data.size(); ++k) {
                const int dst = prepend + k;
                if (dst >= 0 && dst < shared_w) a[dst] = b.data[k];
            }
            out.beats.push_back(std::move(a));
            out.peak_cols.push_back(prepend + b.peak);   // scatters across beats
            out.foot_cols.push_back(prepend + b.foot);   // scatters across beats
        }
        out.up50_aligned_col = up50_anchor;

        // ---- Pass 2: foot-baseline vertical DC match (mirrors ECG) ------
        // Match each beat's mean over a small window around its OWN foot column
        // to the reference beat's. Feet no longer share a column (we aligned on the
        // upslope), so the window is centered per beat at foot_cols[i]. A
        // windowed mean, rather than the single argmin sample, avoids an
        // order-statistic bias at the foot. A constant vertical shift moves no
        // column, so it cannot affect the (cusp-free) horizontal alignment.
        if (out.ref_beat_index >= 0
            && out.ref_beat_index < static_cast<int>(out.beats.size()))
        {
            const int fb_w = std::max(1, out.median_length / 50);
            auto foot_baseline = [&](size_t i) -> double {
                const int fc = out.foot_cols[i];
                const int lo = std::max(0, fc - fb_w);
                const int hi = std::min(shared_w, fc + fb_w + 1);
                double sum = 0.0; int n = 0;
                for (int k = lo; k < hi; ++k) {
                    const double v = out.beats[i][k];
                    if (!std::isnan(v)) { sum += v; ++n; }
                }
                return n >= 1 ? sum / n
                    : std::numeric_limits<double>::quiet_NaN();
                };

            const double target = foot_baseline(static_cast<size_t>(out.ref_beat_index));
            if (!std::isnan(target)) {
                for (size_t i = 0; i < out.beats.size(); ++i) {
                    const double b_base = foot_baseline(i);
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