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
        size_t total_beats = 0;
        int r_aligned_col = -1;
        int q_aligned_col = -1;
        int ref_beat_index = -1;
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

        // ---- Pass 3: PQ-baseline vertical DC alignment ------------------
        if (out.ref_beat_index >= 0 && out.ref_beat_index < static_cast<int>(out.beats.size()))
        {
            const int pq_min_w = std::max(1, static_cast<int>(std::lround(0.010 * fs)));  // 10 ms window
            const int pq_ceiling = static_cast<int>(std::lround(0.100 * fs));               // 100 ms max lookback

            const int pq_hi_bound = R_anchor;
            const int pq_lo_bound = std::max(0, R_anchor - pq_ceiling);

            auto find_flat_baseline = [&](const std::vector<double>& beat)
                -> std::pair<double, int>   // {mean, window_start}
                {
                    const int hi = std::min(pq_hi_bound,
                        static_cast<int>(beat.size()) - pq_min_w);
                    if (hi < pq_lo_bound)
                        return { std::numeric_limits<double>::quiet_NaN(), -1 };

                    double best_var = std::numeric_limits<double>::infinity();
                    double best_mean = std::numeric_limits<double>::quiet_NaN();
                    int best_start = -1;

                    for (int s = pq_lo_bound; s <= hi; ++s) {
                        double sum = 0.0, sumsq = 0.0; int n = 0;
                        for (int k = s; k < s + pq_min_w; ++k) {
                            const double v = beat[k];
                            if (std::isnan(v)) continue;
                            sum += v; sumsq += v * v; ++n;
                        }
                        if (n < static_cast<int>(pq_min_w * 0.7)) continue;
                        const double mean = sum / n;
                        const double var = sumsq / n - mean * mean;
                        if (var < best_var) { best_var = var; best_mean = mean; best_start = s; }
                    }
                    return { best_mean, best_start };
                };

            const auto [target, target_start] =
                find_flat_baseline(out.beats[out.ref_beat_index]);

            if (!std::isnan(target)) {
                for (auto& beat : out.beats) {
                    const auto [b_base, b_start] = find_flat_baseline(beat);
                    if (std::isnan(b_base)) continue;
                    const double d = b_base - target;
                    if (d == 0.0) continue;
                    for (double& v : beat)
                        if (!std::isnan(v)) v -= d;
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
        size_t total_beats = 0;
        int    up50_aligned_col = -1;             // shared column all half-height points land on
        int    foot_aligned_col = -1;             // feet scatter (per-beat); not a shared column
        int    peak_aligned_col = -1;             // peaks scatter (per-beat); not a shared column
        int    ref_beat_index = -1;
    };

    inline PpgBeatSet extract_ppg_beats_and_align(const std::vector<double>& signal, const std::vector<size_t>& rPeaks)
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