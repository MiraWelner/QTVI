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
//        - ECG: anchor R at percent_interval_preceeding_rpeak * max_mode_len.
//        - PPG: anchor the 50%-upslope (half-height) point at the largest
//               such column among survivors. The half-height point sits on
//               the steep systolic upstroke, so it is well-localized (like
//               ECG's R) and pinning it does NOT create an apex cusp; peaks
//               scatter and the apex averages honestly.
//   4) Vertical DC shift:
//        - ECG: match each beat's PR-baseline mean to the mode beat's.
//        - PPG: match each beat's foot-baseline mean (a small window around
//               its own foot column) to the mode beat's.
//

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>
namespace alignment {
    constexpr double percent_interval_preceeding_rpeak = 0.3; //how far before the R peak the snip goes, in terms of percent of the RR interval length
    constexpr double percent_interval_following_rpeak = 1.4;   //how far after the R peak the snip goes, in terms of percent of the RR interval length
	constexpr double mode_error = 0.05; // RMS tolerance for mode beat clustering (5% of the mode beat's amplitude)

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

    // =========================================================================
    // ECG
    // =========================================================================
    struct BeatSet {
        std::vector<std::vector<double>> beats;   // NaN-padded, R-aligned
        std::vector<size_t> r_indices;
        std::vector<int>    mode_lens;

        int    mode_length = -1;
        size_t mode_count = 0;
        size_t total_beats = 0;

        size_t mode_group_size = 0;
        std::vector<size_t> mode_group_indices;
        double mode_group_rms_tol = 0.0;

        int r_aligned_col = -1;
        int mode_beat_index = -1;
    };

    inline BeatSet extract_beats_and_mode(const std::vector<double>& signal,
        const std::vector<size_t>& rPeaks,
        double rms_tol = mode_error)
    {
        BeatSet out;
        const int64_t N = static_cast<int64_t>(signal.size());
        if (N == 0 || rPeaks.size() < 2) {
            fprintf(stderr, "[align] no beats: signal=%lld rPeaks=%zu\n",
                static_cast<long long>(N), rPeaks.size());
            return out;
        }

        // Compact the parallel (beats, r_indices, mode_lens) vectors, keeping
        // only entries where keep[i] is true.
        auto apply_mask = [&](const std::vector<bool>& keep) {
            std::vector<std::vector<double>> kb;
            std::vector<size_t> kr;
            std::vector<int>    km;
            for (size_t i = 0; i < keep.size(); ++i) {
                if (!keep[i]) continue;
                kb.push_back(std::move(out.beats[i]));
                kr.push_back(out.r_indices[i]);
                km.push_back(out.mode_lens[i]);
            }
            out.beats = std::move(kb);
            out.r_indices = std::move(kr);
            out.mode_lens = std::move(km);
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
            out.mode_lens.push_back(static_cast<int>(rr));
        }
        if (out.beats.empty()) return out;

        // ---- Tukey rejection: length (1.5*IQR) --------------------------
        {
            std::vector<double> lens_d(out.mode_lens.begin(), out.mode_lens.end());
            apply_mask(keep_within_tukey(lens_d, 1.5));
        }
        if (out.beats.empty()) return out;

        // ---- Tukey rejection: R amplitude (1.5*IQR) ---------------------
        // Position pass is skipped for ECG: R sits at column
        // rr_before_samples(mode_lens[i]) in every beat by construction, so
        // intra-beat position is fixed.
        {
            std::vector<double> amps;
            amps.reserve(out.beats.size());
            for (size_t i = 0; i < out.beats.size(); ++i) {
                const int r_col = static_cast<int>(rr_before_samples(out.mode_lens[i]));
                amps.push_back(
                    (r_col < (int)out.beats[i].size() && !std::isnan(out.beats[i][r_col]))
                    ? out.beats[i][r_col]
                    : std::numeric_limits<double>::quiet_NaN());
            }
            auto keepA = keep_within_tukey(amps, 1.5);
            for (size_t i = 0; i < keepA.size(); ++i)
                if (std::isnan(amps[i])) keepA[i] = false;
            apply_mask(keepA);
        }
        out.total_beats = out.beats.size();
        if (out.beats.empty()) return out;

        // ---- mode length -----------------------------------------------
        std::unordered_map<int, size_t> hist;
        hist.reserve(out.beats.size() * 2);
        for (int L : out.mode_lens) hist[L]++;
        for (const auto& kv : hist) {
            if (kv.second > out.mode_count
                || (kv.second == out.mode_count && kv.first > out.mode_length)) {
                out.mode_length = kv.first;
                out.mode_count = kv.second;
            }
        }

        // ---- shape cluster among mode-length beats ---------------------
        out.mode_group_rms_tol = rms_tol;
        if (out.mode_length > 0 && out.mode_count > 0) {
            std::vector<size_t> modeIdx;
            modeIdx.reserve(out.mode_count);
            for (size_t i = 0; i < out.beats.size(); ++i)
                if (out.mode_lens[i] == out.mode_length)
                    modeIdx.push_back(i);

            struct Cluster { size_t seed; std::vector<size_t> members; };
            std::vector<Cluster> clusters;

            const int L = out.mode_length;
            const double L_inv = 1.0 / static_cast<double>(L);
            auto rms_diff = [&](const std::vector<double>& a,
                const std::vector<double>& b) {
                    double ss = 0.0;
                    for (int k = 0; k < L; ++k) {
                        const double d = a[k] - b[k];
                        ss += d * d;
                    }
                    return std::sqrt(ss * L_inv);
                };

            for (size_t idx : modeIdx) {
                const auto& beat = out.beats[idx];
                bool placed = false;
                for (auto& c : clusters) {
                    if (rms_diff(beat, out.beats[c.seed]) <= rms_tol) {
                        c.members.push_back(idx);
                        placed = true;
                        break;
                    }
                }
                if (!placed) {
                    Cluster c;
                    c.seed = idx;
                    c.members.push_back(idx);
                    clusters.push_back(std::move(c));
                }
            }

            for (const auto& c : clusters) {
                if (c.members.size() > out.mode_group_size) {
                    out.mode_group_size = c.members.size();
                    out.mode_group_indices = c.members;
                }
            }
            if (!out.mode_group_indices.empty()) {
                out.mode_beat_index =
                    static_cast<int>(out.mode_group_indices.front());
            }
        }

        // ---- Pass 1: R-align on shared axis ----------------------------
        int max_mode_len = 0;
        for (int L : out.mode_lens)
            if (L > max_mode_len) max_mode_len = L;
        if (max_mode_len <= 0) return out;

        const int R_anchor = static_cast<int>(rr_before_samples(max_mode_len));
        const int shared_w = R_anchor + static_cast<int>(rr_after_samples(max_mode_len));

        const double NaND = std::numeric_limits<double>::quiet_NaN();
        std::vector<std::vector<double>> aligned;
        aligned.reserve(out.beats.size());
        for (size_t i = 0; i < out.beats.size(); ++i) {
            const auto& b = out.beats[i];
            const int L = static_cast<int>(b.size());
            const int r_in_beat = static_cast<int>(rr_before_samples(out.mode_lens[i]));
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
        if (out.mode_beat_index >= 0
            && out.mode_beat_index < static_cast<int>(out.beats.size()))
        {
            const int pr_w = std::max(3, out.mode_length / 20);
            const int pr_gap = std::max(1, out.mode_length / 50);
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
                pr_baseline(out.beats[out.mode_beat_index]);

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

        return out;
    }

    // =========================================================================
    // PPG
    // =========================================================================
    struct PpgBeatSet {
        std::vector<std::vector<double>> beats;   // NaN-padded, 50%-upslope-aligned
        std::vector<int> peak_cols;               // per-beat systolic peak column (varies)
        std::vector<int> foot_cols;               // per-beat foot column (varies)
        int    mode_length = -1;
        size_t mode_count = 0;
        size_t total_beats = 0;
        int    up50_aligned_col = -1;             // shared column all half-height points land on
        int    foot_aligned_col = -1;             // feet scatter (per-beat); not a shared column
        int    peak_aligned_col = -1;             // peaks scatter (per-beat); not a shared column
        int    mode_beat_index = -1;
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
        std::vector<int> mode_lens;
        raw.reserve(rPeaks.size());

        // Compact the parallel (raw, mode_lens) vectors, keeping only entries
        // where keep[i] is true.
        auto apply_mask = [&](const std::vector<bool>& keep) {
            std::vector<Raw> filt;
            std::vector<int> filt_lens;
            filt.reserve(raw.size());
            filt_lens.reserve(raw.size());
            for (size_t i = 0; i < raw.size(); ++i)
                if (keep[i]) {
                    filt.push_back(std::move(raw[i]));
                    filt_lens.push_back(mode_lens[i]);
                }
            raw.swap(filt);
            mode_lens.swap(filt_lens);
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
            mode_lens.push_back(static_cast<int>(rr));
        }
        if (raw.empty()) return out;

        // ---- Tukey rejection: length (1.5*IQR) --------------------------
        {
            std::vector<double> lens_d(mode_lens.begin(), mode_lens.end());
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

        // ---- Tukey rejection: 50%-upslope position (3.0*IQR) ------------
        // up50 is the horizontal alignment fiducial and is searched per beat
        // (unlike ECG's analytic R column), so we guard its position here.
        {
            std::vector<double> poss;
            poss.reserve(raw.size());
            for (const auto& r : raw) poss.push_back(static_cast<double>(r.up50));
            apply_mask(keep_within_tukey(poss, 3.0));
        }
        out.total_beats = raw.size();
        if (raw.empty()) return out;

        // ---- mode length -----------------------------------------------
        std::unordered_map<int, size_t> hist;
        for (int L : mode_lens) hist[L]++;
        for (const auto& kv : hist) {
            if (kv.second > out.mode_count
                || (kv.second == out.mode_count && kv.first > out.mode_length)) {
                out.mode_length = kv.first;
                out.mode_count = kv.second;
            }
        }
        for (size_t i = 0; i < raw.size(); ++i) {
            if (mode_lens[i] == out.mode_length) {
                out.mode_beat_index = static_cast<int>(i);
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
        // to the mode beat's. Feet no longer share a column (we aligned on the
        // upslope), so the window is centered per beat at foot_cols[i]. A
        // windowed mean, rather than the single argmin sample, avoids an
        // order-statistic bias at the foot. A constant vertical shift moves no
        // column, so it cannot affect the (cusp-free) horizontal alignment.
        if (out.mode_beat_index >= 0
            && out.mode_beat_index < static_cast<int>(out.beats.size()))
        {
            const int fb_w = std::max(1, out.mode_length / 50);
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

            const double target = foot_baseline(static_cast<size_t>(out.mode_beat_index));
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