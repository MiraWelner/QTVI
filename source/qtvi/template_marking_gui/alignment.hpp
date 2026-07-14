#pragma once
//
// alignment.hpp
//
// Per-bin beat alignment for the ECG. Runs BEFORE any normalization or
// template averaging: takes the raw ECG signal + the ch1 R-peak train
// for one bin and slices individual beats, then reports the "mode beat"
// -- the beat length that occurs most often (exact integer equality,
// no tolerance).
//
// Step 1 (this file, so far):
//   * Slice every beat as [R_i - 0.25 * RR, R_i + 0.75 * RR], where
//     RR = R_{i+1} - R_i. Length = 1.0 * RR by construction. R sits at
//     column 0.25 * RR inside the snip.
//   * The last R-peak has no following R, so it has no RR and no beat.
//   * Count how many beats share each integer length. The most common
//     length is the "mode length"; its count is printed to stderr as
//     a debug statement.
//
// Later steps (not yet implemented) will use the mode as the reference
// for vertical + horizontal alignment across beats.
//
// Usage:
//   auto res = alignment::extract_beats_and_mode(ecgSignal, rPeaks);
//   // res.beats[i] is one beat snip; res.mode_length / res.mode_count
//   // give the mode.
//

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <unordered_map>
#include <vector>

namespace alignment {

    struct BeatSet {
        // Each beat: length 1.0*RR, R at column 0.25*RR.
        // Beats near the edges of the signal that would need samples outside
        // [0, signal.size()) are dropped so every beat is fully populated.
        std::vector<std::vector<double>> beats;

        // The R-peak sample-index (in the ORIGINAL signal frame) that each
        // beat is centered on. Same length as `beats`.
        std::vector<size_t> r_indices;

        // Length statistics.
        int    mode_length = -1;   // most common beat length, samples (-1 if empty)
        size_t mode_count = 0;     // how many beats share that length
        size_t total_beats = 0;    // beats.size(), for convenience

        // Shape-cluster statistics (among beats of `mode_length`):
        //   Beats of the mode length are clustered by pairwise RMS: a beat
        //   joins an existing cluster iff RMS(beat - cluster.first) <= tol.
        //   First-member membership test (order-dependent but transitive-ish
        //   under tight tolerance).
        //   `mode_group_size` = size of the LARGEST such cluster.
        //   `mode_group_indices` = indices (into `beats`) of the beats in
        //   that largest cluster.
        //   `mode_group_rms_tol` = the tolerance used (default 0.001).
        size_t mode_group_size = 0;
        std::vector<size_t> mode_group_indices;
        double mode_group_rms_tol = 0.0;

        // After pass 1 (R-align) and pass 2 (Q-align), `beats` is REPLACED
        // by NaN-padded aligned versions. Every beat has the same width;
        // real samples occupy a subrange, NaN elsewhere.
        //   r_aligned_col  = column where R sits on the shared axis
        //   q_aligned_col  = column where Q sits on the shared axis (-1 if
        //                    no Q detected on mode beat)
        //   q_cols[i]      = column where Q sits for beat i (pass 1 axis
        //                    coordinates BEFORE pass 2 shifting; -1 if not
        //                    detected)
        //   mode_beat_index = index of the mode beat used for Q anchoring
        //                     (first beat of the winning shape cluster; -1
        //                     if no cluster)
        int r_aligned_col = -1;
        int q_aligned_col = -1;
        std::vector<int> q_cols;
        int mode_beat_index = -1;
    };

    // Slice every beat using the [-0.25 RR, +0.75 RR] rule. `rPeaks` must be
    // sorted, in the same sample frame as `signal`. Beats whose full extent
    // would run past the ends of the signal are skipped.
    //
    // Also clusters the mode-length beats by shape (RMS <= rms_tol, first-
    // member test). Reports the largest cluster's size / indices.
    inline BeatSet extract_beats_and_mode(const std::vector<double>& signal,
        const std::vector<size_t>& rPeaks,
        double rms_tol = 0.05)
    {
        BeatSet out;
        const int64_t N = static_cast<int64_t>(signal.size());
        if (N == 0 || rPeaks.size() < 2) {
            fprintf(stderr, "[align] no beats: signal=%lld rPeaks=%zu\n",
                static_cast<long long>(N), rPeaks.size());
            return out;
        }

        // First pass: gather RRs so we can reject outliers before slicing.
        // A single missed R-peak makes one RR huge (~2x median or worse); a
        // spurious R inserted between real ones makes RRs tiny. Both wreck
        // the shared axis width. Median +/- 2x is a permissive but robust
        // gate.
        std::vector<int64_t> rrs;
        rrs.reserve(rPeaks.size());
        for (size_t i = 0; i + 1 < rPeaks.size(); ++i) {
            const int64_t rr = static_cast<int64_t>(rPeaks[i + 1])
                - static_cast<int64_t>(rPeaks[i]);
            if (rr > 3) rrs.push_back(rr);
        }
        double rr_median = 0.0;
        if (!rrs.empty()) {
            std::vector<int64_t> tmp = rrs;
            std::nth_element(tmp.begin(), tmp.begin() + tmp.size() / 2, tmp.end());
            rr_median = static_cast<double>(tmp[tmp.size() / 2]);
        }
        const int64_t rr_min = static_cast<int64_t>(rr_median * 0.5);
        const int64_t rr_max = static_cast<int64_t>(rr_median * 2.0);

        size_t rejected_rr = 0;
        for (size_t i = 0; i + 1 < rPeaks.size(); ++i) {
            const int64_t r0 = static_cast<int64_t>(rPeaks[i]);
            const int64_t r1 = static_cast<int64_t>(rPeaks[i + 1]);
            const int64_t rr = r1 - r0;
            if (rr <= 3) continue;                             // degenerate
            if (rr_median > 0 && (rr < rr_min || rr > rr_max)) {
                ++rejected_rr; continue;                       // RR outlier
            }

            const int64_t before = rr / 4;                    // 0.25 RR
            const int64_t after = rr - before;                // 0.75 RR
            const int64_t len = before + after;               // 1.0 RR (== rr, integer safe)
            const int64_t start = r0 - before;
            const int64_t end = r0 + after;                   // exclusive

            if (start < 0 || end > N) continue;               // partial beat -> drop

            std::vector<double> beat(static_cast<size_t>(len));
            std::memcpy(beat.data(), signal.data() + start,
                static_cast<size_t>(len) * sizeof(double));

            out.beats.push_back(std::move(beat));
            out.r_indices.push_back(rPeaks[i]);
        }
        out.total_beats = out.beats.size();
        if (rejected_rr > 0) {
            fprintf(stderr, "[align] rejected %zu RR-outlier beats "
                "(median=%.0f, keep range=[%lld, %lld])\n",
                rejected_rr, rr_median, static_cast<long long>(rr_min),
                static_cast<long long>(rr_max));
        }

        // Mode length: integer equality, no tolerance.
        std::unordered_map<int, size_t> hist;
        hist.reserve(out.beats.size() * 2);
        for (const auto& b : out.beats)
            hist[static_cast<int>(b.size())]++;

        for (const auto& kv : hist) {
            if (kv.second > out.mode_count
                || (kv.second == out.mode_count && kv.first > out.mode_length)) {
                out.mode_length = kv.first;
                out.mode_count = kv.second;
            }
        }

        // Step 2: cluster the mode-length beats by shape (RMS <= rms_tol,
        // first-member membership test). The winning "mode ECG waveform"
        // group is the largest cluster.
        out.mode_group_rms_tol = rms_tol;
        if (out.mode_length > 0 && out.mode_count > 0) {
            // Gather the indices of beats with the mode length.
            std::vector<size_t> modeIdx;
            modeIdx.reserve(out.mode_count);
            for (size_t i = 0; i < out.beats.size(); ++i)
                if (static_cast<int>(out.beats[i].size()) == out.mode_length)
                    modeIdx.push_back(i);

            // First-member cluster assignment. Each cluster stores the
            // index of its "seed" (first) beat plus every subsequent beat
            // whose RMS-diff to the seed is <= tol.
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
                        break;   // first-member test: stop at first match
                    }
                }
                if (!placed) {
                    Cluster c;
                    c.seed = idx;
                    c.members.push_back(idx);
                    clusters.push_back(std::move(c));
                }
            }

            // Largest cluster wins; ties go to the first-formed cluster.
            for (const auto& c : clusters) {
                if (c.members.size() > out.mode_group_size) {
                    out.mode_group_size = c.members.size();
                    out.mode_group_indices = c.members;
                }
            }
        }

        // ================================================================
        // Pass 1: R alignment (in place, NaN-pad prepend, replace beats).
        // Anchor R at column A = 0.25 * max(RR) on a shared axis of width
        // max(RR). Each beat's R sits at column 0.25 * RR_i within its own
        // vector; prepending (A - 0.25*RR_i) NaN samples slides it to A.
        // Then append NaN so every beat has the same total width.
        // ================================================================
        int max_len = 0;
        for (const auto& b : out.beats)
            if (static_cast<int>(b.size()) > max_len) max_len = static_cast<int>(b.size());
        if (max_len > 0) {
            const int R_anchor = max_len / 4;                 // 0.25 * max_RR
            const int shared_w = max_len;                     // R at A, data ends at A + 0.75*RR_i <= max_RR

            std::vector<std::vector<double>> aligned;
            aligned.reserve(out.beats.size());
            const double NaND = std::numeric_limits<double>::quiet_NaN();
            for (const auto& b : out.beats) {
                const int L = static_cast<int>(b.size());
                const int r_in_beat = L / 4;                  // 0.25 * L
                const int prepend = R_anchor - r_in_beat;     // shift right by prepend
                std::vector<double> a(shared_w, NaND);
                // Place b[0..L) into a[prepend..prepend+L). Both bounds
                // stay inside [0, shared_w) by construction.
                for (int k = 0; k < L; ++k) a[prepend + k] = b[k];
                aligned.push_back(std::move(a));
            }
            out.beats = std::move(aligned);
            out.r_aligned_col = R_anchor;

            // Remap mode_group_indices? The indices refer to the same beat
            // slots (only vector contents changed, not order), so no remap.

            // Pass 2: Q detection + Q alignment.
            // Q lives in [0, R_anchor) of the mode beat's now-aligned data.
            // Detection: last local min (first-derivative sign change - -> +)
            // before R, in the beat's REAL (non-NaN) pre-R region. Fallback:
            // steepest upstroke (argmax of first derivative) in the pre-R
            // region.
            auto detectQ = [&](const std::vector<double>& beat, int R_col) -> int {
                // Find the beat's real pre-R region [firstReal, R_col).
                int firstReal = 0;
                while (firstReal < R_col && std::isnan(beat[firstReal])) ++firstReal;
                if (R_col - firstReal < 3) return -1;   // too short to search

                // Look for the LAST local min in [firstReal+1, R_col-1).
                int last_min = -1;
                for (int i = firstReal + 1; i < R_col - 1; ++i) {
                    if (std::isnan(beat[i - 1]) || std::isnan(beat[i]) || std::isnan(beat[i + 1])) continue;
                    if (beat[i] <= beat[i - 1] && beat[i] <= beat[i + 1]) last_min = i;
                }
                if (last_min >= 0) return last_min;

                // Fallback: argmax of first derivative in [firstReal+1, R_col).
                int best = -1;
                double best_slope = -std::numeric_limits<double>::infinity();
                for (int i = firstReal + 1; i < R_col; ++i) {
                    if (std::isnan(beat[i - 1]) || std::isnan(beat[i])) continue;
                    const double slope = beat[i] - beat[i - 1];
                    if (slope > best_slope) { best_slope = slope; best = i; }
                }
                return best;
                };

            // Q column for every beat + the mode beat.
            out.q_cols.assign(out.beats.size(), -1);
            for (size_t i = 0; i < out.beats.size(); ++i)
                out.q_cols[i] = detectQ(out.beats[i], R_anchor);

            int Q_mode = -1;
            if (!out.mode_group_indices.empty()) {
                const size_t mIdx = out.mode_group_indices.front();
                Q_mode = out.q_cols[mIdx];
            }
            out.mode_beat_index = out.mode_group_indices.empty()
                ? -1 : static_cast<int>(out.mode_group_indices.front());
            out.q_aligned_col = Q_mode;

            // Q-align: shift each beat so its Q lands at Q_mode. Beats with
            // no Q detected (q_cols[i] < 0), or whose required shift exceeds
            // the sanity cap, are left as-is (R-aligned, no Q shift). Cap
            // is 1/8 of max_len (~roughly 100 ms at typical HR). Q-detection
            // failures on noisy beats used to blow up the shared width by
            // dragging one beat's data far out; the cap absorbs those.
            if (Q_mode >= 0) {
                const int q_shift_cap = std::max(3, max_len / 8);
                size_t rejected_q = 0;
                // Compute the new shared width from max prepend + longest data.
                int max_prepend = 0;
                int max_end = 0;
                std::vector<int> shifts(out.beats.size(), 0);
                for (size_t i = 0; i < out.beats.size(); ++i) {
                    if (out.q_cols[i] < 0) { shifts[i] = 0; continue; }
                    const int delta = Q_mode - out.q_cols[i];   // shift right by delta
                    if (std::abs(delta) > q_shift_cap) {
                        shifts[i] = 0;    // sanity-reject: no Q shift for this beat
                        out.q_cols[i] = -1;
                        ++rejected_q;
                        continue;
                    }
                    shifts[i] = delta;
                    // Find last real sample of this beat to know its extent.
                    int lastReal = static_cast<int>(out.beats[i].size()) - 1;
                    while (lastReal >= 0 && std::isnan(out.beats[i][lastReal])) --lastReal;
                    if (lastReal < 0) continue;
                    int firstReal = 0;
                    while (firstReal < static_cast<int>(out.beats[i].size())
                        && std::isnan(out.beats[i][firstReal])) ++firstReal;
                    const int newStart = firstReal + delta;
                    const int newEnd = lastReal + delta;
                    if (-newStart > max_prepend) max_prepend = -newStart;
                    if (newEnd > max_end) max_end = newEnd;
                }
                if (rejected_q > 0)
                    fprintf(stderr, "[align] rejected %zu Q-shifts > cap=%d\n",
                        rejected_q, q_shift_cap);
                if (max_prepend < 0) max_prepend = 0;
                const int new_width = max_prepend + max_end + 1;

                std::vector<std::vector<double>> q_aligned;
                q_aligned.reserve(out.beats.size());
                for (size_t i = 0; i < out.beats.size(); ++i) {
                    std::vector<double> a(new_width, NaND);
                    const int shift = shifts[i] + max_prepend;   // total left-pad on new axis
                    for (int k = 0; k < static_cast<int>(out.beats[i].size()); ++k) {
                        const double v = out.beats[i][k];
                        if (std::isnan(v)) continue;
                        const int dst = k + shift;
                        if (dst >= 0 && dst < new_width) a[dst] = v;
                    }
                    q_aligned.push_back(std::move(a));
                }
                out.beats = std::move(q_aligned);
                out.q_aligned_col = Q_mode + max_prepend;
                out.r_aligned_col = R_anchor + max_prepend;
            }
        }

        fprintf(stderr,
            "[align] beats=%zu unique_lengths=%zu mode_length=%d mode_count=%zu "
            "mode_group_size=%zu rms_tol=%.4f R_col=%d Q_col=%d\n",
            out.total_beats, hist.size(), out.mode_length, out.mode_count,
            out.mode_group_size, out.mode_group_rms_tol,
            out.r_aligned_col, out.q_aligned_col);

        return out;
    }

}   // namespace alignment