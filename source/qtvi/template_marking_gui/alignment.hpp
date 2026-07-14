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

        // Per-beat 1.0*RR value — the "mode-relevant length". Beat vectors
        // are 1.25*RR long (they include the tail up to the next R), so
        // beat.size() alone can't drive the mode calc. Same length as
        // `beats`.
        std::vector<int> mode_lens;

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

        for (size_t i = 0; i + 1 < rPeaks.size(); ++i) {
            const int64_t r0 = static_cast<int64_t>(rPeaks[i]);
            const int64_t r1 = static_cast<int64_t>(rPeaks[i + 1]);
            const int64_t rr = r1 - r0;
            if (rr <= 3) continue;                             // degenerate only

            const int64_t before = rr / 4;                    // 0.25 RR
            const int64_t after = rr + rr / 10;                // 1.1 RR: extend past next R
            const int64_t len = before + after;               // 1.25 RR total
            const int64_t start = r0 - before;
            const int64_t end = r0 + after;                   // exclusive, at next R

            // Slice from signal; anything outside [0, N) is NaN.
            std::vector<double> beat(static_cast<size_t>(len),
                std::numeric_limits<double>::quiet_NaN());
            const int64_t copyStart = std::max<int64_t>(0, start);
            const int64_t copyEnd = std::min<int64_t>(N, end);
            for (int64_t k = copyStart; k < copyEnd; ++k) {
                beat[static_cast<size_t>(k - start)] = signal[static_cast<size_t>(k)];
            }

            out.beats.push_back(std::move(beat));
            out.r_indices.push_back(rPeaks[i]);
            out.mode_lens.push_back(static_cast<int>(rr));    // 1.0*RR for mode calc
        }
        out.total_beats = out.beats.size();

        // Mode length: integer equality on the 1.0*RR value (mode_lens),
        // NOT on beat.size() which is 1.25*RR and would let the tail
        // influence the histogram.
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

        // Step 2: cluster the mode-length beats by shape (RMS <= rms_tol,
        // first-member membership test). The winning "mode ECG waveform"
        // group is the largest cluster.
        out.mode_group_rms_tol = rms_tol;
        if (out.mode_length > 0 && out.mode_count > 0) {
            // Gather the indices of beats with the mode length (using
            // the 1.0*RR value in mode_lens, not beat.size()).
            std::vector<size_t> modeIdx;
            modeIdx.reserve(out.mode_count);
            for (size_t i = 0; i < out.beats.size(); ++i)
                if (out.mode_lens[i] == out.mode_length)
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
        // Anchor R at column A = 0.25 * max(mode_lens) on a shared axis of
        // width 1.25 * max(mode_lens). Each beat's R sits at column
        // 0.25 * RR_i within its own vector (which is 1.25*RR_i long);
        // prepending (A - 0.25*RR_i) NaN samples slides R to A. Beats
        // whose 1.0*RR is smaller than the max still fit inside the
        // shared axis with NaN tail.
        // ================================================================
        int max_mode_len = 0;
        for (int L : out.mode_lens)
            if (L > max_mode_len) max_mode_len = L;
        if (max_mode_len > 0) {
            const int R_anchor = max_mode_len / 4;              // 0.25 * max_RR
            const int shared_w = R_anchor + max_mode_len + max_mode_len / 10;  // 1.35 * max_RR

            std::vector<std::vector<double>> aligned;
            aligned.reserve(out.beats.size());
            const double NaND = std::numeric_limits<double>::quiet_NaN();
            for (size_t i = 0; i < out.beats.size(); ++i) {
                const auto& b = out.beats[i];
                const int L = static_cast<int>(b.size());       // 1.25 * RR_i
                const int r_in_beat = out.mode_lens[i] / 4;     // 0.25 * RR_i
                const int prepend = R_anchor - r_in_beat;       // shift right by prepend
                std::vector<double> a(shared_w, NaND);
                for (int k = 0; k < L; ++k) {
                    const int dst = prepend + k;
                    if (dst >= 0 && dst < shared_w) a[dst] = b[k];
                }
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
            // is 1/8 of max_mode_len (~roughly 100 ms at typical HR). Q-detection
            // failures on noisy beats used to blow up the shared width by
            // dragging one beat's data far out; the cap absorbs those.
            if (Q_mode >= 0) {
                const int q_shift_cap = std::max(3, max_mode_len / 8);
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

                // ============================================================
                // Pass 3: PR-baseline vertical alignment.
                //
                // Every Q-aligned beat has its Q at the same column
                // (out.q_aligned_col). Pick a PR-baseline window ending just
                // before Q on that shared axis: width = max_mode_len / 20
                // (~5% of RR), ending max_mode_len / 50 (~2% of RR) before Q.
                // Both scale with sample rate.
                //
                // Compute the mode beat's PR-baseline mean = target level.
                // For every other beat: compute its own PR baseline over the
                // SAME column range, then subtract (beat_baseline - target)
                // from every non-NaN sample. Beats without a Q shift (q_cols
                // < 0), or whose PR window has too few valid samples, are
                // left untouched (still R-aligned, just no vertical shift).
                // ============================================================
                if (out.mode_beat_index >= 0
                    && out.mode_beat_index < static_cast<int>(out.beats.size())) {
                    const int Qc = out.q_aligned_col;
                    const int pr_w = std::max(3, max_mode_len / 20);
                    const int pr_gap = std::max(1, max_mode_len / 50);
                    const int pr_lo = std::max(0, Qc - pr_gap - pr_w);
                    const int pr_hi = std::max(pr_lo, Qc - pr_gap);   // exclusive

                    auto pr_baseline = [&](const std::vector<double>& beat) {
                        double sum = 0.0; int n = 0;
                        for (int k = pr_lo; k < pr_hi && k < static_cast<int>(beat.size()); ++k) {
                            const double v = beat[k];
                            if (!std::isnan(v)) { sum += v; ++n; }
                        }
                        return (n >= 3) ? std::make_pair(sum / n, n)
                            : std::make_pair(std::numeric_limits<double>::quiet_NaN(), n);
                        };

                    const auto [target, target_n] =
                        pr_baseline(out.beats[out.mode_beat_index]);

                    if (std::isnan(target)) {
                        fprintf(stderr, "[align] pass3: mode beat has no PR baseline "
                            "(window=[%d,%d) valid=%d) -- skipping DC shift\n",
                            pr_lo, pr_hi, target_n);
                    }
                    else {
                        size_t shifted = 0, skipped = 0;
                        for (size_t i = 0; i < out.beats.size(); ++i) {
                            // Beats that weren't Q-shifted (Q undetected or shift
                            // capped) don't have their PR region on the shared
                            // Q-axis, so the window doesn't correspond to their
                            // real PR interval -- skip them.
                            if (out.q_cols[i] < 0) { ++skipped; continue; }
                            const auto [b_base, b_n] = pr_baseline(out.beats[i]);
                            if (std::isnan(b_base)) { ++skipped; continue; }
                            const double d = b_base - target;
                            if (d == 0.0) { ++shifted; continue; }
                            for (double& v : out.beats[i])
                                if (!std::isnan(v)) v -= d;
                            ++shifted;
                        }
                        fprintf(stderr, "[align] pass3: PR-baseline vertical align "
                            "window=[%d,%d) target=%.4g shifted=%zu skipped=%zu\n",
                            pr_lo, pr_hi, target, shifted, skipped);
                    }
                }
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