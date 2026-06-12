// ============================================================================
// File: pairRtoPPGBeat.hpp
// Pair ECG R-peaks to PPG pulse valleys
// Returns matrix where each row is [ppg_valley_idx, ecg_R_idx]
// -1.0 indicates unpaired
// ============================================================================
#pragma once

#include "RunLength.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>

namespace pairrtoppg_detail {
    static constexpr int MAX_CONFLICT_ITERATIONS = 100;
} // namespace pairrtoppg_detail

inline vector<vector<double>> pairRtoPPGBeat(
    const vector<double>& ecg,
    const vector<double>& ppg,
    const vector<size_t>& ecgRIndex,
    const vector<size_t>& ppgMinAmps,
    double ecgRate,
    double ppgRate)
{
    if (ecg.empty() || ppg.empty()) return {};

    // 1. Create time vectors (in minutes)
    vector<double> ecgtime(ecg.size());
    vector<double> ppgtime(ppg.size());

    for (size_t i = 0; i < ecg.size(); ++i)
        ecgtime[i] = (double)i / ecgRate / 60.0;
    for (size_t i = 0; i < ppg.size(); ++i)
        ppgtime[i] = (double)i / ppgRate / 60.0;

    // 2. Initialize pairs: [ppg_valley_idx, ecg_R_idx]
    vector<vector<double>> pairs;
    pairs.reserve(ecgRIndex.size());

    for (size_t i = 0; i < ecgRIndex.size(); ++i) {
        if (ecgRIndex[i] < ecg.size())
            pairs.push_back({ std::numeric_limits<double>::quiet_NaN(), (double)ecgRIndex[i] });
    }

    // 3. Pair every R-peak to the nearest PPG valley within the neighboring R-R interval.
    //    Performance: ppgMinAmps is sorted in temporal order, and ppgtime is monotone
    //    in index, so a per-R-peak linear scan over ppgMinAmps was O(R*P). We replace
    //    it with two binary searches per R-peak (O(log P)). The original bounds check
    //    (p_idx < ppgtime.size()) is preserved by trimming ppgMinAmps to a valid prefix.
    const size_t ppg_valid_end = static_cast<size_t>(
        std::lower_bound(ppgMinAmps.begin(), ppgMinAmps.end(), ppgtime.size())
        - ppgMinAmps.begin());

    for (size_t i = 0; i < pairs.size(); ++i) {
        size_t beginidx, endidx;

        if (i == 0) {
            beginidx = 0;
            endidx = (pairs.size() > 1) ? static_cast<size_t>(pairs[i + 1][1]) : ecg.size() - 1;
        }
        else if (i == pairs.size() - 1) {
            beginidx = static_cast<size_t>(pairs[i - 1][1]);
            endidx = ecg.size() - 1;
        }
        else {
            beginidx = static_cast<size_t>(pairs[i - 1][1]);
            endidx = static_cast<size_t>(pairs[i + 1][1]);
        }

        beginidx = std::min(beginidx, ecg.size() - 1);
        endidx = std::min(endidx, ecg.size() - 1);

        double begtime = ecgtime[beginidx];
        double endtime = ecgtime[endidx];

        // Range of ppgMinAmps[j] entries whose time lies in [begtime, endtime].
        // Comparators read ppgtime at the (in-bounds) PPG index.
        auto lo_it = std::lower_bound(
            ppgMinAmps.begin(), ppgMinAmps.begin() + ppg_valid_end, begtime,
            [&](size_t pidx, double t) { return ppgtime[pidx] < t; });
        auto hi_it = std::upper_bound(
            lo_it, ppgMinAmps.begin() + ppg_valid_end, endtime,
            [&](double t, size_t pidx) { return t < ppgtime[pidx]; });

        if (lo_it != hi_it) {
            double minError = Inf;
            size_t best_ppg_idx = 0;
            size_t ecg_r_idx = static_cast<size_t>(pairs[i][1]);

            if (ecg_r_idx < ecgtime.size()) {
                const double r_t = ecgtime[ecg_r_idx];
                for (auto it = lo_it; it != hi_it; ++it) {
                    size_t p_idx = *it;
                    double error = std::abs(ppgtime[p_idx] - r_t);
                    if (error < minError) {
                        minError = error;
                        best_ppg_idx = p_idx;
                    }
                }
            }
            pairs[i][0] = (double)best_ppg_idx;
        }
        else {
            pairs[i][0] = -1.0;
        }
    }

    // 4. Add unpaired PPG valleys.
    //    Original: O(P*R) by linear-scanning pairs for each candidate.
    //    Optimized: O(R+P) using a hash set of already-paired PPG indices.
    //    The original compared raw doubles for equality; pairs[k][0] at this
    //    point is either a value that came directly from (double)ppgMinAmps[i]
    //    (i.e., an exactly representable size_t cast) or -1.0 / NaN. So
    //    casting back to size_t for the "paired" set is safe and exact.
    //    We update the set as we push new rows so that any (pathological)
    //    duplicate in ppgMinAmps is still de-duplicated as the original did.
    {
        std::unordered_set<size_t> paired_ppg;
        paired_ppg.reserve(pairs.size() * 2 + ppgMinAmps.size());
        for (const auto& row : pairs) {
            double v = row[0];
            if (!std::isnan(v) && v >= 0.0) {
                paired_ppg.insert(static_cast<size_t>(v));
            }
        }
        for (size_t i = 0; i < ppgMinAmps.size(); ++i) {
            auto ins = paired_ppg.insert(ppgMinAmps[i]);
            if (ins.second) {
                pairs.push_back({ (double)ppgMinAmps[i], -1.0 });
            }
        }
    }

    // 5. Sort by PPG index
    std::sort(pairs.begin(), pairs.end(), [](const vector<double>& a, const vector<double>& b) {
        if (std::isnan(a[0]) || std::isnan(b[0])) return false;
        if (a[0] != b[0]) return a[0] < b[0];
        return a[1] < b[1];
        });

    // 6. Conflict resolution: when multiple R-peaks map to the same PPG valley,
    //    keep the closest one and unlink the rest.
    //
    //    Original: after fixing one run, break and restart -- O(K*P) for K runs.
    //    Optimized: process every run in one RLE pass. Runs in the snapshot are
    //    disjoint (different PPG values or separated by other rows), and we only
    //    mutate column 0 of *loser* rows within each run, never touching column
    //    0 of any other run. So fixing all runs in one pass yields the same
    //    final state as the original. The outer while is kept (with its safety
    //    bound) but in practice exits after one iteration.
    bool changed = true;
    int iteration = 0;
    while (changed && iteration++ < pairrtoppg_detail::MAX_CONFLICT_ITERATIONS) {
        changed = false;

        vector<double> first_col;
        first_col.reserve(pairs.size());
        for (const auto& row : pairs) first_col.push_back(row[0]);

        vector<double> B, N, BI;
        RunLength(first_col, B, N, BI);

        for (size_t i = 0; i < N.size(); ++i) {
            if (N[i] >= 2 && B[i] >= 0.0) {
                size_t start_idx = static_cast<size_t>(BI[i]);
                size_t end_idx = start_idx + static_cast<size_t>(N[i]) - 1;

                if (end_idx >= pairs.size()) continue;

                double min_err = Inf;
                size_t winner_row = start_idx;
                size_t ppg_idx = static_cast<size_t>(pairs[start_idx][0]);

                for (size_t j = start_idx; j <= end_idx; ++j) {
                    double r_idx = pairs[j][1];
                    if (r_idx >= 0 && r_idx < (double)ecg.size()) {
                        double err = std::abs(ppgtime[ppg_idx] - ecgtime[(size_t)r_idx]);
                        if (err < min_err) {
                            min_err = err;
                            winner_row = j;
                        }
                    }
                }

                for (size_t j = start_idx; j <= end_idx; ++j) {
                    if (j == winner_row) continue;
                    pairs[j][0] = -1.0;
                    changed = true;
                }
                // No break: process every run in this pass. See note above.
            }
        }
    }

    // 7. Final bounds cleanup
    for (auto& row : pairs) {
        if (row[0] >= (double)ppg.size()) row[0] = (double)ppg.size() - 1;
        if (row[1] >= (double)ecg.size()) row[1] = (double)ecg.size() - 1;
        if (std::isnan(row[0]) || row[0] < 0) row[0] = -1.0;
        if (std::isnan(row[1]) || row[1] < 0) row[1] = -1.0;
    }

    // 8. Orphaned ECG peaks with no PPG match: use ECG index as fallback
    for (size_t i = 0; i < pairs.size(); ++i) {
        if (pairs[i][0] == -1.0 && pairs[i][1] >= 0.0)
            pairs[i][0] = pairs[i][1];
    }

    return pairs;
}