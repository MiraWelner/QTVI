// ============================================================================
// File: pairRtoPPGBeat.hpp
// Pair ECG R-peaks to PPG pulse valleys
// Returns matrix where each row is [ppg_valley_idx, ecg_R_idx]
// -1.0 indicates unpaired
// ============================================================================
#pragma once

#include "SignalProcessingTypes.h"
#include "RunLength.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace pairrtoppg_detail {
    static constexpr int MAX_CONFLICT_ITERATIONS = 100;
} // namespace pairrtoppg_detail

inline vector<vector<double>> pairRtoPPGBeat(
    const vector<double>& ecg,
    const vector<double>& ppg,
    const vector<size_t>& ecgRIndex,
    const vector<size_t>& ppgMinAmps)
{
    if (ecg.empty() || ppg.empty()) return {};

    // 1. Create time vectors (in minutes)
    vector<double> ecgtime(ecg.size());
    vector<double> ppgtime(ppg.size());

    for (size_t i = 0; i < ecg.size(); ++i)
        ecgtime[i] = (double)i / ECG_SAMPLE_RATE / 60.0;
    for (size_t i = 0; i < ppg.size(); ++i)
        ppgtime[i] = (double)i / PPG_SAMPLE_RATE / 60.0;

    // 2. Initialize pairs: [ppg_valley_idx, ecg_R_idx]
    vector<vector<double>> pairs;
    pairs.reserve(ecgRIndex.size());

    for (size_t i = 0; i < ecgRIndex.size(); ++i) {
        if (ecgRIndex[i] < ecg.size())
            pairs.push_back({ NaN, (double)ecgRIndex[i] });
    }

    // 3. Pair every R-peak to the nearest PPG valley within the neighboring R-R interval
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

        vector<size_t> candidates;
        for (size_t j = 0; j < ppgMinAmps.size(); ++j) {
            size_t p_idx = ppgMinAmps[j];
            if (p_idx < ppgtime.size()) {
                double p_time = ppgtime[p_idx];
                if (p_time >= begtime && p_time <= endtime)
                    candidates.push_back(j);
            }
        }

        if (!candidates.empty()) {
            double minError = Inf;
            size_t best_ppg_idx = 0;
            size_t ecg_r_idx = static_cast<size_t>(pairs[i][1]);

            for (size_t p_j : candidates) {
                size_t p_idx = ppgMinAmps[p_j];
                if (p_idx < ppgtime.size() && ecg_r_idx < ecgtime.size()) {
                    double error = std::abs(ppgtime[p_idx] - ecgtime[ecg_r_idx]);
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

    // 4. Add unpaired PPG valleys
    for (size_t i = 0; i < ppgMinAmps.size(); ++i) {
        bool found = false;
        double p_val = (double)ppgMinAmps[i];
        for (size_t k = 0; k < pairs.size(); ++k) {
            if (pairs[k][0] == p_val) { found = true; break; }
        }
        if (!found)
            pairs.push_back({ p_val, -1.0 });
    }

    // 5. Sort by PPG index
    std::sort(pairs.begin(), pairs.end(), [](const vector<double>& a, const vector<double>& b) {
        if (std::isnan(a[0]) || std::isnan(b[0])) return false;
        if (a[0] != b[0]) return a[0] < b[0];
        return a[1] < b[1];
        });

    // 6. Conflict resolution: when multiple R-peaks map to the same PPG valley,
    //    keep the closest one and unlink the rest.
    bool changed = true;
    int iteration = 0;
    while (changed && iteration++ < pairrtoppg_detail::MAX_CONFLICT_ITERATIONS) {
        changed = false;

        vector<double> first_col;
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

                if (changed) break;
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