#include "pairRtoPPGBeat.h"
#include "RunLength.h"
#include "nanfastsmooth.h"
#include "PeakFinder.h"
#include <algorithm>
#include <cmath>
#include <limits>

// Safety definitions
#ifndef NaN
#define NaN std::numeric_limits<double>::quiet_NaN()
#endif
#ifndef Inf
#define Inf std::numeric_limits<double>::infinity()
#endif

vector<vector<double>> pairRtoPPGBeat(const vector<double>& ecg, const vector<double>& ppg,
    double ecgSamplingRate, double ppgSamplingRate,
    const vector<size_t>& ecgRIndex, const vector<size_t>& ppgMinAmps) {

    // 0. Safety check for empty signals
    if (ecg.empty() || ppg.empty()) return {};

    // 1. Create time vectors
    vector<double> ecgtime(ecg.size());
    vector<double> ppgtime(ppg.size());

    for (size_t i = 0; i < ecg.size(); ++i) {
        ecgtime[i] = (double)i / ecgSamplingRate / 60.0;
    }

    for (size_t i = 0; i < ppg.size(); ++i) {
        ppgtime[i] = (double)i / ppgSamplingRate / 60.0;
    }

    // 2. Initialize pairs matrix: [ppg_valley_idx, ecg_R_idx]
    vector<vector<double>> pairs;
    pairs.reserve(ecgRIndex.size());

    for (size_t i = 0; i < ecgRIndex.size(); ++i) {
        // Validation: Only accept indices that actually exist in the current signal
        if (ecgRIndex[i] < ecg.size()) {
            pairs.push_back({ (double)NaN, (double)ecgRIndex[i] });
        }
    }

    // 3. Pair every R to a PPG idx
    for (size_t i = 0; i < pairs.size(); ++i) {
        size_t beginidx = 0, endidx = 0;

        // CRITICAL BOUNDS CHECK: Prevent 0xc0000005 by clamping indices
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

        // Clamp to signal range
        beginidx = std::min(beginidx, ecg.size() - 1);
        endidx = std::min(endidx, ecg.size() - 1);

        double begtime = ecgtime[beginidx];
        double enbtime = ecgtime[endidx];

        // Find possible PPG pairings
        vector<size_t> possible_ppg_parings;
        for (size_t j = 0; j < ppgMinAmps.size(); ++j) {
            size_t p_idx = ppgMinAmps[j];
            if (p_idx < ppgtime.size()) {
                double p_time = ppgtime[p_idx];
                if (p_time >= begtime && p_time <= enbtime) {
                    possible_ppg_parings.push_back(j);
                }
            }
        }

        if (!possible_ppg_parings.empty()) {
            double minError = 1e18;
            size_t best_ppg_idx = 0;
            size_t ecg_r_idx = static_cast<size_t>(pairs[i][1]);

            for (size_t p_j : possible_ppg_parings) {
                size_t p_idx = ppgMinAmps[p_j];
                // Safety check before error calculation
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
            if (pairs[k][0] == p_val) {
                found = true;
                break;
            }
        }
        if (!found) {
            pairs.push_back({ p_val, -1.0 });
        }
    }

    // 5. Sort by PPG index
    std::sort(pairs.begin(), pairs.end(), [](const vector<double>& a, const vector<double>& b) {
        if (std::isnan(a[0]) || std::isnan(b[0])) return false;
        if (a[0] != b[0]) return a[0] < b[0];
        return a[1] < b[1];
        });

    // 6. Conflict Resolution Loop (Resolve multiple Rs assigned to one PPG)
    bool changed = true;
    int safety_counter = 0;
    while (changed && safety_counter++ < 100) {
        changed = false;

        // RE-CALCULATE summary every iteration to avoid using stale indices (BI)
        vector<double> first_col;
        for (const auto& row : pairs) first_col.push_back(row[0]);

        vector<double> B, N, BI;
        RunLength(first_col, B, N, BI);

        for (size_t i = 0; i < N.size(); ++i) {
            // If N[i] > 1, the same PPG index (B[i]) is assigned to multiple ECG peaks
            if (N[i] >= 2 && B[i] >= 0.0) {
                size_t start_idx = static_cast<size_t>(BI[i]);
                size_t end_idx = start_idx + static_cast<size_t>(N[i]) - 1;

                if (end_idx >= pairs.size()) continue;

                // Pick the ECG peak with the smallest time error to the PPG valley
                double min_err = 1e18;
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

                // Keep the winner, set others to their own ECG index as a fallback
                for (size_t j = start_idx; j <= end_idx; ++j) {
                    if (j == winner_row) continue;
                    pairs[j][0] = pairs[j][1]; // Use ECG index as temporary match
                    changed = true;
                }

                if (changed) break; // Re-run RunLength on updated data
            }
        }
    }

    // 7. Final Safety Cleanup
    for (auto& row : pairs) {
        // Final bounds check to ensure output values are within signal ranges
        if (row[0] >= (double)ppg.size()) row[0] = (double)ppg.size() - 1;
        if (row[1] >= (double)ecg.size()) row[1] = (double)ecg.size() - 1;

        if (std::isnan(row[0]) || row[0] < 0) row[0] = -1.0;
        if (std::isnan(row[1]) || row[1] < 0) row[1] = -1.0;
    }

    return pairs;
}
