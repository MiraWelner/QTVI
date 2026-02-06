// ============================================================================
// File: pairRtoPPGBeat.cpp
// ============================================================================
#include "pairRtoPPGBeat.h"
#include "RunLength.h"
#include "nanfastsmooth.h"
#include "PeakFinder.h"

vector<vector<double>> pairRtoPPGBeat(const vector<double>& ecg, const vector<double>& ppg,
    double ecgSamplingRate, double ppgSamplingRate,
    const vector<size_t>& ecgRIndex, const vector<size_t>& ppgMinAmps) {
    // Create time vectors
    vector<double> ecgtime(ecg.size());
    vector<double> ppgtime(ppg.size());

    for (size_t i = 0; i < ecg.size(); ++i) {
        ecgtime[i] = i / ecgSamplingRate / 60.0;
    }

    for (size_t i = 0; i < ppg.size(); ++i) {
        ppgtime[i] = i / ppgSamplingRate / 60.0;
    }

    // Initialize pairs matrix: [ppg_valley_idx, ecg_R_idx]
    vector<vector<double>> pairs(ecgRIndex.size(), vector<double>(2, NaN));

    for (size_t i = 0; i < ecgRIndex.size(); ++i) {
        pairs[i][1] = ecgRIndex[i];
    }

    // Pair every R to a PPG idx
    for (size_t i = 0; i < ecgRIndex.size(); ++i) {
        size_t beginidx, endidx;

        if (i == 0) {
            beginidx = 0;
            endidx = ecgRIndex.size() > 1 ? ecgRIndex[i + 1] : ecgtime.size() - 1;
        }
        else if (i == ecgRIndex.size() - 1) {
            beginidx = ecgRIndex[i - 1];
            endidx = ecgtime.size() - 1;
        }
        else {
            beginidx = ecgRIndex[i - 1];
            endidx = ecgRIndex[i + 1];
        }

        double begtime = ecgtime[beginidx];
        double enbtime = ecgtime[endidx];

        // Find possible PPG pairings
        vector<size_t> possible_ppg_parings;
        for (size_t j = 0; j < ppgMinAmps.size(); ++j) {
            if (ppgMinAmps[j] < ppgtime.size()) {
                double ppg_time = ppgtime[ppgMinAmps[j]];
                if (begtime <= ppg_time && enbtime >= ppg_time) {
                    possible_ppg_parings.push_back(j);
                }
            }
        }

        // Find minimum error
        if (possible_ppg_parings.empty()) {
            pairs[i][0] = -1;
        }
        else {
            double minError = Inf;
            size_t minIdx = 0;

            for (size_t j = 0; j < possible_ppg_parings.size(); ++j) {
                size_t ppg_idx = ppgMinAmps[possible_ppg_parings[j]];
                if (ppg_idx < ppgtime.size() && ecgRIndex[i] < ecgtime.size()) {
                    double error = std::abs(ppgtime[ppg_idx] - ecgtime[ecgRIndex[i]]);
                    if (error < minError) {
                        minError = error;
                        minIdx = j;
                    }
                }
            }

            pairs[i][0] = ppgMinAmps[possible_ppg_parings[minIdx]];
        }
    }

    // Add unpaired PPG valleys
    for (size_t i = 0; i < ppgMinAmps.size(); ++i) {
        bool found = false;
        for (const auto& pair : pairs) {
            if (pair[0] == ppgMinAmps[i]) {
                found = true;
                break;
            }
        }

        if (!found) {
            pairs.push_back({ static_cast<double>(ppgMinAmps[i]), -1 });
        }
    }

    // Sort by PPG index
    std::sort(pairs.begin(), pairs.end(),
        [](const vector<double>& a, const vector<double>& b) {
            return a[0] < b[0];
        });

    // Process duplicates using RunLength
    vector<double> first_column;
    for (const auto& pair : pairs) {
        first_column.push_back(pair[0]);
    }

    vector<double> B, N, BI;
    RunLength(first_column, B, N, BI);

    // Correct R's assigned to more than 1 PPG
    size_t i = 0;
    while (i < N.size()) {
        if (N[i] == 2 && B[i] != -1) {
            size_t idx1 = static_cast<size_t>(BI[i]);
            size_t idx2 = idx1 + static_cast<size_t>(N[i]) - 1;

            if (idx2 >= pairs.size()) {
                i++;
                continue;
            }

            // Handle duplicate pairing
            if (idx2 + 1 < pairs.size() && pairs[idx2 + 1][1] == -1) {
                // Next R is unpaired
                pairs[idx2 + 1][1] = pairs[idx2][1];
                pairs.erase(pairs.begin() + idx2);
                RunLength(first_column, B, N, BI);
            }
            else {
                // Choose the one with smaller error
                double error1 = Inf, error2 = Inf;

                if (pairs[idx1][0] < ppgtime.size() && pairs[idx1][1] >= 0 &&
                    pairs[idx1][1] < ecgtime.size()) {
                    error1 = std::abs(ppgtime[static_cast<size_t>(pairs[idx1][0])] -
                        ecgtime[static_cast<size_t>(pairs[idx1][1])]);
                }

                if (pairs[idx2][0] < ppgtime.size() && pairs[idx2][1] >= 0 &&
                    pairs[idx2][1] < ecgtime.size()) {
                    error2 = std::abs(ppgtime[static_cast<size_t>(pairs[idx2][0])] -
                        ecgtime[static_cast<size_t>(pairs[idx2][1])]);
                }

                if (error1 < error2) {
                    // Keep first, need to find new PPG for second
                    pairs[idx2][0] = pairs[idx2][1];  // Use ECG index as temporary
                }
                else {
                    // Keep second, need to find new PPG for first
                    pairs[idx1][0] = pairs[idx1][1];  // Use ECG index as temporary
                }

                RunLength(first_column, B, N, BI);
            }
        }
        else if (N[i] > 2 && B[i] != -1) {
            // Multiple Rs for one PPG
            size_t idx1 = static_cast<size_t>(BI[i]);
            size_t idx2 = idx1 + static_cast<size_t>(N[i]) - 1;

            // Check for exact match
            bool found = false;
            for (size_t j = idx1; j <= idx2 && j < pairs.size(); ++j) {
                if (pairs[j][0] == pairs[j][1]) {
                    found = true;
                    // Update others
                    for (size_t k = idx1; k <= idx2 && k < pairs.size(); ++k) {
                        if (k != j) {
                            pairs[k][0] = pairs[k][1];
                        }
                    }
                    break;
                }
            }

            if (!found) {
                throw std::runtime_error("implement me: multiple Rs for one PPG without exact match");
            }

            RunLength(first_column, B, N, BI);
        }

        i++;

        if (i > N.size()) {
            break;
        }
    }

    // Final clean: propagate R index to PPG when R is not paired
    for (auto& pair : pairs) {
        if (pair[0] == -1) {
            pair[0] = pair[1];
        }
    }

    // Sort again
    std::sort(pairs.begin(), pairs.end(),
        [](const vector<double>& a, const vector<double>& b) {
            return a[0] < b[0];
        });

    return pairs;
}