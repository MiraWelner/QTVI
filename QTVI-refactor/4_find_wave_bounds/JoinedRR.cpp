// ============================================================================
// File: JoinedRR.cpp
// ============================================================================
#include <numeric>
#include <cmath>
#include <algorithm>
#include "JoinedRR.h"
#include "rpeakdetect.h"
#include "pan_tompkin.h"
#include "ecgLms.h"
#include "RRsimpleSquared.h"
#include "StatsUtils.h"
#include "FilterUtils.h"

// Helper function to create sorted list
vector<vector<double>> sortedList(const vector<vector<size_t>>& output, const vector<double>& weight) {
    size_t total_length = 0;
    for (const auto& vec : output) {
        total_length += vec.size();
    }

    vector<vector<double>> lst(total_length, vector<double>(2));

    size_t prev = 0;
    for (size_t i = 0; i < output.size(); ++i) {
        for (size_t j = 0; j < output[i].size(); ++j) {
            lst[prev][0] = output[i][j];
            lst[prev][1] = weight[i];
            prev++;
        }
    }

    std::sort(lst.begin(), lst.end(),
        [](const vector<double>& a, const vector<double>& b) {
            return a[0] < b[0];
        });

    return lst;
}

vector<size_t> RPeakfromRWave(const vector<double>& ecg, const vector<size_t>& rWaveIdx) {
    if (rWaveIdx.size() <= 1) {
        return rWaveIdx;
    }

    // Calculate median difference to determine window size
    vector<double> diffs;
    for (size_t i = 1; i < rWaveIdx.size(); ++i) {
        diffs.push_back(static_cast<double>(rWaveIdx[i]) - rWaveIdx[i - 1]);
    }
    std::sort(diffs.begin(), diffs.end());
    double avg_diff = diffs[diffs.size() / 2];
    int half_window_size = static_cast<int>(std::round(avg_diff / 6.0));

    vector<size_t> ridxs = rWaveIdx;

    for (size_t i = 0; i < rWaveIdx.size(); ++i) {
        // Define the search window boundaries
        int start_idx = static_cast<int>(rWaveIdx[i]) - half_window_size;
        int end_idx = static_cast<int>(rWaveIdx[i]) + half_window_size;

        // Clip to vector bounds
        size_t winstart = std::max(0, start_idx);
        size_t winend = std::min(static_cast<int>(ecg.size()) - 1, end_idx);

        if (winstart >= winend) continue;

        // Find the local maximum within the window
        double max_val = ecg[winstart];
        size_t max_id = winstart;

        for (size_t j = winstart + 1; j <= winend; ++j) {
            if (ecg[j] > max_val) {
                max_val = ecg[j];
                max_id = j;
            }
        }
        ridxs[i] = max_id;
    }

    return ridxs;
}


vector<size_t> JoinedRR(const vector<double>& ecgSeg, double ecgSamplingRate, double diff_range) {
    if (std_dev(ecgSeg) == 0) {
        return vector<size_t>();
    }

    // Run multiple algorithms
    vector<vector<size_t>> output(6);
    vector<double> weights = { 0.75, 0.25, 0.25, 1.25, 1.5, 0.75 };

    try {
        auto result1 = rpeakdetect(ecgSeg, ecgSamplingRate);
        output[0] = result1.R_index;
    }
    catch (...) {
        output[0].clear();
    }

    try {
        auto result2 = rpeakdetect(ecgSeg, ecgSamplingRate, 0.1);
        output[1] = result2.R_index;
    }
    catch (...) {
        output[1].clear();
    }

    try {
        auto result3 = rpeakdetect(ecgSeg, ecgSamplingRate, 0.4);
        output[2] = result3.R_index;
    }
    catch (...) {
        output[2].clear();
    }

    try {
        auto result4 = pan_tompkin(ecgSeg, ecgSamplingRate);
        output[3] = result4.qrs_i_raw;
    }
    catch (...) {
        output[3].clear();
    }

    try {
        // Create butter filter coefficients
        vector<double> b, a;
        vector<double> Wn = { 5.0 * 2.0 / ecgSamplingRate, 12.0 * 2.0 / ecgSamplingRate };
        butter(3, Wn, b, a);

        vector<double> ecg_centered = ecgSeg;
        double meanVal = mean(ecgSeg);
        for (auto& val : ecg_centered) {
            val -= meanVal;
        }

        output[4] = ecgLms(ecg_centered, static_cast<int>(ecgSamplingRate), b, a, 0);

    }
    catch (...) {
        output[4].clear();
    }

    // Get initial potential peaks
    auto potentialPeaks = sortedList(output, weights);

    // Calculate median distance
    vector<double> median_dists;
    for (const auto& vec : output) {
        if (vec.size() > 1) {
            vector<double> diffs = diff(vector<double>(vec.begin(), vec.end()));
            double med = median(diffs);
            if (!std::isnan(med)) {
                median_dists.push_back(med);
            }
        }
    }

    double median_dist = median(median_dists);
    if (std::isnan(median_dist)) {
        median_dist = ecgSamplingRate * 0.6;  // default ~60 BPM
    }

    // Add simple squared method
    auto simpleResult = RRsimpleSquared(ecgSeg, median_dist / 2.0);
    output[5] = simpleResult.first;

    // Refine peaks for some algorithms
    for (size_t r = 3; r < 6; ++r) {
        output[r] = RPeakfromRWave(ecgSeg, output[r]);
    }

    // Get updated potential peaks
    potentialPeaks = sortedList(output, weights);

    // Shift peaks within diff_range to largest value
        // 1. Get unique peak locations
    std::map<size_t, double> weighted_map;
    for (const auto& p : potentialPeaks) {
        weighted_map[static_cast<size_t>(p[0])] += p[1];
    }

    // 2. Iterative Merging (Matches MATLAB diff_range logic)
    // We group peaks that are within diff_range and move them to the strongest local ECG point
    bool changed = true;
    while (changed) {
        changed = false;
        vector<size_t> keys;
        for (auto const& [key, val] : weighted_map) keys.push_back(key);
        std::sort(keys.begin(), keys.end());

        for (size_t i = 0; i + 1 < keys.size(); ++i) {
            if (keys[i + 1] - keys[i] <= diff_range) {
                // Determine which peak is "better" (higher voltage)
                size_t best_idx = (ecgSeg[keys[i]] > ecgSeg[keys[i + 1]]) ? keys[i] : keys[i + 1];
                size_t worse_idx = (best_idx == keys[i]) ? keys[i + 1] : keys[i];

                // Merge weights to the best index
                weighted_map[best_idx] += weighted_map[worse_idx];
                weighted_map.erase(worse_idx);

                changed = true;
                break; // Restart loop to ensure consistency
            }
        }
    }

    // 3. Apply the 2.4 Threshold
    vector<size_t> rr;
    for (auto const& [idx, weight] : weighted_map) {
        if (weight >= 2.4) {
            rr.push_back(idx);
        }
    }
    std::sort(rr.begin(), rr.end());
    return rr;
}