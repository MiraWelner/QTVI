

// ============================================================================
// File: JoinedRR.cpp
// ============================================================================
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

// Helper to get R peak from R wave (simplified)
vector<size_t> RPeakfromRWave(const vector<double>& ecgSeg, const vector<size_t>& rwaves) {
    // This is a placeholder - actual implementation would refine peak locations
    return rwaves;
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
    vector<double> uniq;
    for (const auto& peak : potentialPeaks) {
        if (uniq.empty() || uniq.back() != peak[0]) {
            uniq.push_back(peak[0]);
        }
    }

    vector<bool> mask(uniq.size() - 1, false);
    for (size_t i = 0; i < uniq.size() - 1; ++i) {
        if (uniq[i + 1] - uniq[i] <= diff_range) {
            mask[i] = true;
        }
    }

    for (size_t i = 0; i < mask.size(); ++i) {
        if (mask[i]) {
            size_t idx1 = static_cast<size_t>(uniq[i]);
            size_t idx2 = static_cast<size_t>(uniq[i + 1]);

            if (idx1 < ecgSeg.size() && idx2 < ecgSeg.size()) {
                double target = ecgSeg[idx1] > ecgSeg[idx2] ? uniq[i] : uniq[i + 1];

                for (auto& peak : potentialPeaks) {
                    if (peak[0] == uniq[i] || peak[0] == uniq[i + 1]) {
                        peak[0] = target;
                    }
                }
            }
        }
    }

    // Recalculate unique
    uniq.clear();
    for (const auto& peak : potentialPeaks) {
        if (uniq.empty() || uniq.back() != peak[0]) {
            uniq.push_back(peak[0]);
        }
    }

    // Calculate weighted peaks
    map<double, double> weighted_peaks_map;
    for (const auto& peak : potentialPeaks) {
        weighted_peaks_map[peak[0]] += peak[1];
    }

    vector<vector<double>> weighted_peaks;
    for (const auto& entry : weighted_peaks_map) {
        weighted_peaks.push_back({ entry.first, entry.second });
    }

    // Filter by threshold
    vector<size_t> rr;
    for (const auto& peak : weighted_peaks) {
        if (peak[1] >= 2.4) {
            rr.push_back(static_cast<size_t>(peak[0]));
        }
    }

    return rr;
}