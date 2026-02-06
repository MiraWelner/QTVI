
// ============================================================================
// File: SegmentPPG.cpp
// ============================================================================
#include "SegmentPPG.h"
#include "nanfastsmooth.h"
#include "StatsUtils.h"
#include "stdoutlier.h"
#include "RunLength.h"

// Helper function to segment beats
pair<vector<size_t>, vector<size_t>> segBeats(const vector<double>& ppg, const vector<int>& mask);

// Helper function to find new valley from peaks
size_t newVallyFromPeaks(const vector<double>& ppg, size_t currVal, size_t currPeak,
    const vector<size_t>& peaks_idxs, const vector<size_t>& vally_idxs,
    const vector<bool>& peakoutliers, const vector<bool>& vallyoutlier);

SegmentPPGResult SegmentPPG(const vector<double>& ppg, double sampleRate) {
    SegmentPPGResult result;

    // Aggressive smooth and tight moving mean
    vector<double> ppg_smooth = nanfastsmooth(ppg, sampleRate * 0.25, 3);
    vector<double> M = movmean(ppg_smooth, static_cast<size_t>(sampleRate));

    vector<double> plow = ppg;
    vector<int> p_mask(ppg.size(), 0);

    for (size_t i = 0; i < ppg.size(); ++i) {
        if (ppg_smooth[i] > M[i]) {
            plow[i] = NaN;
            p_mask[i] = 1;
        }
        else {
            p_mask[i] = 0;
        }
    }

    auto beatResult = segBeats(ppg, p_mask);
    vector<size_t> peakidx = beatResult.first;
    vector<size_t> vallyidx = beatResult.second;

    if (std::abs(static_cast<int>(peakidx.size()) - static_cast<int>(vallyidx.size())) > 1) {
        throw std::runtime_error("not expecting this");
    }

    vector<bool> ppg_outliers_peaks = stdoutlier(
        vector<double>(peakidx.begin(), peakidx.end()), 2.5, 100, "both", false);
    vector<bool> vallyoutlier_time = stdoutlier(
        vector<double>(vallyidx.begin(), vallyidx.end()), 2.5, 100, "both", false);

    size_t pidx;
    if (!vallyidx.empty() && !peakidx.empty()) {
        if (vallyidx[0] < peakidx[0]) {
            pidx = 0;  // valley peak valley peak ....
        }
        else {
            pidx = 1;  // peak valley peak valley ....
        }
    }
    else {
        pidx = 0;
    }

    for (size_t vidx = 0; vidx < vallyidx.size() - 1; ++vidx) {
        if (vidx == 0 || vidx == vallyidx.size() - 2) {
            pidx++;
            continue;
        }

        if (vallyoutlier_time[vidx]) {
            vallyidx[vidx] = newVallyFromPeaks(ppg, vidx, pidx, peakidx, vallyidx,
                ppg_outliers_peaks, vallyoutlier_time);
        }

        pidx++;
    }

    result.ppgMinAmps = vallyidx;
    result.maxAmps = peakidx;

    return result;
}

// Helper function implementation
pair<vector<size_t>, vector<size_t>> segBeats(const vector<double>& ppg, const vector<int>& mask) {
    // Run-length encoding
    vector<int> B;
    vector<double> N;
    vector<double> BI;
    RunLength(mask, B, N, BI);

    // Separate high and low regions
    vector<pair<size_t, size_t>> high, low;

    for (size_t i = 0; i < B.size(); ++i) {
        size_t start = static_cast<size_t>(BI[i]);
        size_t end = start + static_cast<size_t>(N[i]) - 1;

        if (B[i] == 0) {
            low.push_back({ start, end });
        }
        else {
            high.push_back({ start, end });
        }
    }

    // Find maxima in high regions
    vector<size_t> peaks;
    for (const auto& region : high) {
        auto maxResult = max_element_index(ppg, region.first, region.second + 1);
        peaks.push_back(region.first + maxResult.second - 1);
    }

    // Find minima in low regions
    vector<size_t> valleys;
    for (const auto& region : low) {
        auto minResult = min_element_index(ppg, region.first, region.second + 1);
        valleys.push_back(region.first + minResult.second - 1);
    }

    return { peaks, valleys };
}

size_t newVallyFromPeaks(const vector<double>& ppg, size_t currVal, size_t currPeak,
    const vector<size_t>& peaks_idxs, const vector<size_t>& vally_idxs,
    const vector<bool>& peakoutliers, const vector<bool>& vallyoutlier) {
    if (currPeak == 0 || currPeak >= peaks_idxs.size()) {
        return vally_idxs[currVal];
    }

    size_t start = peaks_idxs[currPeak - 1];
    size_t end = peaks_idxs[currPeak];

    if (start >= ppg.size() || end >= ppg.size() || start >= end) {
        return vally_idxs[currVal];
    }

    auto minResult = min_element_index(ppg, start, end + 1);
    size_t idx = start + minResult.second;

    // Check if peaks are not outliers
    if (currPeak > 0 && currPeak < peakoutliers.size()) {
        if (!peakoutliers[currPeak - 1] && !peakoutliers[currPeak]) {
            if (idx == vally_idxs[currVal]) {
                return vally_idxs[currVal];
            }
            else {
                return idx;
            }
        }
    }

    return vally_idxs[currVal];
}