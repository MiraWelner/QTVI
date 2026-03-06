#include "SegmentPPG.h"
#include "nanfastsmooth.h"
#include "StatsUtils.h"
#include "stdoutlier.h"
#include "RunLength.h"
#include <stdexcept>
#include <cmath>
#include <algorithm>

// Internal helper for segmentation logic
pair<vector<size_t>, vector<size_t>> segBeats(const vector<double>& ppg, const vector<int>& mask);
size_t newVallyFromPeaks(const vector<double>& ppg, size_t currVal, size_t currPeak,
    const vector<size_t>& peaks_idxs, const vector<size_t>& vally_idxs,
    const vector<bool>& peakoutliers, const vector<bool>& vallyoutlier);

SegmentPPGResult SegmentPPG(const vector<double>& ppg, double sampleRate) {
    SegmentPPGResult result;

    if (ppg.empty()) throw std::runtime_error("Empty signal");

    // Match MATLAB: ppg_smooth = nanfastsmooth(ppg, sampleRate * .25, 3);
    vector<double> ppg_smooth = nanfastsmooth(ppg, sampleRate * 0.25, 3);
    // Match MATLAB: M = movmean(ppg_smooth, sampleRate); 
    // Ensure movmean uses a centered window to match MATLAB's default
    vector<double> M = movmean(ppg_smooth, static_cast<size_t>(sampleRate));

    vector<int> p_mask(ppg.size(), 0);
    for (size_t i = 0; i < ppg.size(); ++i) {
        // Match MATLAB: p_mask(ppg_smooth > M) = 1
        p_mask[i] = (ppg_smooth[i] > M[i]) ? 1 : 0;
    }

    auto beatResult = segBeats(ppg, p_mask);
    vector<size_t> peakidx = beatResult.first;
    vector<size_t> vallyidx = beatResult.second;

    // CRITICAL: Replicate MATLAB crash/error behavior for flat/bad signals
    // MATLAB: if vallyidx(1) < peakidx(1) will throw if either is empty
    if (peakidx.empty() || vallyidx.empty()) {
        throw std::runtime_error("Incomplete beat segmentation (flat signal)");
    }

    if (std::abs(static_cast<int>(peakidx.size()) - static_cast<int>(vallyidx.size())) > 1) {
        throw std::runtime_error("not expecting this");
    }

    vector<bool> ppg_outliers_peaks = stdoutlier(
        vector<double>(peakidx.begin(), peakidx.end()), 2.5, 100, "both", false);
    vector<bool> vallyoutlier_time = stdoutlier(
        vector<double>(vallyidx.begin(), vallyidx.end()), 2.5, 100, "both", false);

    size_t pidx;
    // MATLAB is 1-based, C++ 0-based adjustment
    if (vallyidx[0] < peakidx[0]) {
        pidx = 0;
    }
    else {
        pidx = 1;
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

pair<vector<size_t>, vector<size_t>> segBeats(const vector<double>& ppg, const vector<int>& mask) {
    vector<int> B;
    vector<double> N, BI;
    RunLength(mask, B, N, BI);

    vector<size_t> peaks, valleys;
    for (size_t i = 0; i < B.size(); ++i) {
        size_t start = static_cast<size_t>(BI[i]);
        size_t end = start + static_cast<size_t>(N[i]) - 1;
        if (end >= ppg.size()) end = ppg.size() - 1;

        if (B[i] == 0) { // Low region (Valleys)
            auto res = min_element_index(ppg, start, end + 1);
            valleys.push_back(start + res.second);
        }
        else { // High region (Peaks)
            auto res = max_element_index(ppg, start, end + 1);
            peaks.push_back(start + res.second);
        }
    }
    return { peaks, valleys };
}

size_t newVallyFromPeaks(const vector<double>& ppg, size_t currVal, size_t currPeak,
    const vector<size_t>& peaks_idxs, const vector<size_t>& vally_idxs,
    const vector<bool>& peakoutliers, const vector<bool>& vallyoutlier) {

    if (currPeak == 0 || currPeak >= peaks_idxs.size()) return vally_idxs[currVal];

    size_t start = peaks_idxs[currPeak - 1];
    size_t end = peaks_idxs[currPeak];
    auto res = min_element_index(ppg, start, end + 1);
    size_t new_idx = start + res.second;

    if (!peakoutliers[currPeak - 1] && !peakoutliers[currPeak]) {
        return new_idx;
    }
    return vally_idxs[currVal];
}