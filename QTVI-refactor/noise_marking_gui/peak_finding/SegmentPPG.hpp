/**
 * @file   SegmentPPG.hpp
 * @brief  Segment PPG signal into individual pulses (minima and maxima).
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-05-17
 */
#pragma once

#include "nanfastsmooth.hpp"
#include "StatsUtils.hpp"
#include "stdoutlier.hpp"
#include "RunLength.hpp"
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <vector>

struct SegmentPPGResult {
    vector<size_t> minAmps;  // valley (trough) indices
    vector<size_t> maxAmps;  // peak indices
};

namespace segmentppg_detail {

    // Segment signal into peaks and valleys based on above/below baseline mask
    inline pair<vector<size_t>, vector<size_t>> segBeats(const vector<double>& ppg, const vector<int>& mask)
    {
        vector<int> B;
        vector<double> N, BI;
        RunLength(mask, B, N, BI);

        vector<size_t> peaks, valleys;
        peaks.reserve(B.size() / 2);
        valleys.reserve(B.size() / 2);

        for (size_t i = 0; i < B.size(); ++i) {
            size_t start = static_cast<size_t>(BI[i]);
            size_t end = std::min(start + static_cast<size_t>(N[i]) - 1, ppg.size() - 1);

            if (B[i] == 0)
                valleys.push_back(start + min_element_index(ppg, start, end + 1).second);
            else
                peaks.push_back(start + max_element_index(ppg, start, end + 1).second);
        }

        return { peaks, valleys };
    }

    // Recompute a valley location from surrounding peaks when flagged as outlier
    inline size_t newVallyFromPeaks(
        const vector<double>& ppg, size_t currVal, size_t currPeak,
        const vector<size_t>& peaks_idxs, const vector<size_t>& vally_idxs,
        const vector<bool>& peakoutliers, const vector<bool>& vallyoutlier)
    {
        if (currPeak == 0 || currPeak >= peaks_idxs.size())
            return vally_idxs[currVal];

        if (!peakoutliers[currPeak - 1] && !peakoutliers[currPeak]) {
            size_t start = peaks_idxs[currPeak - 1];
            size_t end = peaks_idxs[currPeak];
            return start + min_element_index(ppg, start, end + 1).second;
        }

        return vally_idxs[currVal];
    }

} // namespace segmentppg_detail

inline SegmentPPGResult SegmentPPG(const vector<double>& ppg, double ppgRate) {
    if (ppg.empty()) throw std::runtime_error("Empty PPG signal");

    // 1. Smooth and compute baseline
    vector<double> ppg_smooth = nanfastsmooth(ppg, ppgRate * 0.25, 3);
    vector<double> baseline = movmean(ppg_smooth, static_cast<size_t>(ppgRate));

    // 2. Build mask: 1 where smoothed signal is above baseline, 0 below
    vector<int> mask(ppg.size());
    for (size_t i = 0; i < ppg.size(); ++i)
        mask[i] = (ppg_smooth[i] > baseline[i]) ? 1 : 0;

    // 3. Segment into peaks and valleys
    auto [peakidx, vallyidx] = segmentppg_detail::segBeats(ppg, mask);

    if (peakidx.empty() || vallyidx.empty())
        throw std::runtime_error("Incomplete beat segmentation (flat signal)");
    if (std::abs(static_cast<int>(peakidx.size()) - static_cast<int>(vallyidx.size())) > 1)
        throw std::runtime_error("Peak/valley count mismatch exceeds tolerance");

    // 4. Detect outliers in timing
    vector<bool> peak_outliers = stdoutlier(
        vector<double>(peakidx.begin(), peakidx.end()), 2.5, 100, "both", false);
    vector<bool> valley_outliers = stdoutlier(
        vector<double>(vallyidx.begin(), vallyidx.end()), 2.5, 100, "both", false);

    // 5. Correct outlier valleys using neighboring peaks
    size_t pidx = (vallyidx[0] < peakidx[0]) ? 0 : 1;

    for (size_t vidx = 0; vidx < vallyidx.size() - 1; ++vidx) {
        if (vidx == 0 || vidx == vallyidx.size() - 2) {
            pidx++;
            continue;
        }

        if (valley_outliers[vidx]) {
            vallyidx[vidx] = segmentppg_detail::newVallyFromPeaks(
                ppg, vidx, pidx, peakidx, vallyidx, peak_outliers, valley_outliers);
        }
        pidx++;
    }

    return { vallyidx, peakidx };
}