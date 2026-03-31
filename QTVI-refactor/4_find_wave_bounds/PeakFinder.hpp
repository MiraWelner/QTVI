// ============================================================================
// File: PeakFinder.hpp
// Find peaks in signal data with minimum distance constraint
// ============================================================================
#pragma once

#include <vector>
#include <algorithm>
#include <cstdlib>

namespace peakfinder_detail {

    struct PeakCandidate {
        double val;
        size_t pos;
    };

} // namespace peakfinder_detail

inline void findpeaks(const std::vector<double>& data,
    std::vector<double>& pks,
    std::vector<size_t>& locs,
    double minPeakDistance = 0)
{
    using peakfinder_detail::PeakCandidate;

    pks.clear();
    locs.clear();
    if (data.size() < 3) return;

    std::vector<PeakCandidate> candidates;

    // 1. Identify all local maxima (and center of plateaus)
    for (size_t i = 1; i < data.size() - 1; ++i) {
        if (data[i] > data[i - 1]) {
            size_t j = i;
            while (j < data.size() - 1 && data[j] == data[j + 1]) j++;

            if (j < data.size() - 1 && data[j] > data[j + 1]) {
                candidates.push_back({ data[i], i + (j - i) / 2 });
                i = j;
            }
        }
    }

    // 2. Priority sort (tallest peaks first — essential for MATLAB matching)
    std::sort(candidates.begin(), candidates.end(), [](const PeakCandidate& a, const PeakCandidate& b) {
        if (a.val != b.val) return a.val > b.val;
        return a.pos < b.pos;
        });

    // 3. Elimination based on MinPeakDistance
    std::vector<bool> keep(candidates.size(), true);
    std::vector<PeakCandidate> final_peaks;

    for (size_t i = 0; i < candidates.size(); ++i) {
        if (!keep[i]) continue;
        final_peaks.push_back(candidates[i]);

        for (size_t j = i + 1; j < candidates.size(); ++j) {
            if (keep[j]) {
                long long dist = std::abs((long long)candidates[i].pos - (long long)candidates[j].pos);
                if (dist < (long long)minPeakDistance) {
                    keep[j] = false;
                }
            }
        }
    }

    // 4. Sort back to temporal order
    std::sort(final_peaks.begin(), final_peaks.end(), [](const PeakCandidate& a, const PeakCandidate& b) {
        return a.pos < b.pos;
        });

    for (const auto& p : final_peaks) {
        pks.push_back(p.val);
        locs.push_back(p.pos);
    }
}