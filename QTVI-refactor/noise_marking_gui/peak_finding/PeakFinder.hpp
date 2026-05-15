// ============================================================================
// File: PeakFinder.hpp
// Find peaks in signal data with minimum distance constraint
// ============================================================================
#pragma once

#include <vector>
#include <algorithm>
#include <cstdlib>
#include <set>

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

    // 2. Priority sort (tallest peaks first -- essential for MATLAB matching)
    std::sort(candidates.begin(), candidates.end(), [](const PeakCandidate& a, const PeakCandidate& b) {
        if (a.val != b.val) return a.val > b.val;
        return a.pos < b.pos;
        });

    // 3. Elimination based on MinPeakDistance.
    //    Original: O(N^2) by marking every survivor that lies within range of
    //    each accepted peak. Optimized: O(N log N) by keeping accepted peak
    //    positions in a sorted set, and for each candidate checking only the
    //    two neighbors closest in position. Equivalent semantics: a candidate
    //    is rejected iff some earlier-accepted peak is at strict distance
    //    < minPeakDistance, exactly as before.
    const long long minDist = (long long)minPeakDistance;
    std::set<size_t> accepted_positions;
    std::vector<PeakCandidate> final_peaks;
    final_peaks.reserve(candidates.size());

    for (size_t i = 0; i < candidates.size(); ++i) {
        const size_t cpos = candidates[i].pos;
        bool ok = true;

        if (minDist > 0 && !accepted_positions.empty()) {
            auto it = accepted_positions.lower_bound(cpos);
            // Neighbor at or after cpos
            if (it != accepted_positions.end()) {
                long long d = (long long)*it - (long long)cpos;
                if (d < 0) d = -d;
                if (d < minDist) ok = false;
            }
            // Neighbor strictly before cpos
            if (ok && it != accepted_positions.begin()) {
                auto prev = std::prev(it);
                long long d = (long long)cpos - (long long)*prev;
                if (d < 0) d = -d;
                if (d < minDist) ok = false;
            }
        }

        if (ok) {
            final_peaks.push_back(candidates[i]);
            accepted_positions.insert(cpos);
        }
    }

    // 4. Sort back to temporal order
    std::sort(final_peaks.begin(), final_peaks.end(), [](const PeakCandidate& a, const PeakCandidate& b) {
        return a.pos < b.pos;
        });

    pks.reserve(final_peaks.size());
    locs.reserve(final_peaks.size());
    for (const auto& p : final_peaks) {
        pks.push_back(p.val);
        locs.push_back(p.pos);
    }
}