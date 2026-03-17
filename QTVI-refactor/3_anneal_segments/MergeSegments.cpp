#include "MergeSegments.hpp"

// MATLAB: MergeSegments (called inline in AnnealSegments.m)
// Merges overlapping or touching [start, end] pairs.

std::vector<std::pair<uint64_t, uint64_t>> MergeSegments(
    std::vector<std::pair<uint64_t, uint64_t>> segs)
{
    if (segs.empty()) return {};
    std::sort(segs.begin(), segs.end());
    std::vector<std::pair<uint64_t, uint64_t>> merged;
    merged.push_back(segs[0]);
    for (size_t i = 1; i < segs.size(); ++i) {
        if (segs[i].first <= merged.back().second)
            merged.back().second = std::max(merged.back().second, segs[i].second);
        else
            merged.push_back(segs[i]);
    }
    return merged;
}