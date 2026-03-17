#include "getExclusionIntervals.hpp"

// MATLAB: getExclusionIntervals.m
// Given interval (a, b) and a set of break pairs, return the good intervals.
// MATLAB: exclusions = [a, reshape(breaks', 1, []), b]
//         exclusions = reshape(exclusions', 2, [])' -> Nx2 pairs
//
// Example: a=0, breaks=[[10,20],[30,40]], b=100
//   flat = [0, 10, 20, 30, 40, 100]
//   pairs = [[0,10], [20,30], [40,100]]

std::vector<std::pair<uint64_t, uint64_t>> getExclusionIntervals(
    uint64_t a,
    uint64_t b,
    const std::vector<std::pair<uint64_t, uint64_t>>& breaks)
{
    std::vector<uint64_t> flat;
    flat.push_back(a);
    for (const auto& br : breaks) {
        flat.push_back(br.first);
        flat.push_back(br.second);
    }
    flat.push_back(b);

    std::vector<std::pair<uint64_t, uint64_t>> good;
    for (size_t j = 0; j + 1 < flat.size(); j += 2) {
        uint64_t len = flat[j + 1] - flat[j];
        if (len > 0)
            good.push_back({ flat[j], flat[j + 1] });
    }
    return good;
}