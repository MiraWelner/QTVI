#include "RoundToClosestBin.hpp"

// MATLAB: RoundToClosestBin (called in AnnealSegments.m)
// For each index, find the first bin_break >= index. Return 1-based bin number.
// If index exceeds all breaks, returns the last bin.

std::vector<int> RoundToClosestBin(
    const std::vector<uint64_t>& bin_breaks,
    const std::vector<uint64_t>& indices)
{
    std::vector<int> result(indices.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        int best_bin = static_cast<int>(bin_breaks.size()); // default: last bin
        for (size_t b = 0; b < bin_breaks.size(); ++b) {
            if (indices[i] <= bin_breaks[b]) {
                best_bin = static_cast<int>(b) + 1;
                break;
            }
        }
        result[i] = best_bin;
    }
    return result;
}
