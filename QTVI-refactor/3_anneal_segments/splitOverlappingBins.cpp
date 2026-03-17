#include "splitOverlappingBins.hpp"

// MATLAB: splitOverLappingBins.m
// Split any exclusion that spans multiple bins so each piece lies
// within a single bin. Sort by idx_start afterward.
//
// IMPORTANT: MATLAB iterates over the original size only; new entries
// are appended at the end and not revisited.

void splitOverlappingBins(
    std::vector<Exclusion>& exclusions,
    const std::vector<uint64_t>& bin_breaks)
{
    size_t original_size = exclusions.size();
    for (size_t i = 0; i < original_size; ++i) {
        if (exclusions[i].bin_start != exclusions[i].bin_end) {
            uint64_t temp_end = exclusions[i].idx_end;
            int temp_bin_end = exclusions[i].bin_end;

            for (int bin = exclusions[i].bin_start; bin <= temp_bin_end; ++bin) {
                if (bin == exclusions[i].bin_start) {
                    exclusions[i].idx_end = bin_breaks[bin - 1];
                    exclusions[i].bin_end = bin;
                }
                else if (bin == temp_bin_end) {
                    exclusions.push_back({
                        bin_breaks[bin - 2], temp_end, bin, bin
                        });
                }
                else {
                    exclusions.push_back({
                        bin_breaks[bin - 2], bin_breaks[bin - 1], bin, bin
                        });
                }
            }
        }
    }

    std::sort(exclusions.begin(), exclusions.end(),
        [](const Exclusion& a, const Exclusion& b) {
            return a.idx_start < b.idx_start;
        });
}
