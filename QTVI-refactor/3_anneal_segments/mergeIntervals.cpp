#include "mergeIntervals.hpp"

// MATLAB: mergeIntervals.m
// Merge adjacent good_sections that share a border and are both flagged
// for movement.
//
// IMPORTANT: MATLAB increments i unconditionally (i = i + 1 at end of while),
// even after a merge+delete. This means after merging, the next pair is
// skipped. We replicate that behavior exactly.

void mergeIntervals(
    std::vector<Section>& good_sections,
    double ppgSR,
    double min_bin_size_mins)
{
    size_t i = 0;
    while (i + 1 < good_sections.size()) {
        if (good_sections[i].end == good_sections[i + 1].begin &&
            good_sections[i].flag != 0 && good_sections[i + 1].flag != 0)
        {
            double seg1_min = ((double)(good_sections[i].end - good_sections[i].begin) / ppgSR) / 60.0;
            double seg2_min = ((double)(good_sections[i + 1].end - good_sections[i + 1].begin) / ppgSR) / 60.0;

            if (seg1_min + seg2_min >= min_bin_size_mins) {
                good_sections[i].dir = 0;
                good_sections[i].flag = 0;
            }
            else {
                // MATLAB: [~, idx] = max([seg1 seg2]); idx = idx - 1;
                int idx = (seg1_min >= seg2_min) ? 0 : 1;
                good_sections[i].dir = good_sections[i + idx].dir;
                good_sections[i].flag = 1;
            }

            good_sections[i].end = good_sections[i + 1].end;
            good_sections.erase(good_sections.begin() + i + 1);
        }

        // MATLAB: i = i + 1 always (even after merge/delete)
        ++i;
    }
}
