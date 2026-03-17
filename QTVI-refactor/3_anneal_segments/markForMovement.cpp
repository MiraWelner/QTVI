#include "markForMovement.hpp"
#include "getExclusionIntervals.hpp"

// MATLAB: markForMovement.m
// For each bin that contains exclusions, determine the good intervals
// and assign movement direction (left/right/none) and a move flag.

std::vector<Section> markForMovement(
    const std::vector<Exclusion>& exclusions,
    const std::vector<uint64_t>& bin_breaks,
    int bin_count,
    double ppgSR,
    uint64_t bin_size_samples,
    double min_bin_size_mins)
{
    std::set<int> update_bins_set;
    for (const auto& ex : exclusions) update_bins_set.insert(ex.bin_start);

    std::vector<Section> marked;

    for (int cur_bin : update_bins_set) {
        uint64_t bin_begin = bin_breaks[cur_bin - 1] - bin_size_samples;
        uint64_t bin_end = bin_breaks[cur_bin - 1];
        double bin_half = (double)(bin_end - (bin_size_samples / 2));

        // Get exclusions for this bin
        std::vector<std::pair<uint64_t, uint64_t>> bin_excl;
        for (const auto& ex : exclusions) {
            if (ex.bin_start == cur_bin)
                bin_excl.push_back({ ex.idx_start, ex.idx_end });
        }

        // MATLAB: good = getExclusionIntervals(bin_begin, bin_end, bin_exclusions)
        auto good = getExclusionIntervals(bin_begin, bin_end, bin_excl);

        // For each good interval, compute movement direction
        // MATLAB: [m, movement_dir] = max((good - bin_half)');
        struct PotentialSection {
            uint64_t begin, end;
            double m;
            int dir;
            int flag;
        };
        std::vector<PotentialSection> potential;

        for (size_t gi = 0; gi < good.size(); ++gi) {
            double val_begin = (double)good[gi].first - bin_half;
            double val_end = (double)good[gi].second - bin_half;
            double m = std::max(val_begin, val_end);

            double good_time_mins = ((double)(good[gi].second - good[gi].first) / ppgSR) / 60.0;
            bool too_small = good_time_mins < min_bin_size_mins;

            int movement_dir;
            if (cur_bin == 1) {
                movement_dir = (m <= 0) ? 2 : 0;
            }
            else if (cur_bin == bin_count) {
                movement_dir = (m >= 0) ? 1 : 0;
            }
            else {
                movement_dir = (m <= 0) ? 1 : 0;
            }

            int move_flag;
            if (!too_small) {
                movement_dir = 0;
                move_flag = 0;
            }
            else {
                move_flag = 1;
            }

            potential.push_back({ good[gi].first, good[gi].second, m, movement_dir, move_flag });
        }

        // MATLAB: filter out small edge segments for first/last bin
        for (const auto& p : potential) {
            if (cur_bin == 1 && p.m <= 0 && p.flag == 1) continue;
            if (cur_bin == bin_count && p.m >= 0 && p.flag == 1) continue;
            marked.push_back({ p.begin, p.end, p.dir, p.flag });
        }
    }

    return marked;
}