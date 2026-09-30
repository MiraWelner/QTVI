#pragma once
/**
 * @file   beat_substitute.hpp
 * @brief  when a beat is premature (ectopic) or voted ectopic due to being surrounded by ectopic beats, it cannot
 *         be left as a hole, thus it is substituted with a blend of the previous and subsequent beats. The pvcs are determined
 *         by pcv_filter.hpp.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace beat_substitute {
    inline constexpr double ewma_alpha = 0.125; //1/8

    inline std::vector<double> substituteBeat(const std::vector<double>& avgOld, const std::vector<double>& current) {
        //substitute premature beat (current) with the exponentially weighted moving average of the average old beat and the current beat
        std::vector<double> out(avgOld.size());
        for (size_t j = 0; j < out.size(); ++j)
            out[j] = (1 - ewma_alpha) * avgOld[j] + ewma_alpha * current[j];
        return out;   // preserves temporal continuity
    }


    inline std::vector<double> substituteBeatNaNSafe(const std::vector<double>& avgOld, const std::vector<double>& current)
    {
        //wraps above but NAN safe
        const double NaN = std::numeric_limits<double>::quiet_NaN();
        const size_t n = std::max(avgOld.size(), current.size());
        std::vector<double> a(n, NaN), c(n, NaN);
        for (size_t j = 0; j < n; ++j) {
            const double av = (j < avgOld.size()) ? avgOld[j] : NaN;
            const double cv = (j < current.size()) ? current[j] : NaN;
            const bool na = std::isnan(av), nc = std::isnan(cv);
            if (na && nc) continue;          // both NaN -> blend yields NaN
            a[j] = na ? cv : av;
            c[j] = nc ? av : cv;
        }
        return substituteBeat(a, c);
    }
}  // namespace beat_substitute