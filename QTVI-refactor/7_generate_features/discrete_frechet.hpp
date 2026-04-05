#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>
#include <limits>

/**
 * @file discrete_frechet.hpp
 * @brief Discrete Fréchet distance between two polygonal curves.
 *
 * Port of DiscreteFrechetDist.m (ZCD, 2011/2013).
 * Reference: Eiter & Mannila, "Computing Discrete Fréchet Distance",
 * Technical Report 94/64, Vienna University of Technology, 1994.
 *
 * Uses iterative bottom-up DP with two-row rolling buffer (O(n) memory)
 * when only the distance is needed.
 */

namespace ppg {

    /**
     * @brief Compute the discrete Fréchet distance between two 1-D curves.
     *
     * Uses O(min(|P|,|Q|)) memory via a two-row rolling buffer.
     * No coupling sequence is computed — only the scalar distance.
     *
     * @param P     First curve values.
     * @param Q     Second curve values.
     * @param dfcn  Optional distance function (defaults to absolute difference).
     * @return The discrete Fréchet distance (coupling measure).
     */
    inline double discrete_frechet_dist(
        const std::vector<double>& P,
        const std::vector<double>& Q,
        std::function<double(double, double)> dfcn = nullptr)
    {
        if (P.empty() || Q.empty()) return 0.0;

        if (dfcn == nullptr)
            dfcn = [](double a, double b) { return std::abs(a - b); };

        // Ensure P is the longer curve so we iterate over rows (P)
        // and the rolling buffer has length cols (Q) = the shorter dimension.
        const std::vector<double>* rows = &P;
        const std::vector<double>* cols = &Q;
        if (P.size() < Q.size()) std::swap(rows, cols);

        const int nr = static_cast<int>(rows->size());
        const int nc = static_cast<int>(cols->size());

        // Two-row rolling buffer
        constexpr double INF = std::numeric_limits<double>::infinity();
        std::vector<double> prev(nc, INF);
        std::vector<double> curr(nc);

        // Row 0
        prev[0] = dfcn((*rows)[0], (*cols)[0]);
        for (int j = 1; j < nc; ++j)
            prev[j] = std::max(prev[j - 1], dfcn((*rows)[0], (*cols)[j]));

        // Rows 1..nr-1
        for (int i = 1; i < nr; ++i) {
            curr[0] = std::max(prev[0], dfcn((*rows)[i], (*cols)[0]));
            for (int j = 1; j < nc; ++j) {
                curr[j] = std::max(
                    std::min({ prev[j], prev[j - 1], curr[j - 1] }),
                    dfcn((*rows)[i], (*cols)[j]));
            }
            std::swap(prev, curr);
        }

        return prev[nc - 1];
    }

} // namespace ppg