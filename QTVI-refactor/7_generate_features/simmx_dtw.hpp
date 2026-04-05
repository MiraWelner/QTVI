#pragma once

#include "common.hpp"
#include <vector>
#include <cmath>

/**
 * @file simmx_dtw.hpp
 * @brief Compute a similarity (cost) matrix between two PLA-decomposed signals
 *        using the absolute difference of piecewise slopes.
 *
 * Port of simmx_dtw.m.
 */

namespace ppg {

    /**
     * @brief Result of the DTW similarity matrix computation.
     */
    struct SimmxResult {
        Matrix w;                ///< Cost matrix (|slope1[i] - slope2[j]|).
        std::vector<double> ta;  ///< Time-axis weights for signal 1 (segment lengths).
        std::vector<double> tb;  ///< Time-axis weights for signal 2 (segment lengths).
    };

    /**
     * @brief Build the slope-difference cost matrix for DTW alignment.
     * @param y1    Signal 1 values.
     * @param pla1  PLA breakpoint indices for signal 1 (0-based).
     * @param y2    Signal 2 values.
     * @param pla2  PLA breakpoint indices for signal 2 (0-based).
     * @return SimmxResult containing cost matrix w and time-axis vectors ta, tb.
     */
    inline SimmxResult simmx_dtw(
        const std::vector<double>& y1, const std::vector<int>& pla1,
        const std::vector<double>& y2, const std::vector<int>& pla2)
    {
        const int n1 = static_cast<int>(pla1.size());
        const int n2 = static_cast<int>(pla2.size());

        // Compute slopes and time-axis weights for signal 1
        std::vector<double> slope1(n1, 0.0);
        std::vector<double> ta(n1, 0.0);
        ta[0] = 1.0;
        for (int i = 1; i < n1; ++i) {
            int denom = pla1[i] - pla1[i - 1];
            slope1[i] = (denom != 0)
                ? (y1[pla1[i]] - y1[pla1[i - 1]]) / static_cast<double>(denom)
                : 0.0;
            ta[i] = static_cast<double>(denom);
        }

        // Compute slopes and time-axis weights for signal 2
        std::vector<double> slope2(n2, 0.0);
        std::vector<double> tb(n2, 0.0);
        tb[0] = 1.0;
        for (int i = 1; i < n2; ++i) {
            int denom = pla2[i] - pla2[i - 1];
            slope2[i] = (denom != 0)
                ? (y2[pla2[i]] - y2[pla2[i - 1]]) / static_cast<double>(denom)
                : 0.0;
            tb[i] = static_cast<double>(denom);
        }

        // Build cost matrix: w(i,j) = |slope1[i] - slope2[j]|
        Matrix w(n1, n2);
        for (int i = 0; i < n1; ++i)
            for (int j = 0; j < n2; ++j)
                w(i, j) = std::abs(slope1[i] - slope2[j]);

        return { std::move(w), std::move(ta), std::move(tb) };
    }

} // namespace ppg