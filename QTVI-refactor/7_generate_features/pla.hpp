#pragma once

#include <vector>
#include <cmath>
#include <algorithm>

/**
 * @file pla.hpp
 * @brief Piecewise Linear Approximation (PLA) for ECG/PPG segmentation.
 *
 * Port of PLA.m (Vullings, Automated ECG Segmentation with DTW).
**/

namespace ppg {

    /**
     * @brief Result of piecewise linear approximation.
     */
    struct PLAResult {
        std::vector<double> signal;  ///< The (unmodified) input signal.
        std::vector<int>    breaks;  ///< Indices of PLA breakpoints (0-based).
    };

    /**
     * @brief Compute piecewise linear approximation breakpoints.
     * @param input  Input signal.
     * @param step   Initial step size (default 10).
     * @param thresh Distance threshold for segment acceptance (default 10).
     * @return PLAResult containing the signal and breakpoint indices.
     */
    inline PLAResult pla(const std::vector<double>& input, int step = 10, double thresh = 10.0) {
        const int n = static_cast<int>(input.size());
        PLAResult result;
        result.signal = input;

        if (n == 0) return result;

        result.breaks.push_back(0); // first breakpoint (0-based)

        int s1 = step;
        int i = 0;

        while (i < n) {
            int i_plus_s = std::min(i + s1, n - 1);
            bool interrupt = false;

            while (!interrupt) {
                int j = i + 1;
                while (j <= i_plus_s) {
                    double distance = input[i_plus_s] - input[i];
                    double denom = static_cast<double>(i_plus_s - i);
                    double dcur = input[j] - input[i] - (distance * (j - i) / denom);

                    if (std::abs(dcur) > thresh) {
                        s1 = j - i;
                        i_plus_s = i + s1;
                        j = i + 1;
                        interrupt = true;
                        continue;
                    }
                    ++j;
                }

                if (interrupt) {
                    result.breaks.push_back(j - 2); // 0-based
                    i = j - 3;
                    s1 = step;
                    break;
                }

                // Distance within threshold — expand segment
                if (i_plus_s >= n - 1) {
                    break;
                }
                i_plus_s = std::min(i_plus_s + s1, n - 1);
            }

            ++i;
        }

        // Ensure last point is a breakpoint
        if (result.breaks.empty() || result.breaks.back() != n - 1) {
            result.breaks.push_back(n - 1);
        }

        return result;
    }

} // namespace ppg