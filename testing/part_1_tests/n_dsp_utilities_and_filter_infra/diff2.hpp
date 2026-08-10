/**
 * @file   diff2.hpp
 * @brief  Diff2.hpp computes a smoothed numerical derivative. Instead of a simple first difference (x[i+1] - x[i]), 
           it averages four slopes around each point:

                 Forward slope: x[i+1] - x[i] (weight 2)
                 Backward slope: x[i] - x[i-1] (weight 2)
                 Two-back slope: x[i-1] - x[i-2] (weight 1)
                 Two-forward slope: x[i+2] - x[i+1] (weight 1)

           The result is their weighted average: (2 x forward + 2 x backward + 1 x twoBack + 1 x twoForward) / 6
           This is essentially a Savitzky-Golay-style smoothed derivative. By blending nearby slopes with heavier 
           weight on the immediate neighbors and lighter weight on the more distant ones, it estimates the local rate of 
           change while suppressing high-frequency noise that a simple difference would amplify. The nd parameter lets you 
           apply this repeatedly for higher-order derivative
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-30
 */
#pragma once
#include <vector>
using std::vector;
inline vector<double> diff2(const vector<double>& X, int nd = 1) {
    if (X.empty()) return {};

    vector<double> result = X;

    for (int k = 0; k < nd; ++k) {
        size_t n = result.size();
        if (n <= 1) return {};

        vector<double> slopeForward(n - 1);
        for (size_t i = 0; i < n - 1; ++i)
            slopeForward[i] = result[i + 1] - result[i];

        vector<double> slopeBack(n - 1, 0.0);
        for (size_t i = 1; i < n - 1; ++i)
            slopeBack[i] = slopeForward[i - 1];

        vector<double> slopeTwoBack(n - 1, 0.0);
        for (size_t i = 2; i < n - 1; ++i)
            slopeTwoBack[i] = slopeForward[i - 2];

        vector<double> slopeTwoForward(n - 1, 0.0);
        for (size_t i = 0; i < n - 2; ++i)
            slopeTwoForward[i] = slopeForward[i + 1];

        result.resize(n - 1);
        for (size_t i = 0; i < n - 1; ++i) {
            result[i] = (slopeForward[i] * 2 + slopeBack[i] * 2 +
                slopeTwoBack[i] + slopeTwoForward[i]) / 6.0;
        }
    }

    return result;
}