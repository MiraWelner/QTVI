// ============================================================================
// File: diff2.cpp
// ============================================================================
#include "diff2.h"

vector<double> diff2(const vector<double>& X, int nd) {
    if (X.empty()) return {};

    vector<double> result = X;

    for (int k = 0; k < nd; ++k) {
        size_t n = result.size();
        if (n <= 1) return {};

        // Forward slope: result[i+1] - result[i]
        vector<double> slopeForward(n - 1);
        for (size_t i = 0; i < n - 1; ++i) {
            slopeForward[i] = result[i + 1] - result[i];
        }

        // Backward slope: slope of index point and the point before
        vector<double> slopeBack(n - 1, 0.0);
        for (size_t i = 1; i < n - 1; ++i) {
            slopeBack[i] = slopeForward[i - 1];
        }

        // Slope two points back
        vector<double> slopeTwoBack(n - 1, 0.0);
        for (size_t i = 2; i < n - 1; ++i) {
            slopeTwoBack[i] = slopeForward[i - 2];
        }

        // Slope two points forward
        vector<double> slopeTwoForward(n - 1, 0.0);
        for (size_t i = 0; i < n - 2; ++i) {
            slopeTwoForward[i] = slopeForward[i + 1];
        }

        result.resize(n - 1);
        for (size_t i = 0; i < n - 1; ++i) {
            result[i] = (slopeForward[i] * 2 + slopeBack[i] * 2 +
                slopeTwoBack[i] + slopeTwoForward[i]) / 6.0;
        }
    }

    return result;
}