// ============================================================================
// File: diff2.cpp
// ============================================================================
#include "diff2.h"

vector<double> diff2(const vector<double>& X, int nd) {
    if (X.empty()) {
        return vector<double>();
    }

    vector<double> result = X;

    for (int k = 0; k < nd; ++k) {
        size_t n = result.size();

        if (n <= 1) {
            return vector<double>();
        }

        vector<double> X1(n - 1);
        for (size_t i = 0; i < n - 1; ++i) {
            X1[i] = result[i + 1] - result[i];  // slope of the index point and the point after
        }

        vector<double> X2(n - 1, 0.0);
        for (size_t i = 1; i < n - 1; ++i) {
            X2[i] = X1[i - 1];  // slope of the index point and the point before
        }

        vector<double> X3(n - 1, 0.0);
        for (size_t i = 2; i < n - 1; ++i) {
            X3[i] = X1[i - 2];  // slope of the two points before the index point
        }

        vector<double> X4(n - 1, 0.0);
        for (size_t i = 0; i < n - 2; ++i) {
            X4[i] = X1[i + 1];  // slope of the two points after the index point
        }

        result.resize(n - 1);
        for (size_t i = 0; i < n - 1; ++i) {
            result[i] = (X1[i] * 2 + X2[i] * 2 + X3[i] + X4[i]) / 6.0;
        }
    }

    return result;
}
