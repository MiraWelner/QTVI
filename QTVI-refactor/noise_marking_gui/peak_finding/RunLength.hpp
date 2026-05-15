// ============================================================================
// File: RunLength.h
// Run-length encoding and decoding (header-only, templated)
// ============================================================================
#pragma once

#include "SignalProcessingTypes.hpp"
#include <type_traits>

// Encode: consecutive equal values are grouped into (value, count, start_index).
//
// Two paths via if constexpr:
//   - Floating-point T: original behavior, treating NaN==NaN as equal. Tracks
//     currentIndex so the NaN check reads the original value at that index.
//   - Non-floating T (e.g. int): plain equality with B.back(). Equivalent
//     because std::isnan on a non-floating-point operand is always false, so
//     the original code path reduced to "areEqual = (X[i] == X[currentIndex])".
template<typename T>
void RunLength(const vector<T>& X, vector<T>& B, vector<double>& N, vector<double>& BI) {
    B.clear();
    N.clear();
    BI.clear();
    if (X.empty()) return;

    B.push_back(X[0]);
    N.push_back(1);
    BI.push_back(0);

    if constexpr (std::is_floating_point_v<T>) {
        size_t currentIndex = 0;
        for (size_t i = 1; i < X.size(); ++i) {
            bool areEqual = false;
            const bool i_nan = std::isnan(X[i]);
            const bool c_nan = std::isnan(X[currentIndex]);
            if (i_nan && c_nan) {
                areEqual = true;
            }
            else if (!i_nan && !c_nan) {
                areEqual = (X[i] == X[currentIndex]);
            }

            if (areEqual) {
                N.back()++;
            }
            else {
                B.push_back(X[i]);
                N.push_back(1);
                BI.push_back(static_cast<double>(i));
                currentIndex = i;
            }
        }
    }
    else {
        // Non-floating-point fast path: no NaN, just compare to last run value.
        for (size_t i = 1; i < X.size(); ++i) {
            if (X[i] == B.back()) {
                N.back()++;
            }
            else {
                B.push_back(X[i]);
                N.push_back(1);
                BI.push_back(static_cast<double>(i));
            }
        }
    }
}

// Convenience overload returning a tuple
template<typename T>
tuple<vector<T>, vector<double>, vector<double>> RunLength(const vector<T>& X) {
    vector<T> B;
    vector<double> N, BI;
    RunLength(X, B, N, BI);
    return std::make_tuple(B, N, BI);
}

// Decode: expand (values, counts) back into a flat vector
template<typename T>
vector<T> RunLength(const vector<T>& B, const vector<double>& N) {
    vector<T> X;
    for (size_t i = 0; i < B.size(); ++i) {
        for (int j = 0; j < static_cast<int>(N[i]); ++j) {
            X.push_back(B[i]);
        }
    }
    return X;
}