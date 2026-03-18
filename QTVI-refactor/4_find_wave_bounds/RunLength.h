// ============================================================================
// File: RunLength.h
// Run-length encoding and decoding (header-only, templated)
// ============================================================================
#pragma once

#include "SignalProcessingTypes.h"

// Encode: consecutive equal values are grouped into (value, count, start_index)
template<typename T>
void RunLength(const vector<T>& X, vector<T>& B, vector<double>& N, vector<double>& BI) {
    B.clear();
    N.clear();
    BI.clear();
    if (X.empty()) return;

    B.push_back(X[0]);
    N.push_back(1);
    BI.push_back(0);

    size_t currentIndex = 0;

    for (size_t i = 1; i < X.size(); ++i) {
        bool areEqual = false;
        if (std::isnan(static_cast<double>(X[i])) && std::isnan(static_cast<double>(X[currentIndex]))) {
            areEqual = true;
        }
        else if (!std::isnan(static_cast<double>(X[i])) && !std::isnan(static_cast<double>(X[currentIndex]))) {
            areEqual = (X[i] == X[currentIndex]);
        }

        if (areEqual) {
            N.back()++;
        }
        else {
            B.push_back(X[i]);
            N.push_back(1);
            BI.push_back(i);
            currentIndex = i;
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