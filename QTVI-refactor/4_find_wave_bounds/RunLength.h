// ============================================================================
// File: RunLength.h
// Run-length encoding and decoding
// ============================================================================
#ifndef RUNLENGTH_H
#define RUNLENGTH_H

#include "SignalProcessingTypes.h"

// RunLength encoding function
template<typename T>
void RunLength(const vector<T>& X, vector<T>& B, vector<double>& N, vector<double>& BI) {
    if (X.empty()) {
        B.clear();
        N.clear();
        BI.clear();
        return;
    }

    B.clear();
    N.clear();
    BI.clear();

    // First element
    B.push_back(X[0]);
    N.push_back(1);
    BI.push_back(0);  // 0-indexed in C++

    size_t currentIndex = 0;

    for (size_t i = 1; i < X.size(); ++i) {
        // Compare with NaN handling
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

// Overload for common case with just 3 outputs
template<typename T>
tuple<vector<T>, vector<double>, vector<double>> RunLength(const vector<T>& X) {
    vector<T> B;
    vector<double> N;
    vector<double> BI;
    RunLength(X, B, N, BI);
    return std::make_tuple(B, N, BI);
}

// RunLength decoding function
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

#endif // RUNLENGTH_H