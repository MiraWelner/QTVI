// ============================================================================
// File: stdoutlier.h
// Detect outliers using standard deviation
// ============================================================================
#ifndef STDOUTLIER_H
#define STDOUTLIER_H

#include "SignalProcessingTypes.h"

// Detect outliers in data using moving mean and standard deviation
vector<bool> stdoutlier(const vector<double>& data,
    double multiplier,
    size_t mean_window,
    const string& direction,
    bool debug_plot = false);

#endif // STDOUTLIER_H