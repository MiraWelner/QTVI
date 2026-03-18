// ============================================================================
// File: stdoutlier.h
// Detect outliers using moving mean and standard deviation of differences
// ============================================================================
#pragma once

#include "SignalProcessingTypes.h"

vector<bool> stdoutlier(const vector<double>& data,
    double multiplier,
    size_t mean_window,
    const string& direction,
    bool debug_plot = false);