// ============================================================================
// File: nanfastsmooth.h
// Fast smoothing with moving average, ignoring NaN values
// ============================================================================
#pragma once

#include "SignalProcessingTypes.h"

vector<double> nanfastsmooth(const vector<double>& Y, double w, int type = 1, double tol = 0.5);