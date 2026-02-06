// ============================================================================
// File: nanfastsmooth.h
// Fast smoothing function that handles NaN values
// ============================================================================
#ifndef NANFASTSMOOTH_H
#define NANFASTSMOOTH_H

#include "SignalProcessingTypes.h"

// nanfastsmooth - smooths vector Y with moving average of width w ignoring NaNs
vector<double> nanfastsmooth(const vector<double>& Y, double w, int type = 1, double tol = 0.5);

#endif // NANFASTSMOOTH_H