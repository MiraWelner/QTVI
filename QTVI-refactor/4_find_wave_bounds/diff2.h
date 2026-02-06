// ============================================================================
// File: diff2.h
// Differentiate or difference by 4 points
// ============================================================================
#ifndef DIFF2_H
#define DIFF2_H

#include "SignalProcessingTypes.h"

// diff2 - Differentiate or difference by 4 points
// Modified from DIFF by Chen-Huan Chen 6/4/96
vector<double> diff2(const vector<double>& X, int nd = 1);

#endif // DIFF2_H