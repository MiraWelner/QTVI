// ============================================================================
// File: RRsimpleSquared.h
// Simple R-R peak detection using squared signal
// ============================================================================
#ifndef RRSIMPLESQUARED_H
#define RRSIMPLESQUARED_H

#include "SignalProcessingTypes.h"

// Simple R-R peak detection
pair<vector<size_t>, vector<double>> RRsimpleSquared(const vector<double>& ecg, double minDist);

#endif // RRSIMPLESQUARED_H