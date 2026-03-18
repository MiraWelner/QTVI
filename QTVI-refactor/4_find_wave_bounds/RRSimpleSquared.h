// ============================================================================
// File: RRSimpleSquared.h
// R-R peak detection using squared signal with threshold
// ============================================================================
#pragma once

#include "SignalProcessingTypes.h"

pair<vector<size_t>, vector<double>> RRsimpleSquared(const vector<double>& ecg, double minDist);