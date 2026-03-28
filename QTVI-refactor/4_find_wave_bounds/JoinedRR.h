// ============================================================================
// File: JoinedRR.h
// Ensemble R-R detection using multiple weighted algorithms
// ============================================================================
#pragma once

#include "SignalProcessingTypes.h"

vector<size_t> JoinedRR(const vector<double>& ecgSeg, const std::string fileID);