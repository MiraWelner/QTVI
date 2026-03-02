// ============================================================================
// File: JoinedRR.h
// Joined R-R detection using multiple algorithms
// ============================================================================
#ifndef JOINEDRR_H
#define JOINEDRR_H

#include "SignalProcessingTypes.h"

// Joined R-R detection
vector<size_t> JoinedRR(const vector<double>& ecgSeg, double ecgSamplingRate, double diff_range, const std::string fileID);

#endif // JOINEDRR_H