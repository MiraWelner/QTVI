// ============================================================================
// File: PeakFinder.h
// Find peaks in signal data
// ============================================================================
#ifndef PEAKFINDER_H
#define PEAKFINDER_H

#include "SignalProcessingTypes.h"

// Find peaks with minimum peak distance
void findpeaks(const vector<double>& data,
    vector<double>& pks,
    vector<size_t>& locs,
    double minPeakDistance = 0);

#endif // PEAKFINDER_H