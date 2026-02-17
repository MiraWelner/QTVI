// ============================================================================
// File: PeakFinder.h
// Find peaks in signal data
// ============================================================================
#ifndef PEAKFINDER_H
#define PEAKFINDER_H

#include "SignalProcessingTypes.h"
#include <vector>

// Find peaks with minimum peak distance
void findpeaks(const std::vector<double>& data,
    std::vector<double>& pks,
    std::vector<size_t>& locs,
    double minPeakDistance = 0);

#endif // PEAKFINDER_H
