// ============================================================================
// File: PeakFinder.h
// Find peaks in signal data with minimum distance constraint
// ============================================================================
#pragma once

#include <vector>

void findpeaks(const std::vector<double>& data,
    std::vector<double>& pks,
    std::vector<size_t>& locs,
    double minPeakDistance = 0);