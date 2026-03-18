// ============================================================================
// File: SegmentPPG.h
// Segment PPG signal into individual pulses (minima and maxima)
// ============================================================================
#pragma once

#include "SignalProcessingTypes.h"

struct SegmentPPGResult {
    vector<size_t> ppgMinAmps;  // valley (trough) indices
    vector<size_t> maxAmps;     // peak indices
};

SegmentPPGResult SegmentPPG(const vector<double>& ppg, double sampleRate);