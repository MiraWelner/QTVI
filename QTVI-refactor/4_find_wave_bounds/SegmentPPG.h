// ============================================================================
// File: SegmentPPG.h
// Segment PPG signal to find minima and maxima
// ============================================================================
#ifndef SEGMENTPPG_H
#define SEGMENTPPG_H

#include "SignalProcessingTypes.h"

// Segment PPG signal
struct SegmentPPGResult {
    vector<size_t> ppgMinAmps;
    vector<size_t> maxAmps;
};

SegmentPPGResult SegmentPPG(const vector<double>& ppg, double sampleRate);

#endif // SEGMENTPPG_H
