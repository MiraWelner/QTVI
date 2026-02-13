// ============================================================================
// File: FindWaveBounds_EKGandPPG.h
// Find wave boundaries in ECG and PPG data
// ============================================================================
#ifndef FINDWAVEBOUNDS_EKGANDPPG_H
#define FINDWAVEBOUNDS_EKGANDPPG_H

#include <vector>
#include "FindWaveBounds.h"  // This brings in AnnealedSegment and WaveData

// Find wave bounds in ECG and PPG
std::vector<WaveData> FindWaveBounds_EKGandPPG(const std::vector<AnnealedSegment>& annealedSegments,
    int dbg_plot,
    bool use_R_algorithms);

#endif // FINDWAVEBOUNDS_EKGANDPPG_H
