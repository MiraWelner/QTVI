// ============================================================================
// File: pairRtoPPGBeat.h
// Pair ECG R-peaks to PPG pulse valleys
// ============================================================================
#pragma once

#include "SignalProcessingTypes.h"

// Returns matrix where each row is [ppg_valley_idx, ecg_R_idx]
// -1.0 indicates unpaired
vector<vector<double>> pairRtoPPGBeat(const vector<double>& ecg,
    const vector<double>& ppg,
    double ecgSamplingRate,
    double ppgSamplingRate,
    const vector<size_t>& ecgRIndex,
    const vector<size_t>& ppgMinAmps);