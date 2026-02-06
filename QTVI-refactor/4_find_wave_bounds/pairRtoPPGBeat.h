// ============================================================================
// File: pairRtoPPGBeat.h
// Pair ECG R-peaks to PPG beats
// ============================================================================
#ifndef PAIRRTOPPGBEAT_H
#define PAIRRTOPPGBEAT_H

#include "SignalProcessingTypes.h"

// Pair R-peaks to PPG beats
// Returns matrix where each row is [ppg_valley_idx, ecg_R_idx]
vector<vector<double>> pairRtoPPGBeat(const vector<double>& ecg,
    const vector<double>& ppg,
    double ecgSamplingRate,
    double ppgSamplingRate,
    const vector<size_t>& ecgRIndex,
    const vector<size_t>& ppgMinAmps);

#endif // PAIRRTOPPGBEAT_H#pragma once
