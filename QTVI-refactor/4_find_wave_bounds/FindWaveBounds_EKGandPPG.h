// ============================================================================
// File: FindWaveBounds_EKGandPPG.h
// Find wave boundaries in ECG and PPG data
// ============================================================================
#ifndef FINDWAVEBOUNDS_EKGANDPPG_H
#define FINDWAVEBOUNDS_EKGANDPPG_H

#include "SignalProcessingTypes.h"

// Segment data structure
struct AnnealedSegment {
    vector<double> ecg;
    vector<double> po;  // PPG/pulse oximetry
    double ecgSampleRate;
    double ppgSampleRate;
    vector<size_t> r_peaks;
    vector<size_t> ppg_bin_indexs;
    vector<size_t> ecg_bin_indexs;
};

// Wave data structure
struct WaveData {
    vector<double> ecgSeg;
    vector<double> ppgSeg;
    vector<size_t> ecgRIndex;
    vector<size_t> ppgMinAmps;
    vector<size_t> ppgMaxAmps;
    vector<vector<double>> pairs;
    size_t index;
    double ecgSamplingRate;
    double ppgSamplingRate;
    vector<size_t> ppg_bin_indexs;
    vector<size_t> ecg_bin_indexs;
    bool bad_segment;
};

// Find wave bounds in ECG and PPG
vector<WaveData> FindWaveBounds_EKGandPPG(const vector<AnnealedSegment>& annealedSegments,
    int dbg_plot,
    bool use_R_algorithms);

#endif // FINDWAVEBOUNDS_EKGANDPPG_H