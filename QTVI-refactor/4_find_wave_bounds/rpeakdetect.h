// ============================================================================
// File: rpeakdetect.h
// QRS detector based on Pan-Tompkins algorithm
// ============================================================================
#ifndef RPEAKDETECT_H
#define RPEAKDETECT_H

#include "SignalProcessingTypes.h"

// R peak detection - batch QRS detector
struct RPeakDetectResult {
    vector<size_t> R_index;
    vector<double> hrv;
    vector<double> R_t;
    vector<double> R_amp;
    vector<double> S_t;
    vector<double> S_amp;
};

RPeakDetectResult rpeakdetect(const vector<double>& data,
    double samp_freq = 256.0,
    double thresh = 0.2,
    int testmode = 0);

#endif // RPEAKDETECT_H
