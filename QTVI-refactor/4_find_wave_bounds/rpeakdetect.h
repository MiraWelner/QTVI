// ============================================================================
// File: rpeakdetect.h
// R-peak detection using bandpass filtering, differentiation, and thresholding
// ============================================================================
#pragma once

#include <vector>
#include <cstddef>
#include <string>

using namespace std;

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
    int testmode = 0,
    std::string fileID = "");