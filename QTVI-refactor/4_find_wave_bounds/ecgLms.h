// ============================================================================
// File: ecgLms.h
// ECG LMS-based R-wave detection
// ============================================================================
#pragma once

#include "SignalProcessingTypes.h"

vector<size_t> ecgLms(const vector<double>& ecg,
    int sampling,
    const vector<double>& b_butter_ecg4mwi,
    const vector<double>& a_butter_ecg4mwi,
    int dbg = 0,
    const std::string& fileID = "");