// ============================================================================
// File: ecgLms.h
// ECG LMS-based R-wave detection
// ============================================================================
#ifndef ECGLMS_H
#define ECGLMS_H

#include "SignalProcessingTypes.h"

// ECG LMS R-wave detection
vector<size_t> ecgLms(const vector<double>& ecg,
    int sampling,
    const vector<double>& b_butter_ecg4mwi,
    const vector<double>& a_butter_ecg4mwi,
    int dbg = 0);

#endif // ECGLMS_H
#pragma once
