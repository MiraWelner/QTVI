// ============================================================================
// File: pan_tompkin.h
// Complete implementation of Pan-Tompkins algorithm
// ============================================================================
#ifndef PAN_TOMPKIN_H
#define PAN_TOMPKIN_H

#include "SignalProcessingTypes.h"

// Pan-Tompkins QRS detection result
struct PanTompkinResult {
    vector<size_t> qrs_i_raw;
    vector<double> qrs_amp_raw;
    int delay;
};

PanTompkinResult pan_tompkin(const vector<double>& ecg, double fs, int gr = 0);

#endif // PAN_TOMPKIN_H

// ============================================================================
// File: pan_tompkin.cpp   `  
// Complete implement