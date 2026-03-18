// ============================================================================
// File: ecglaux.h
// ECG R-wave detection using Moving Window Integration
// ============================================================================
#pragma once

#include "SignalProcessingTypes.h"

tuple<vector<size_t>, double, double> ecglaux(
    const vector<double>& ecg,
    const vector<double>& mwisignal,
    int sampling,
    double mwithold,
    double mvimaxval,
    int mwiwidthpts,
    int refractpts,
    double mwitholdfract,
    double mwitholdff);