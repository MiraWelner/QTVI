// ============================================================================
// File: ecglaux.h
// ECG R-wave detection auxiliary function
// ============================================================================
#ifndef ECGLAUX_H
#define ECGLAUX_H

#include "SignalProcessingTypes.h"

// ECG R-wave detection using Moving Window Integration
tuple<vector<size_t>, double, double> ecglaux(
    const vector<double>& ecg,
    const vector<double>& mwisignal,
    int sampling,
    double mwithold,
    double mvimaxval,
    int mwiwidthpts,
    int refractpts,
    double mwitholdfract,
    double mwitholdff
);
#endif // ECGLAUX_H
