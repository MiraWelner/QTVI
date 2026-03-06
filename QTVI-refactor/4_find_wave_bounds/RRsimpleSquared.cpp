// ============================================================================
// File: RRsimpleSquared.cpp
// Faithful translation of RRsimpleSquared.m
//
// MATLAB:
//   ecgSigSq = ecg.^2;
//   [rramps,rridx] = findpeaks(ecgSigSq, ...
//       'MinPeakHeight', mean(ecgSigSq)+std(ecgSigSq)*2, ...
//       'MinPeakDistance', minDist);
//   rridx = rridx';
// ============================================================================
#include "RRsimpleSquared.h"
#include "StatsUtils.h"
#include "PeakFinder.h"
#include <algorithm>

pair<vector<size_t>, vector<double>> RRsimpleSquared(const vector<double>& ecg, double minDist) {
    if (ecg.size() < 3) return { {}, {} };

    // 1. Square the signal
    vector<double> ecgSigSq(ecg.size());
    for (size_t i = 0; i < ecg.size(); ++i) {
        ecgSigSq[i] = ecg[i] * ecg[i];
    }

    // 2. Calculate threshold: mean + 2*std
    double meanVal = mean(ecgSigSq);
    double stdVal = std_dev(ecgSigSq);
    double threshold = meanVal + stdVal * 2.0;

    // 3. Use findpeaks (matches MATLAB's findpeaks with MinPeakDistance)
    vector<double> pks;
    vector<size_t> locs;

    try {
        findpeaks(ecgSigSq, pks, locs, minDist);
    }
    catch (...) {
        return { {}, {} };
    }

    // 4. Filter by MinPeakHeight threshold
    vector<size_t> rridx;
    vector<double> rramps;
    for (size_t i = 0; i < locs.size(); ++i) {
        if (pks[i] >= threshold) {
            rridx.push_back(locs[i]);
            rramps.push_back(pks[i]);
        }
    }

    return { rridx, rramps };
}