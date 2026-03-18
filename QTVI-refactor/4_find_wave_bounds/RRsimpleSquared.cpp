// ============================================================================
// File: RRsimpleSquared.cpp
// Squares the ECG, then finds peaks above mean + 2*std with minimum distance.
// ============================================================================
#include "RRSimpleSquared.h"
#include "StatsUtils.h"
#include "PeakFinder.h"

pair<vector<size_t>, vector<double>> RRsimpleSquared(const vector<double>& ecg, double minDist) {
    if (ecg.size() < 3) return { {}, {} };

    // 1. Square the signal
    vector<double> ecgSigSq(ecg.size());
    for (size_t i = 0; i < ecg.size(); ++i) {
        ecgSigSq[i] = ecg[i] * ecg[i];
    }

    // 2. Threshold: mean + 2*std
    double threshold = mean(ecgSigSq) + std_dev(ecgSigSq) * 2.0;

    // 3. Find peaks with minimum distance
    vector<double> pks;
    vector<size_t> locs;
    findpeaks(ecgSigSq, pks, locs, minDist);

    // 4. Filter by threshold
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