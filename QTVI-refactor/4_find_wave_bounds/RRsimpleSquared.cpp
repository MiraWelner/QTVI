// ============================================================================
// File: RRsimpleSquared.cpp
// ============================================================================
#include "RRsimpleSquared.h"
#include "StatsUtils.h"
#include "PeakFinder.h"

pair<vector<size_t>, vector<double>> RRsimpleSquared(const vector<double>& ecg, double minDist) {
    // Square the signal
    vector<double> ecgSigSq(ecg.size());
    for (size_t i = 0; i < ecg.size(); ++i) {
        ecgSigSq[i] = ecg[i] * ecg[i];
    }

    // Calculate threshold
    double meanVal = mean(ecgSigSq);
    double stdVal = std_dev(ecgSigSq);
    double threshold = meanVal + stdVal * 2;

    // Find peaks
    vector<double> rramps;
    vector<size_t> rridx;

    try {
        // Find all peaks above threshold
        for (size_t i = 1; i < ecgSigSq.size() - 1; ++i) {
            if (ecgSigSq[i] > ecgSigSq[i - 1] &&
                ecgSigSq[i] > ecgSigSq[i + 1] &&
                ecgSigSq[i] > threshold) {

                // Check minimum distance
                if (rridx.empty() || (i - rridx.back()) >= minDist) {
                    rramps.push_back(ecgSigSq[i]);
                    rridx.push_back(i);
                }
            }
        }
    }
    catch (...) {
        rridx.clear();
        rramps.clear();
    }

    return { rridx, rramps };
}
