// ============================================================================
// File: RRSimpleSquared.hpp
// R-R peak detection using squared signal with threshold
// Squares the ECG, then finds peaks above mean + 2*std with minimum distance.
// ============================================================================
#pragma once

#include "StatsUtils.hpp"
#include "PeakFinder.hpp"

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <stdexcept>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <functional>
#include <tuple>

using std::vector;
using std::string;
using std::map;
using std::pair;
using std::tuple;

inline pair<vector<size_t>, vector<double>> RRsimpleSquared(const vector<double>& ecg, double minDist) {
    if (ecg.size() < 3) return { {}, {} };

    // 1. Square the signal
    vector<double> ecgSigSq(ecg.size());
    for (size_t i = 0; i < ecg.size(); ++i)
        ecgSigSq[i] = ecg[i] * ecg[i];

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