// ============================================================================
// File: PeakFinder.cpp
// ============================================================================
#include "PeakFinder.h"

void findpeaks(const vector<double>& data,
    vector<double>& pks,
    vector<size_t>& locs,
    double minPeakDistance) {
    pks.clear();
    locs.clear();

    if (data.size() < 3) return;

    // Find all local maxima
    for (size_t i = 1; i < data.size() - 1; ++i) {
        if (data[i] > data[i - 1] && data[i] > data[i + 1]) {
            pks.push_back(data[i]);
            locs.push_back(i);
        }
    }

    // Apply minimum peak distance filter
    if (minPeakDistance > 0 && !locs.empty()) {
        vector<double> filtered_pks;
        vector<size_t> filtered_locs;

        filtered_pks.push_back(pks[0]);
        filtered_locs.push_back(locs[0]);

        for (size_t i = 1; i < locs.size(); ++i) {
            if (locs[i] - filtered_locs.back() >= minPeakDistance) {
                filtered_pks.push_back(pks[i]);
                filtered_locs.push_back(locs[i]);
            }
            else {
                // Keep the larger peak
                if (pks[i] > filtered_pks.back()) {
                    filtered_pks.back() = pks[i];
                    filtered_locs.back() = locs[i];
                }
            }
        }

        pks = filtered_pks;
        locs = filtered_locs;
    }
}