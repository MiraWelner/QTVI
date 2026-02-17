#include "PeakFinder.h"
#include <algorithm>

using namespace std;

struct PeakCandidate {
    double val;
    size_t pos;
};

void findpeaks(const vector<double>& data, vector<double>& pks, vector<size_t>& locs, double minPeakDistance) {
    pks.clear(); locs.clear();
    if (data.size() < 3) return;

    vector<PeakCandidate> candidates;

    // 1. Identify all local maxima (and center of plateaus)
    for (size_t i = 1; i < data.size() - 1; ++i) {
        if (data[i] > data[i - 1]) {
            size_t j = i;
            while (j < data.size() - 1 && data[j] == data[j + 1]) j++;

            // It's a peak if the next point after the plateau is lower
            if (j < data.size() - 1 && data[j] > data[j + 1]) {
                candidates.push_back({ data[i], i + (j - i) / 2 });
                i = j;
            }
        }
    }

    // 2. Priority Sort (Tallest peaks first - Essential for MATLAB matching)
    sort(candidates.begin(), candidates.end(), [](const PeakCandidate& a, const PeakCandidate& b) {
        if (a.val != b.val) return a.val > b.val;
        return a.pos < b.pos;
        });

    // 3. Elimination based on MinPeakDistance
    vector<bool> keep(candidates.size(), true);
    vector<PeakCandidate> final_peaks;

    for (size_t i = 0; i < candidates.size(); ++i) {
        if (!keep[i]) continue;
        final_peaks.push_back(candidates[i]);

        // Mark neighbors within distance as removed
        for (size_t j = i + 1; j < candidates.size(); ++j) {
            if (keep[j]) {
                long long dist = std::abs((long long)candidates[i].pos - (long long)candidates[j].pos);
                if (dist < (long long)minPeakDistance) {
                    keep[j] = false;
                }
            }
        }
    }

    // 4. Sort back to temporal order
    sort(final_peaks.begin(), final_peaks.end(), [](const PeakCandidate& a, const PeakCandidate& b) {
        return a.pos < b.pos;
        });

    for (const auto& p : final_peaks) {
        pks.push_back(p.val);
        locs.push_back(p.pos);
    }
}
