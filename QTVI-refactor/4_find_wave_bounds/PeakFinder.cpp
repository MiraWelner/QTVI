#include "PeakFinder.h"
#include <algorithm>

struct Peak {
    double val;
    size_t pos;
};

void findpeaks(const vector<double>& data, vector<double>& pks, vector<size_t>& locs, double minPeakDistance) {
    pks.clear(); locs.clear();
    if (data.size() < 3) return;

    vector<Peak> candidates;
    // 1. Identify all local maxima (handling plateaus)
    for (size_t i = 1; i < data.size() - 1; ++i) {
        if (data[i] > data[i - 1]) {
            size_t j = i;
            while (j < data.size() - 1 && data[j] == data[j + 1]) j++;
            if (j < data.size() - 1 && data[j] > data[j + 1]) {
                candidates.push_back({ data[i], i + (j - i) / 2 });
                i = j;
            }
        }
    }

    // 2. Priority Sort (Tallest peaks first) - This is what MATLAB does
    sort(candidates.begin(), candidates.end(), [](const Peak& a, const Peak& b) {
        return a.val > b.val;
        });

    // 3. Elimination based on MinPeakDistance
    vector<bool> keep(candidates.size(), true);
    vector<Peak> final_peaks;

    for (size_t i = 0; i < candidates.size(); ++i) {
        if (!keep[i]) continue;
        final_peaks.push_back(candidates[i]);

        for (size_t j = i + 1; j < candidates.size(); ++j) {
            if (keep[j] && abs((long long)candidates[i].pos - (long long)candidates[j].pos) < minPeakDistance) {
                keep[j] = false;
            }
        }
    }

    // 4. Sort back to temporal order
    sort(final_peaks.begin(), final_peaks.end(), [](const Peak& a, const Peak& b) {
        return a.pos < b.pos;
        });

    for (auto& p : final_peaks) {
        pks.push_back(p.val);
        locs.push_back(p.pos);
    }
}
