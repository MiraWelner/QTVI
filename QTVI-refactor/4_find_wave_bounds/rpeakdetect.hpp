// ============================================================================
// File: rpeakdetect.hpp
// R-peak detection using bandpass filtering, differentiation, and thresholding
// ============================================================================
#pragma once

#include <vector>
#include <cstddef>
#include <string>
#include <cmath>
#include <algorithm>
#include <numeric>

#include "bandpass.h"
#include "FilterUtils.h"
#include "StatsUtils.h"

using namespace std;

struct RPeakDetectResult {
    vector<size_t> R_index;
    vector<double> hrv;
    vector<double> R_t;
    vector<double> R_amp;
    vector<double> S_t;
    vector<double> S_amp;
};

inline RPeakDetectResult rpeakdetect(const vector<double>& data,
    double samp_freq = 256.0,
    double thresh = 0.2,
    int testmode = 0,
    std::string fileID = "") {
    RPeakDetectResult result;
    if (data.size() < 10) return result;

    size_t n_orig = data.size();
    std::vector<double> x = data;

    // 1. Detrend: x = x - mean(x)
    double mu = mean(x);
    for (auto& val : x) val -= mu;

    // 2. Bandpass filtering
    std::vector<double> bpf = bandpass_filtfilt(2, 0.05, 100.0, 2000.0, x);

    // 3. Differentiation
    std::vector<double> dff = diff(bpf);

    // 4. Squaring
    std::vector<double> sqr(dff.size());
    for (size_t i = 0; i < dff.size(); ++i) {
        sqr[i] = dff[i] * dff[i];
    }
    size_t len = sqr.size();

    // 5. Integration (moving sum + median filter)
    int win_size = (samp_freq >= 256.0) ? (int)std::round(7.0 * samp_freq / 256.0) : 7;
    std::vector<double> d_kernel(win_size, 1.0);
    std::vector<double> filtered_sqr = filter(d_kernel, { 1.0 }, sqr);
    std::vector<double> mdfint = medfilt1(filtered_sqr, 10);

    // 6. Remove filter delay
    int delay = (int)std::ceil((double)win_size / 2.0);
    if (delay > 1 && (size_t)delay <= mdfint.size()) {
        mdfint.erase(mdfint.begin(), mdfint.begin() + (delay - 1));
    }

    // 7. Threshold from middle 50% of signal
    int start_search = std::max(0, std::min((int)std::round((double)len / 4.0) - 1, (int)mdfint.size() - 1));
    int end_search = std::max(start_search, std::min((int)std::round(3.0 * (double)len / 4.0) - 1, (int)mdfint.size() - 1));

    double max_h = 0.0;
    for (int i = start_search; i <= end_search; ++i) {
        if (std::isfinite(mdfint[i]) && mdfint[i] > max_h) {
            max_h = mdfint[i];
        }
    }

    // 8. Identify candidate regions above threshold
    std::vector<int> poss_reg(mdfint.size(), 0);
    double limit = thresh * max_h;
    for (size_t i = 0; i < mdfint.size(); ++i) {
        if (!std::isnan(mdfint[i]) && mdfint[i] > limit) {
            poss_reg[i] = 1;
        }
    }

    // 9. Find region boundaries (rising/falling edges)
    std::vector<int> left, right;
    if (!poss_reg.empty()) {
        if (poss_reg[0] == 1) left.push_back(1);
        for (size_t i = 1; i < poss_reg.size(); ++i) {
            if (poss_reg[i] == 1 && poss_reg[i - 1] == 0) left.push_back(i + 1);
            if (poss_reg[i] == 0 && poss_reg[i - 1] == 1) right.push_back(i);
        }
        if (poss_reg.back() == 1) right.push_back((int)poss_reg.size());
    }

    // 10. Find local max within each region (R-peaks)
    size_t num_segs = std::min(left.size(), right.size());
    std::vector<size_t> maxloc;
    std::vector<double> maxval;

    for (size_t i = 0; i < num_segs; ++i) {
        int l = left[i] - 1;
        int r = right[i] - 1;

        double cur_max = -1e30;
        int cur_max_i = l;

        for (int k = l; k <= r; ++k) {
            if (k >= 0 && k < (int)bpf.size() && !std::isnan(bpf[k])) {
                if (bpf[k] > cur_max) { cur_max = bpf[k]; cur_max_i = k; }
            }
        }
        maxval.push_back(cur_max);
        maxloc.push_back((size_t)cur_max_i);
    }

    // 11. Assign results (polarity is handled upstream, always use maxima as R-peaks)
    result.R_index = maxloc;
    result.R_amp = maxval;

    auto get_time = [&](size_t idx) { return (double)(idx + 1) / samp_freq; };
    for (auto loc : maxloc) result.R_t.push_back(get_time(loc));

    // 12. HRV: diff of R-peak times
    if (result.R_t.size() > 1) {
        for (size_t i = 0; i < result.R_t.size() - 1; ++i) {
            result.hrv.push_back(result.R_t[i + 1] - result.R_t[i]);
        }
    }

    return result;
}