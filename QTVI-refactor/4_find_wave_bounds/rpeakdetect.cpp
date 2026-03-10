#include "rpeakdetect.h"
#include "FilterUtils.h"
#include "StatsUtils.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include <numeric>
#include <fstream>


// Declarations to resolve C3861 (ensure these match the signatures in filterECG*.cpp)
std::vector<double> filterECG128Hz(const std::vector<double>& x);
std::vector<double> filterECG256Hz(const std::vector<double>& x);

// Internal helper for MATLAB-style diff (result size = N-1)
static std::vector<double> get_diff(const std::vector<double>& v) {
    if (v.size() < 2) return {};
    std::vector<double> res;
    res.reserve(v.size() - 1);
    for (size_t i = 0; i < v.size() - 1; ++i) {
        res.push_back(v[i + 1] - v[i]);
    }
    return res;
}

RPeakDetectResult rpeakdetect(const std::vector<double>& data, double samp_freq, double thresh, int testmode, std::string fileID) {
    RPeakDetectResult result;
    if (data.size() < 10) return result; // Safety check for very short data

    size_t n_orig = data.size();
    std::vector<double> x = data;

    // 1. Detrend: x = x - mean(x)
    double mu = mean(x); // Using mean from StatsUtils.h
    for (auto& val : x) val -= mu;

    // 2. Bandpass Filtering (Replicating MATLAB exist() logic)
    std::vector<double> bpf;
    if (std::abs(samp_freq - 128.0) < 0.1) {
        bpf = filterECG128Hz(x);
    }
    else if (std::abs(samp_freq - 256.0) < 0.1) {
        bpf = filterECG256Hz(x);
    }
    else {
        bpf = x;
    }

    // 3. Differentiation: dff = diff(bpf)
    std::vector<double> dff = get_diff(bpf);

    // 4. Squaring: sqr = dff .* dff
    std::vector<double> sqr(dff.size());
    for (size_t i = 0; i < dff.size(); ++i) {
        sqr[i] = dff[i] * dff[i];
    }
    size_t len = sqr.size(); // Corresponds to 'len = len - 1' in MATLAB

    // 5. Integration (Moving Sum)
    int win_size = (samp_freq >= 256.0) ? (int)std::round(7.0 * samp_freq / 256.0) : 7;
    std::vector<double> d_kernel(win_size, 1.0);
    // filter(d, 1, sqr)
    std::vector<double> filtered_sqr = filter(d_kernel, { 1.0 }, sqr);
    // mdfint = medfilt1(..., 10)
    std::vector<double> mdfint = medfilt1(filtered_sqr, 10);

    // 6. Remove Filter Delay: mdfint = mdfint(delay:length(mdfint))
    int delay = (int)std::ceil((double)win_size / 2.0);
    if (delay > 1 && (size_t)delay <= mdfint.size()) {
        mdfint.erase(mdfint.begin(), mdfint.begin() + (delay - 1));
    }

    // 7. Thresholding: NaN-Safe max_h calculation
    int start_search = (int)std::round((double)len / 4.0) - 1;
    int end_search = (int)std::round(3.0 * (double)len / 4.0) - 1;

    start_search = std::max(0, std::min(start_search, (int)mdfint.size() - 1));
    end_search = std::max(start_search, std::min(end_search, (int)mdfint.size() - 1));

    double max_h = -1e30; // Start with very small value
    bool found_valid = false;
    for (int i = start_search; i <= end_search; ++i) {
        if (std::isfinite(mdfint[i])) {
            if (mdfint[i] > max_h) {
                max_h = mdfint[i];
                found_valid = true;
            }
        }
    }
    // Safety: If the whole search range was invalid, default to a sensible threshold
    if (!found_valid) max_h = 0.0;

    // 8. Identify possible regions (Handle NaN by treating as 0)
    std::vector<int> poss_reg(mdfint.size(), 0);
    double limit = thresh * max_h;
    for (size_t i = 0; i < mdfint.size(); ++i) {
        if (!std::isnan(mdfint[i])) {
            poss_reg[i] = (mdfint[i] > limit) ? 1 : 0;
        }
    }

    // 9. Segment detection
    std::vector<int> left, right;
    if (!poss_reg.empty()) {
        if (poss_reg[0] == 1) left.push_back(1);
        for (size_t i = 1; i < (int)poss_reg.size(); ++i) {
            if (poss_reg[i] == 1 && poss_reg[i - 1] == 0) left.push_back(i + 1);
            if (poss_reg[i] == 0 && poss_reg[i - 1] == 1) right.push_back(i);
        }
        if (poss_reg.back() == 1) right.push_back((int)poss_reg.size());
    }

    // 10. Search for local Max/Min
    size_t num_segs = std::min(left.size(), right.size());
    std::vector<size_t> maxloc, minloc;
    std::vector<double> maxval, minval;

    for (size_t i = 0; i < num_segs; ++i) {
        int l = left[i] - 1;
        int r = right[i] - 1;

        double cur_max = -1e30, cur_min = 1e30;
        int cur_max_i = l, cur_min_i = l;

        for (int k = l; k <= r; ++k) {
            if (k >= 0 && k < (int)bpf.size()) {
                if (!std::isnan(bpf[k])) {
                    if (bpf[k] > cur_max) { cur_max = bpf[k]; cur_max_i = k; }
                    if (bpf[k] < cur_min) { cur_min = bpf[k]; cur_min_i = k; }
                }
            }
        }
        maxval.push_back(cur_max);
        maxloc.push_back((size_t)cur_max_i);
        minval.push_back(cur_min);
        minloc.push_back((size_t)cur_min_i);
    }

    // 11. Initial Assignment (MATLAB lines 162-167)
    result.R_index = maxloc; // R_index = maxloc
    result.R_amp = maxval;
    result.S_amp = minval;

    // Construct time vectors (t = 1/fs : 1/fs : end)
    auto get_time = [&](size_t idx) { return (double)(idx + 1) / samp_freq; };
    for (auto loc : maxloc) result.R_t.push_back(get_time(loc));
    for (auto loc : minloc) result.S_t.push_back(get_time(loc));


    // 12. Lead Inversion Check
    if (!maxloc.empty() && !minloc.empty()) {
        if (minloc.back() < maxloc.back()) {
            std::swap(result.R_t, result.S_t);
            std::swap(result.R_amp, result.S_amp);
        }
    }

    // 13. HRV calculation: hrv = diff(R_t)
    if (result.R_t.size() > 1) {
        for (size_t i = 0; i < result.R_t.size() - 1; ++i) {
            result.hrv.push_back(result.R_t[i + 1] - result.R_t[i]);
        }
    }

    return result;

}
