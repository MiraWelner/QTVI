#include "rpeakdetect.h"
#include "FilterUtils.h"
#include "StatsUtils.h"
#include <cmath>
#include <algorithm>
#include <vector>

using namespace std;

RPeakDetectResult rpeakdetect(const vector<double>& data, double samp_freq, double thresh, int testmode) {
    RPeakDetectResult result;
    if (data.empty()) return result;

    size_t len = data.size();
    vector<double> x = data;

    // 1. Remove Mean
    double mu = mean(x);
    for (auto& val : x) val -= mu;

    // 2. Bandpass Filtering Stage (FIR/IIR)
    // Note: MATLAB script uses external functions filterECG128Hz/256Hz if they exist.
    // Here we use the signal 'x' as 'bpf' as a baseline or apply a 5-15Hz BPF.
    vector<double> bpf = x;
    // (If specific filter functions are needed, they should be called here)

    // 3. Differentiate data
    // MATLAB: dff = diff(bpf); -> result is 1 datum shorter
    vector<double> dff = diff(bpf);

    // 4. Square data
    vector<double> sqr(dff.size());
    for (size_t i = 0; i < dff.size(); ++i) {
        sqr[i] = dff[i] * dff[i];
    }
    size_t new_len = sqr.size();

    // 5. Integrate data over window 'd'
    // MATLAB: d = [ones(1, round(7 * samp_freq / 256))];
    int win_size = (int)std::round(7.0 * samp_freq / 256.0);
    if (win_size < 1) win_size = 1;
    vector<double> d_kernel(win_size, 1.0);

    // MATLAB: filter(d, 1, sqr) -> This is a moving sum (not averaged)
    vector<double> integrated = filter(d_kernel, { 1.0 }, sqr);

    // 6. Median filter
    // MATLAB: mdfint = medfilt1(..., 10);
    vector<double> mdfint = medfilt1(integrated, 10);

    // 7. Remove filter delay
    // MATLAB: delay = ceil(length(d) / 2); mdfint = mdfint(delay:end);
    int delay = (int)std::ceil((double)win_size / 2.0);
    if ((size_t)delay <= mdfint.size()) {
        // Adjust for 0-based indexing: MATLAB delay:end is (delay-1) index in C++
        mdfint.erase(mdfint.begin(), mdfint.begin() + (delay - 1));
    }

    // 8. Find highest bumps in the middle 50%
    // MATLAB: max_h = max(mdfint(round(len / 4):round(3 * len / 4)));
    size_t start_h = (size_t)std::round((double)len / 4.0);
    size_t end_h = (size_t)std::round(3.0 * (double)len / 4.0);
    double max_h = 0.0;
    for (size_t i = start_h; i < end_h && i < mdfint.size(); ++i) {
        if (mdfint[i] > max_h) max_h = mdfint[i];
    }

    // 9. Build segment array (poss_reg)
    // MATLAB: poss_reg = mdfint > (thresh * max_h);
    vector<bool> poss_reg(mdfint.size());
    double threshold_val = thresh * max_h;
    for (size_t i = 0; i < mdfint.size(); ++i) {
        poss_reg[i] = (mdfint[i] > threshold_val);
    }

    // 10. Find boundaries (Left/Right)
    // MATLAB: left = find(diff([0 poss_reg']) == 1);
    // MATLAB: right = find(diff([poss_reg' 0]) == -1);
    vector<size_t> left, right;
    for (size_t i = 0; i < poss_reg.size(); ++i) {
        bool prev = (i == 0) ? false : poss_reg[i - 1];
        bool current = poss_reg[i];
        bool next = (i == poss_reg.size() - 1) ? false : poss_reg[i + 1];

        if (current && !prev) left.push_back(i);
        if (current && !next) right.push_back(i);
    }

    // 11. Loop through segments to find peaks in BPF signal
    vector<double> maxval, minval;
    vector<size_t> maxloc, minloc;

    for (size_t i = 0; i < left.size() && i < right.size(); ++i) {
        size_t L = left[i];
        size_t R = right[i];

        // Bounds check
        if (L >= bpf.size()) continue;
        if (R >= bpf.size()) R = bpf.size() - 1;

        // Find max/min in bpf(left(i):right(i))
        auto max_res = max_element_index(bpf, L, R + 1);
        auto min_res = min_element_index(bpf, L, R + 1);

        maxval.push_back(max_res.first);
        maxloc.push_back(L + max_res.second);
        minval.push_back(min_res.first);
        minloc.push_back(L + min_res.second);
    }

    // 12. Handle Lead Inversion
    // MATLAB: if (minloc(end) < maxloc(end))
    if (!maxloc.empty() && !minloc.empty()) {
        if (minloc.back() < maxloc.back()) {
            result.R_index = minloc;
            result.R_amp = minval;
            result.S_amp = maxval;
            // S_t will use maxloc
            for (size_t idx : maxloc) result.S_t.push_back((double)(idx + 1) / samp_freq);
        }
        else {
            result.R_index = maxloc;
            result.R_amp = maxval;
            result.S_amp = minval;
            // S_t will use minloc
            for (size_t idx : minloc) result.S_t.push_back((double)(idx + 1) / samp_freq);
        }
    }

    // 13. Populate R_t and HRV
    // MATLAB: R_t = t(maxloc); (where t = 1/fs : 1/fs : len/fs)
    for (size_t idx : result.R_index) {
        result.R_t.push_back((double)(idx + 1) / samp_freq);
    }

    // MATLAB: hrv = diff(R_t);
    if (result.R_t.size() > 1) {
        result.hrv = diff(result.R_t);
    }

    return result;
}
