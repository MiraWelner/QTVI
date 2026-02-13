// ============================================================================
// File: rpeakdetect.cpp
// ============================================================================
#include "rpeakdetect.h"
#include "FilterUtils.h"
#include "StatsUtils.h"

RPeakDetectResult rpeakdetect(const vector<double>& data, double samp_freq, double thresh, int testmode) {
    RPeakDetectResult result;

    size_t len = data.size();

    // Make time axis
    vector<double> t(len);
    for (size_t i = 0; i < len; ++i) {
        t[i] = (i + 1) / samp_freq;
    }

    vector<double> x = data;

    // Remove mean
    double meanVal = mean(x);
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] -= meanVal;
    }

    // Bandpass filter data
    vector<double> bpf = x;

    // Simplified filtering for different sampling rates
    if (samp_freq == 256.0 || samp_freq == 128.0) {
        // Low pass filter
        double Wn = 12.0 * 2.0 / samp_freq;
        int N = 3;
        vector<double> a, b;
        butter(N, Wn, "low", b, a);
        bpf = filtfilt(b, a, x);

        // Normalize
        double maxAbsBpf = 0.0;
        for (const auto& val : bpf) {
            maxAbsBpf = std::max(maxAbsBpf, std::abs(val));
        }
        if (maxAbsBpf > 0) {
            for (auto& val : bpf) {
                val /= maxAbsBpf;
            }
        }

        // High pass filter
        Wn = 5.0 * 2.0 / samp_freq;
        butter(N, Wn, "high", b, a);
        bpf = filtfilt(b, a, bpf);

        // Normalize again
        maxAbsBpf = 0.0;
        for (const auto& val : bpf) {
            maxAbsBpf = std::max(maxAbsBpf, std::abs(val));
        }
        if (maxAbsBpf > 0) {
            for (auto& val : bpf) {
                val /= maxAbsBpf;
            }
        }
    }
    else {
        // Bandpass filter for other sampling frequencies
        double f1 = 5.0;
        double f2 = 15.0;
        vector<double> Wn = { f1 * 2.0 / samp_freq, f2 * 2.0 / samp_freq };
        int N = 3;
        vector<double> a, b;
        butter(N, Wn, b, a);
        bpf = filtfilt(b, a, x);

        // Normalize
        double maxAbsBpf = 0.0;
        for (const auto& val : bpf) {
            maxAbsBpf = std::max(maxAbsBpf, std::abs(val));
        }
        if (maxAbsBpf > 0) {
            for (auto& val : bpf) {
                val /= maxAbsBpf;
            }
        }
    }

    // Differentiate
    vector<double> dff = diff(bpf);

    // Square
    vector<double> sqr(dff.size());
    for (size_t i = 0; i < dff.size(); ++i) {
        sqr[i] = dff[i] * dff[i];
    }

    len = sqr.size();

    // Integrate data over window
    int windowSize = 7;
    if (samp_freq >= 256) {
        windowSize = static_cast<int>(std::round(7 * samp_freq / 256.0));
    }

    vector<double> d(windowSize, 1.0);
    vector<double> integrated = filter(d, { 1.0 }, sqr);
    vector<double> mdfint = medfilt1(integrated, 10);

    // Remove filter delay
    int delay = windowSize / 2;
    vector<double> mdfint_delayed(mdfint.begin() + delay, mdfint.end());

    // Find highest bumps
    double max_h = -1.0;
    for (size_t i = 0; i < mdfint_delayed.size(); ++i) {
        if (mdfint_delayed[i] > max_h) {
            max_h = mdfint_delayed[i];
        }
    }

    // Find regions above threshold
    vector<bool> poss_reg(mdfint_delayed.size(), false);
    for (size_t i = 0; i < mdfint_delayed.size(); ++i) {
        poss_reg[i] = mdfint_delayed[i] > (thresh * max_h);
    }

    // Find boundaries
    vector<size_t> left, right;
    bool inRegion = false;

    for (size_t i = 0; i < poss_reg.size(); ++i) {
        if (poss_reg[i] && !inRegion) {
            left.push_back(i);
            inRegion = true;
        }
        else if (!poss_reg[i] && inRegion) {
            right.push_back(i);
            inRegion = false;
        }
    }

    if (inRegion) {
        right.push_back(poss_reg.size());
    }

    // Find peaks in each region
    vector<double> maxval, minval;
    vector<size_t> maxloc, minloc;

    for (size_t i = 0; i < left.size() && i < right.size(); ++i) {
        auto maxResult = max_element_index(bpf, left[i], right[i]);
        auto minResult = min_element_index(bpf, left[i], right[i]);

        maxval.push_back(maxResult.first);
        maxloc.push_back(left[i] + maxResult.second);

        minval.push_back(minResult.first);
        minloc.push_back(left[i] + minResult.second);
    }

    result.R_index = maxloc;
    result.R_amp = maxval;
    result.S_amp = minval;

    result.R_t.resize(maxloc.size());
    result.S_t.resize(minloc.size());

    for (size_t i = 0; i < maxloc.size(); ++i) {
        if (maxloc[i] < t.size()) {
            result.R_t[i] = t[maxloc[i]];
        }
    }

    for (size_t i = 0; i < minloc.size(); ++i) {
        if (minloc[i] < t.size()) {
            result.S_t[i] = t[minloc[i]];
        }
    }

    // Check for lead inversion
    if (!minloc.empty() && !maxloc.empty()) {
        if (minloc.back() < maxloc.back()) {
            result.R_t = result.S_t;
            result.R_amp = result.S_amp;
            result.S_amp = maxval;
            result.S_t.resize(maxloc.size());
            for (size_t i = 0; i < maxloc.size(); ++i) {
                if (maxloc[i] < t.size()) {
                    result.S_t[i] = t[maxloc[i]];
                }
            }
        }
    }

    // Calculate HRV
    result.hrv = diff(result.R_t);

    return result;
}