// ============================================================================
// File: rpeakdetect.cpp
// ============================================================================
#include "rpeakdetect.h"
#include "FilterUtils.h"
#include "StatsUtils.h"
#include <cmath>   // For std::isnan
#include <limits>  // For std::numeric_limits

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
        bpf = filtfilt(a, b, x);

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
        bpf = filtfilt(a, b, bpf);

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
        bpf = filtfilt(a, b, x);

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
    size_t start_idx_for_max_h = len / 4;
    size_t end_idx_for_max_h = 3 * len / 4;
    // Ensure indices are within bounds for mdfint_delayed
    if (start_idx_for_max_h >= mdfint_delayed.size()) start_idx_for_max_h = 0;
    if (end_idx_for_max_h > mdfint_delayed.size()) end_idx_for_max_h = mdfint_delayed.size();


    double max_h = -std::numeric_limits<double>::infinity(); // Use negative infinity
    for (size_t i = start_idx_for_max_h; i < end_idx_for_max_h; ++i) {
        if (mdfint_delayed[i] > max_h) {
            max_h = mdfint_delayed[i];
        }
    }
    // Handle case where max_h might still be -infinity (e.g., empty range or all NaNs)
    if (std::isinf(max_h)) {
        max_h = 0.0; // Or handle as an error, depending on expected behavior
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

    // Corrected section for finding peaks in each region
    // The mdfint_delayed signal has undergone differentiation, squaring, integration,
    // and truncation, all of which introduce a temporal shift/delay relative to the original bpf signal.
    // The calculated offset aims to align the indices from mdfint_delayed (used for 'left' and 'right')
    // with the corresponding time points in the bpf signal.
    //
    // Detailed derivation of temporal_offset:
    // 1. diff operation: introduces ~0.5 sample shift (dff[k] aligns with bpf[k + 0.5])
    // 2. filter(d, {1.0}, sqr) (causal FIR filter of length W=windowSize): group delay is (W-1)/2.
    // 3. medfilt1(integrated, 10): group delay is (10-1)/2 = 4.5.
    // 4. mdfint_delayed truncation: compensates for (windowSize / 2) samples (C++ integer division).
    //
    // Summing these delays:
    // Effective alignment point in bpf for mdfint[k]: k - (W-1)/2 - 4.5 + 0.5 = k - W/2 + 0.5 - 4.5 + 0.5 = k - W/2 - 3.5
    // Since mdfint_delayed[idx] corresponds to mdfint[idx + W/2] (due to truncation):
    // Effective alignment point in bpf for mdfint_delayed[idx]: (idx + W/2) - W/2 - 3.5 = idx - 3.5
    //
    // Thus, an index 'idx' in mdfint_delayed corresponds to an index 'idx - 3.5' in bpf.
    // To get integer indices, we use a floor or ceil. Rounding -3.5 to -4 seems appropriate.
    // This offset shifts the search window in bpf earlier by approximately 4 samples.
    int temporal_offset = -4; // Derived offset

    // Find peaks in each region
    vector<double> maxval_peaks, minval_s_waves; // Renamed for clarity
    vector<size_t> maxloc_r_index, minloc_s_index; // Renamed for clarity

    for (size_t i = 0; i < left.size() && i < right.size(); ++i) {
        // Adjust segment boundaries for bpf
        // Ensure that the calculated segment boundaries are valid and within bpf bounds
        size_t bpf_segment_start_idx = 0;
        if ((int)left[i] + temporal_offset > 0) {
            bpf_segment_start_idx = (size_t)((int)left[i] + temporal_offset);
        }
        else {
            bpf_segment_start_idx = 0; // Clamp to 0
        }

        size_t bpf_segment_end_idx = bpf.size();
        if ((int)right[i] + temporal_offset < (int)bpf.size()) {
            bpf_segment_end_idx = (size_t)((int)right[i] + temporal_offset);
        }
        else {
            bpf_segment_end_idx = bpf.size(); // Clamp to bpf.size()
        }

        // Ensure the segment is valid (start < end)
        if (bpf_segment_start_idx >= bpf_segment_end_idx) {
            continue; // Skip invalid or empty segments
        }

        auto maxResult = max_element_index(bpf, bpf_segment_start_idx, bpf_segment_end_idx);
        auto minResult = min_element_index(bpf, bpf_segment_start_idx, bpf_segment_end_idx);

        // Check for NaN results which might indicate an invalid segment or no valid numbers
        if (std::isnan(maxResult.first) || std::isnan(minResult.first)) {
            continue;
        }

        maxval_peaks.push_back(maxResult.first);
        maxloc_r_index.push_back(maxResult.second); // maxResult.second is the absolute index in bpf

        minval_s_waves.push_back(minResult.first);
        minloc_s_index.push_back(minResult.second); // minResult.second is the absolute index in bpf
    }

    result.R_index = maxloc_r_index;
    result.R_amp = maxval_peaks;
    result.S_amp = minval_s_waves;

    result.R_t.resize(result.R_index.size());
    result.S_t.resize(result.S_amp.size()); // S_amp and S_t should have same size as minloc_s_index

    for (size_t i = 0; i < result.R_index.size(); ++i) {
        if (result.R_index[i] < t.size()) {
            result.R_t[i] = t[result.R_index[i]];
        }
    }

    // Assuming S_t should correspond to minloc_s_index
    for (size_t i = 0; i < minloc_s_index.size(); ++i) {
        if (minloc_s_index[i] < t.size()) {
            result.S_t[i] = t[minloc_s_index[i]];
        }
    }

    // Check for lead inversion
    // The original code uses minloc.back() and maxloc.back()
    // It should now use result.S_amp and result.R_index derived from the new logic
    if (!result.S_amp.empty() && !result.R_index.empty()) {
        // This condition typically checks if the last S-wave occurs after the last R-peak
        // which might indicate an inverted lead.
        // The previous logic was comparing the *indices* of S_amp and R_index.
        // Assuming minloc_s_index is the source for S_t and result.S_amp
        // and maxloc_r_index is the source for R_t and result.R_amp
        if (minloc_s_index.back() < maxloc_r_index.back()) { // Reverted to using derived indices
            // Swap R and S wave information
            vector<double> temp_R_t = result.R_t;
            vector<double> temp_R_amp = result.R_amp;

            result.R_t = result.S_t;
            result.R_amp = result.S_amp;

            result.S_amp = temp_R_amp;
            result.S_t = temp_R_t; // Assign temp_R_t (which holds original R_t) to S_t
        }
    }


    // Calculate HRV
    result.hrv = diff(result.R_t);

    return result;
}
