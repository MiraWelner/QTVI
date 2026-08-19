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

#include "bandpass.hpp"
#include "FilterUtils.hpp"
#include "stats_utils.hpp"

using namespace std;

struct RPeakDetectResult {
    vector<size_t> r_peak_index;
    vector<double> hrv;
    vector<double> R_t;
    vector<double> R_amp;
    vector<double> S_t;
    vector<double> S_amp;
};

/**
 * @brief Threshold-independent preprocessing for rpeakdetect.
 *
 * Steps 1-7 of the original rpeakdetect() (detrend, bandpass, diff, square,
 * integration, medfilt1, delay removal, max_h estimation) are identical
 * across calls that only differ in `thresh`. Splitting them out lets a caller
 * (e.g. JoinedRR, which currently invokes rpeakdetect 3 times with thresh
 * 0.2 / 0.1 / 0.4) do this work once and reuse it.
 *
 * The returned struct carries everything rpeakdetect_apply() needs to finish
 * detection at any threshold.
 */
struct RPeakDetectPrep {
    vector<double> bpf;     // bandpass-filtered, detrended signal
    vector<double> mdfint;  // integration output after delay removal
    double max_h = 0.0;     // max of mdfint over the middle 50% of len
    size_t len = 0;         // original sqr length (used for threshold-region indexing)
    double samp_freq = 0.0;
    bool ok = false;
};

inline RPeakDetectPrep rpeakdetect_prep(const vector<double>& data,
    double samp_freq = 256.0,
    int /*testmode*/ = 0,
    std::string /*fileID*/ = "")
{
    RPeakDetectPrep prep;
    prep.samp_freq = samp_freq;
    if (data.size() < 10) return prep;

    std::vector<double> x = data;

    // 1. Detrend
    double mu = mean(x);
    for (auto& val : x) val -= mu;

    // 2. Bandpass filtering
    prep.bpf = bandpass_filtfilt(2, 0.05, 100.0, 1000.0, x);

    // 3. Differentiation
    std::vector<double> dff = diff(prep.bpf);

    // 4. Squaring
    std::vector<double> sqr(dff.size());
    for (size_t i = 0; i < dff.size(); ++i) sqr[i] = dff[i] * dff[i];
    prep.len = sqr.size();

    // 5. Integration (moving sum + median filter)
    int win_size = (samp_freq >= 256.0) ? (int)std::round(7.0 * samp_freq / 256.0) : 7;
    std::vector<double> d_kernel(win_size, 1.0);
    std::vector<double> filtered_sqr = filter(d_kernel, { 1.0 }, sqr);
    prep.mdfint = medfilt1(filtered_sqr, 10);

    // 6. Remove filter delay
    int delay = (int)std::ceil((double)win_size / 2.0);
    if (delay > 1 && (size_t)delay <= prep.mdfint.size()) {
        prep.mdfint.erase(prep.mdfint.begin(), prep.mdfint.begin() + (delay - 1));
    }

    // 7. Threshold reference from middle 50% of signal
    int start_search = std::max(0, std::min((int)std::round((double)prep.len / 4.0) - 1,
        (int)prep.mdfint.size() - 1));
    int end_search = std::max(start_search,
        std::min((int)std::round(3.0 * (double)prep.len / 4.0) - 1,
            (int)prep.mdfint.size() - 1));
    prep.max_h = 0.0;
    for (int i = start_search; i <= end_search; ++i) {
        if (std::isfinite(prep.mdfint[i]) && prep.mdfint[i] > prep.max_h) {
            prep.max_h = prep.mdfint[i];
        }
    }
    prep.ok = true;
    return prep;
}

/**
 * @brief Threshold-dependent portion: steps 8-12 of the original rpeakdetect.
 *        Uses the prep struct produced by rpeakdetect_prep().
 *
 * @param[in] inverted  True if this channel's lead is inverted (from the
 *                      GUI's "Inverted Lead?" checkbox, persisted through
 *                      the annealed .bin). When true, the true R-peak
 *                      within each candidate region is the local min
 *                      instead of the local max.
 */
inline RPeakDetectResult rpeakdetect_apply(const RPeakDetectPrep& prep, double thresh, bool inverted = false) {
    RPeakDetectResult result;
    if (!prep.ok) return result;

    const auto& mdfint = prep.mdfint;
    const auto& bpf = prep.bpf;
    const double samp_freq = prep.samp_freq;

    // 8. Identify candidate regions above threshold
    std::vector<int> poss_reg(mdfint.size(), 0);
    double limit = thresh * prep.max_h;
    for (size_t i = 0; i < mdfint.size(); ++i) {
        if (!std::isnan(mdfint[i]) && mdfint[i] > limit) poss_reg[i] = 1;
    }

    // 9. Find region boundaries (rising/falling edges)
    std::vector<int> left, right;
    if (!poss_reg.empty()) {
        if (poss_reg[0] == 1) left.push_back(1);
        for (size_t i = 1; i < poss_reg.size(); ++i) {
            if (poss_reg[i] == 1 && poss_reg[i - 1] == 0) left.push_back(static_cast<int>(i + 1));
            if (poss_reg[i] == 0 && poss_reg[i - 1] == 1) right.push_back(static_cast<int>(i));
        }
        if (poss_reg.back() == 1) right.push_back((int)poss_reg.size());
    }

    // 10. Find local max within each region (R-peaks)
    size_t num_segs = std::min(left.size(), right.size());
    std::vector<size_t> maxloc;
    std::vector<double> maxval;
    maxloc.reserve(num_segs);
    maxval.reserve(num_segs);

    for (size_t i = 0; i < num_segs; ++i) {
        int l = left[i] - 1;
        int r = right[i] - 1;
        double cur_best = inverted ? 1e30 : -1e30;
        int cur_best_i = l;
        for (int k = l; k <= r; ++k) {
            if (k >= 0 && k < (int)bpf.size() && !std::isnan(bpf[k])) {
                bool better = inverted ? (bpf[k] < cur_best) : (bpf[k] > cur_best);
                if (better) { cur_best = bpf[k]; cur_best_i = k; }
            }
        }
        maxval.push_back(cur_best);
        maxloc.push_back((size_t)cur_best_i);
    }

    // 11. Assign results
    result.r_peak_index = maxloc;
    result.R_amp = maxval;
    auto get_time = [&](size_t idx) { return (double)(idx + 1) / samp_freq; };
    result.R_t.reserve(maxloc.size());
    for (auto loc : maxloc) result.R_t.push_back(get_time(loc));

    // 12. HRV
    if (result.R_t.size() > 1) {
        result.hrv.reserve(result.R_t.size() - 1);
        for (size_t i = 0; i + 1 < result.R_t.size(); ++i)
            result.hrv.push_back(result.R_t[i + 1] - result.R_t[i]);
    }
    return result;
}

/**
 * @brief Backward-compatible wrapper: identical signature and behavior to the
 *        original rpeakdetect(), plus an inverted-lead flag. Internally
 *        calls prep + apply.
 */
inline RPeakDetectResult rpeakdetect(const vector<double>& data,
    double samp_freq = 256.0,
    double thresh = 0.2,
    int testmode = 0,
    std::string fileID = "",
    bool inverted = false)
{
    auto prep = rpeakdetect_prep(data, samp_freq, testmode, fileID);
    if (!prep.ok) return RPeakDetectResult{};
    return rpeakdetect_apply(prep, thresh, inverted);
}