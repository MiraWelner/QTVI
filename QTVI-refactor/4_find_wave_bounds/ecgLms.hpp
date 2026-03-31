// ============================================================================
// File: ecgLms.hpp
// ECG LMS-based R-wave detection
// ============================================================================
#pragma once

#include "SignalProcessingTypes.h"
#include "FilterUtils.hpp"
#include "StatsUtils.hpp"
#include "diff2.hpp"
#include "ecglaux.hpp"
#include <iomanip>

inline vector<size_t> ecgLms(const vector<double>& ecg, int sampling,
    const vector<double>& b_butter_ecg4mwi,
    const vector<double>& a_butter_ecg4mwi, int dbg = 0, const std::string& fileID = "") {
    double mwitholdfract = 0.25;
    double mwitholdff = 0.80;

    // Baseline wander removal
    double meanVal = mean(ecg);
    vector<double> origecg(ecg.size());
    for (size_t i = 0; i < ecg.size(); ++i) {
        origecg[i] = ecg[i] - meanVal;
    }

    double sl = sampling / 1000.0;

    int mwiwidthmsec = 175;
    int mwiwidthpts = static_cast<int>(std::round(mwiwidthmsec * sl));
    int refractmsec = 250;
    int refractpts = static_cast<int>(std::round(refractmsec * sl));

    // Filter the ECG.
    // When called from JoinedRR with scalar b={5}, a={12}, MATLAB's filtfilt(5,12,x)
    // applies a zero-phase constant gain: y = x * (b/a)^2 = x * (5/12)^2.
    // A general C++ filtfilt may not handle single-element coefficient vectors correctly,
    // so we detect that case and apply the scaling directly.
    vector<double> filtecg;
    if (b_butter_ecg4mwi.size() == 1 && a_butter_ecg4mwi.size() == 1) {
        double gain = b_butter_ecg4mwi[0] / a_butter_ecg4mwi[0];
        double gain_sq = gain * gain; // forward + backward pass
        filtecg.resize(origecg.size());
        for (size_t i = 0; i < origecg.size(); ++i) {
            filtecg[i] = origecg[i] * gain_sq;
        }
    }
    else {
        filtecg = filtfilt(b_butter_ecg4mwi, a_butter_ecg4mwi, origecg);
    }
    vector<double> difffiltecg = diff2(filtecg);

    vector<double> sqdifffiltecg(difffiltecg.size());
    for (size_t i = 0; i < difffiltecg.size(); ++i) {
        sqdifffiltecg[i] = difffiltecg[i] * difffiltecg[i];
    }

    size_t ll = sqdifffiltecg.size();

    // Create Moving-Window-Integration
    vector<double> mwisignal(ll, 0.0);
    for (int cnt = 0; cnt < mwiwidthpts && cnt < static_cast<int>(ll); ++cnt) {
        for (int j = 0; j <= cnt; ++j) {
            mwisignal[cnt] += sqdifffiltecg[j];
        }
    }

    int l = mwiwidthpts;
    vector<double> wholesum(ll);
    wholesum[0] = sqdifffiltecg[0];
    for (size_t i = 1; i < ll; ++i) {
        wholesum[i] = wholesum[i - 1] + sqdifffiltecg[i];
    }

    double maxWhole = *std::max_element(wholesum.begin(), wholesum.end());
    if (maxWhole > std::numeric_limits<double>::max() / 100.0) {
        throw std::runtime_error("wholesum comes close to exceeding max allowed value");
    }

    for (size_t i = l; i < ll; ++i) {
        mwisignal[i] = wholesum[i] - wholesum[i - l];
    }

    // The beginning of mwisignal always starts near zero
    for (int i = 0; i < mwiwidthpts && i < static_cast<int>(ll); ++i) {
        if (mwiwidthpts < static_cast<int>(ll)) {
            mwisignal[i] = mwisignal[mwiwidthpts];
        }
    }

    // Initiate the algorithm with the first 6 seconds
    size_t pt1 = 0;
    size_t pt2 = std::min(static_cast<size_t>(6 * sampling), mwisignal.size());

    // Sort and get mvimaxval from first 6 seconds of MWI
    vector<double> sorted_segment(mwisignal.begin() + pt1, mwisignal.begin() + pt2);
    std::sort(sorted_segment.begin(), sorted_segment.end());
    size_t lx = sorted_segment.size();

    double mvimaxval = 0.0;
    size_t start_idx = static_cast<size_t>(std::round(0.90 * lx)) - 1;
    size_t end_idx = static_cast<size_t>(std::round(0.95 * lx)) - 1;
    double sum_range = 0.0;
    int count = 0;
    for (size_t i = start_idx; i <= end_idx && i < sorted_segment.size(); ++i) {
        sum_range += sorted_segment[i];
        count++;
    }
    mvimaxval = (count > 0) ? sum_range / count : 0.0;

    double mwithold = mwitholdfract * mvimaxval;

    // Run detection on the full signal directly (polarity is handled upstream)
    auto final_result = ecglaux(origecg, mwisignal, sampling, mwithold, mvimaxval,
        mwiwidthpts, refractpts, mwitholdfract, mwitholdff);

    return std::get<0>(final_result);
}