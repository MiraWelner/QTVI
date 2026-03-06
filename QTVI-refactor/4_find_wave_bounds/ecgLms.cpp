// ============================================================================
// File: ecgLms.cpp
// ============================================================================
#include "ecgLms.h"
#include "FilterUtils.h"
#include "StatsUtils.h"
#include "diff2.h"
#include "ecglaux.h"

vector<size_t> ecgLms(const vector<double>& ecg, int sampling,
    const vector<double>& b_butter_ecg4mwi,
    const vector<double>& a_butter_ecg4mwi, int dbg) {
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

    // Initiate the algorithm with the first 6-8 seconds
    size_t pt1 = 0;
    size_t pt2 = std::min(static_cast<size_t>(6 * sampling), mwisignal.size());
    size_t pt3 = std::min(static_cast<size_t>(8 * sampling), mwisignal.size());

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

    // --- Polarity detection on first 8 seconds only ---
    vector<double> ecg_seg(origecg.begin() + pt1, origecg.begin() + pt3);
    vector<double> mwi_seg(mwisignal.begin() + pt1, mwisignal.begin() + pt3);

    vector<double> ecg_seg_neg(ecg_seg.size());
    for (size_t i = 0; i < ecg_seg.size(); ++i) {
        ecg_seg_neg[i] = -ecg_seg[i];
    }

    auto result_pos = ecglaux(ecg_seg, mwi_seg, sampling, mwithold, mvimaxval,
        mwiwidthpts, refractpts, mwitholdfract, mwitholdff);
    vector<size_t> rpos = std::get<0>(result_pos);

    auto result_neg = ecglaux(ecg_seg_neg, mwi_seg, sampling, mwithold, mvimaxval,
        mwiwidthpts, refractpts, mwitholdfract, mwitholdff);
    vector<size_t> rneg = std::get<0>(result_neg);

    // --- Determine polarity, then run on FULL signal (matching MATLAB) ---
    // MATLAB: ecg = origecg or ecg = -1*origecg (full length), then calls ecglaux on full signal.
    // The old C++ code incorrectly passed the 8-second segment to the final detection,
    // which is why Algorithm 4 stopped detecting after ~1826 samples.
    vector<double> ecg_full = origecg;  // default: positive polarity

    if (rpos.size() < 2) {
        rpos = { pt1, pt3 - 1 };
    }
    if (rneg.size() < 2) {
        rneg = { pt1, pt3 - 1 };
    }

    if (rpos.size() < 3) {
        // Not enough positive peaks found — signal is inverted
        for (auto& val : ecg_full) val = -val;
    }
    else if (rneg.size() < 3) {
        // Not enough negative peaks — signal is already correct polarity
        // ecg_full = origecg (already set)
    }
    else {
        // MATLAB: abs(median(origecg(rneg))) > abs(median(origecg(rpos)))
        // Note: MATLAB indexes into origecg (non-negated) for BOTH comparisons
        vector<double> rpos_amps, rneg_amps;
        for (size_t idx : rpos) {
            if (idx < ecg_seg.size()) rpos_amps.push_back(origecg[idx]);
        }
        for (size_t idx : rneg) {
            if (idx < ecg_seg.size()) rneg_amps.push_back(origecg[idx]);
        }

        double med_pos = median(rpos_amps);
        double med_neg = median(rneg_amps);

        if (std::abs(med_neg) > std::abs(med_pos)) {
            for (auto& val : ecg_full) val = -val;
        }
    }

    // --- Final detection on the FULL signal with FULL mwisignal ---
    auto final_result = ecglaux(ecg_full, mwisignal, sampling, mwithold, mvimaxval,
        mwiwidthpts, refractpts, mwitholdfract, mwitholdff);

    return std::get<0>(final_result);
}