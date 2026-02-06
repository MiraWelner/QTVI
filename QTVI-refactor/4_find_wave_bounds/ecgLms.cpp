
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
    double mwitholdfract = 0.25;  // normally 0.30, low # = very tolerant, easy detection
    double mwitholdff = 0.80;     // usually .80

    // Baseline wander removal
    double meanVal = mean(ecg);
    vector<double> ecg_centered = ecg;
    for (auto& val : ecg_centered) {
        val -= meanVal;
    }

    vector<double> origecg = ecg_centered;
    double sl = sampling / 1000.0;

    int mwiwidthmsec = 175;  // usually 175
    int mwiwidthpts = static_cast<int>(std::round(mwiwidthmsec * sl));
    int refractmsec = 250;   // usually 250
    int refractpts = static_cast<int>(std::round(refractmsec * sl));

    // filter, differentiate and square ecg
    vector<double> filtecg = filtfilt(b_butter_ecg4mwi, a_butter_ecg4mwi, origecg);
    vector<double> difffiltecg = diff2(filtecg);

    vector<double> sqdifffiltecg(difffiltecg.size());
    for (size_t i = 0; i < difffiltecg.size(); ++i) {
        sqdifffiltecg[i] = difffiltecg[i] * difffiltecg[i];
    }

    size_t ll = sqdifffiltecg.size();

    // create Moving-Window-Integration
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

    // Check for overflow
    double maxWhole = *std::max_element(wholesum.begin(), wholesum.end());
    if (maxWhole > std::numeric_limits<double>::max() / 100.0) {
        throw std::runtime_error("wholesum comes close to exceeding max allowed value");
    }

    for (size_t i = l; i < ll; ++i) {
        mwisignal[i] = wholesum[i] - wholesum[i - l];
    }

    // The beginning of mwisignal always starts near zero, this is sometimes bad
    for (int i = 0; i < mwiwidthpts && i < static_cast<int>(ll); ++i) {
        if (mwiwidthpts < static_cast<int>(ll)) {
            mwisignal[i] = mwisignal[mwiwidthpts];
        }
    }

    // initiate the algorithm with the first 6-8 seconds
    size_t pt1 = 0;
    size_t pt2 = std::min(static_cast<size_t>(6 * sampling), mwisignal.size());
    size_t pt3 = std::min(static_cast<size_t>(8 * sampling), mwisignal.size());

    // Sort and get mvimaxval
    vector<double> sorted_segment(mwisignal.begin() + pt1, mwisignal.begin() + pt2);
    std::sort(sorted_segment.begin(), sorted_segment.end());
    size_t lx = sorted_segment.size();

    double mvimaxval = 0.0;
    size_t start_idx = static_cast<size_t>(std::round(0.90 * lx));
    size_t end_idx = static_cast<size_t>(std::round(0.95 * lx));
    double sum_range = 0.0;
    size_t count_range = 0;

    for (size_t i = start_idx; i < end_idx && i < sorted_segment.size(); ++i) {
        sum_range += sorted_segment[i];
        count_range++;
    }
    mvimaxval = count_range > 0 ? sum_range / count_range : 0.0;

    double mwithold = mwitholdfract * mvimaxval;

    // Run detection on positive and negative ECG
    vector<double> ecg_segment_pos(origecg.begin() + pt1, origecg.begin() + pt3);
    vector<double> mwi_segment(mwisignal.begin() + pt1, mwisignal.begin() + pt3);

    auto result_pos = ecglaux(ecg_segment_pos, mwi_segment, sampling, mwithold, mvimaxval,
        mwiwidthpts, refractpts, mwitholdfract, mwitholdff);
    vector<size_t> rpos = std::get<0>(result_pos);

    vector<double> ecg_segment_neg = ecg_segment_pos;
    for (auto& val : ecg_segment_neg) {
        val = -val;
    }

    auto result_neg = ecglaux(ecg_segment_neg, mwi_segment, sampling, mwithold, mvimaxval,
        mwiwidthpts, refractpts, mwitholdfract, mwitholdff);
    vector<size_t> rneg = std::get<0>(result_neg);

    // Determine polarity
    vector<double> ecg_to_use;
    if (rpos.size() < 2) {
        rpos = { pt1, pt3 - 1 };
    }
    if (rneg.size() < 2) {
        rneg = { pt1, pt3 - 1 };
    }

    if (rpos.size() < 3) {
        ecg_to_use = ecg_segment_neg;
    }
    else if (rneg.size() < 3) {
        ecg_to_use = ecg_segment_pos;
    }
    else {
        // Calculate median
        double sum_pos = 0.0, sum_neg = 0.0;
        for (size_t idx : rpos) {
            if (idx < ecg_segment_pos.size()) {
                sum_pos += ecg_segment_pos[idx];
            }
        }
        for (size_t idx : rneg) {
            if (idx < ecg_segment_neg.size()) {
                sum_neg += ecg_segment_neg[idx];
            }
        }

        double median_pos = rpos.empty() ? 0.0 : sum_pos / rpos.size();
        double median_neg = rneg.empty() ? 0.0 : sum_neg / rneg.size();

        if (std::abs(median_neg) > std::abs(median_pos)) {
            ecg_to_use = ecg_segment_neg;
        }
        else {
            ecg_to_use = ecg_segment_pos;
        }
    }

    // Run full detection
    auto final_result = ecglaux(ecg_to_use.empty() ? origecg : ecg_to_use,
        mwisignal, sampling, mwithold, mvimaxval,
        mwiwidthpts, refractpts, mwitholdfract, mwitholdff);

    return std::get<0>(final_result);
}
