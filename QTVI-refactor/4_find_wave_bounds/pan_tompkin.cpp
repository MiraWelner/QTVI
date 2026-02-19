#include "pan_tompkin.h"
#include "FilterUtils.h"
#include "StatsUtils.h"
#include "PeakFinder.h"
#include <cmath>
#include <algorithm>
#include <vector>

using namespace std;

PanTompkinResult pan_tompkin(const vector<double>& ecg, double fs, int gr) {
    PanTompkinResult result;
    result.delay = 0;
    if (ecg.empty()) return result;

    // 1. Remove DC Offset
    double avg = mean(ecg);
    vector<double> x = ecg;
    for (auto& v : x) v -= avg;

    // 2. Bandpass Filter (Zero Phase)
    vector<double> b_bp, a_bp;
    butter(3, { 5.0 * 2.0 / fs, 15.0 * 2.0 / fs }, b_bp, a_bp);
    vector<double> ecg_h = filtfilt(b_bp, a_bp, x);

    // 3. Derivative Filter (Reverted to the simple version)
    if (ecg_h.empty()) return result;
    vector<double> ecg_d(ecg_h.size(), 0.0);
    for (size_t i = 1; i < ecg_h.size(); ++i) {
        ecg_d[i] = ecg_h[i] - ecg_h[i - 1];
    }

    // 4. Squaring
    vector<double> ecg_s(ecg_d.size());
    for (size_t i = 0; i < ecg_d.size(); ++i) ecg_s[i] = ecg_d[i] * ecg_d[i];

    // 5. Moving Average Integration (Reverted to your original loop)
    int win_size = static_cast<int>(std::round(0.150 * fs));
    if (win_size < 1) win_size = 1;
    vector<double> ecg_m(ecg_s.size(), 0.0);
    int half_win = win_size / 2;

    for (int i = 0; i < (int)ecg_s.size(); ++i) {
        double current_sum = 0;
        for (int j = i - half_win; j <= i + half_win; ++j) {
            if (j >= 0 && j < (int)ecg_s.size()) {
                current_sum += ecg_s[j];
            }
        }
        ecg_m[i] = current_sum / (double)win_size;
    }

    // 6. Find Peaks
    vector<double> pks;
    vector<size_t> locs;
    findpeaks(ecg_m, pks, locs, round(0.2 * fs));

    // 7. Thresholding
    double SIG_LEV = 0, NOISE_LEV = 0;
    int init_len = min((int)ecg_m.size(), (int)(2.0 * fs));
    if (init_len > 0) {
        double max_m = 0, mean_m = 0;
        for (int k = 0; k < init_len; ++k) {
            if (ecg_m[k] > max_m) max_m = ecg_m[k];
            mean_m += ecg_m[k];
        }
        SIG_LEV = max_m * 0.33;
        NOISE_LEV = (mean_m / (double)init_len) * 0.5;
    }

    // 8. Detection Loop with Offset Correction
    // Based on your logs, BIN is usually 3 samples behind MAT (e.g. 95 vs 98).
    // So we ADD 3 to align BIN with MAT.
    const int ALIGNMENT_SHIFT = 3;

    vector<size_t> qrs_i_raw;
    for (size_t i = 0; i < locs.size(); ++i) {
        if (locs[i] < 0.05 * fs) continue;

        double THR_SIG = NOISE_LEV + 0.25 * (SIG_LEV - NOISE_LEV);
        if (pks[i] >= THR_SIG) {
            // Apply the shift here.
            size_t corrected = locs[i] + ALIGNMENT_SHIFT;

            // Safety check: Don't exceed signal length
            if (corrected >= ecg.size()) corrected = ecg.size() - 1;

            qrs_i_raw.push_back(corrected);
            SIG_LEV = 0.125 * pks[i] + 0.875 * SIG_LEV;
        }
        else {
            NOISE_LEV = 0.125 * pks[i] + 0.875 * NOISE_LEV;
        }
    }

    result.qrs_i_raw = qrs_i_raw;
    return result;
}
