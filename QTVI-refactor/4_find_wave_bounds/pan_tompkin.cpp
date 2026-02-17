// ============================================================================
// File: pan_tompkin.cpp
// FIXED VERSION: Prevents Threshold Deadlock and aligns with MATLAB logic
// ============================================================================
#include "pan_tompkin.h"
#include "FilterUtils.h"
#include "StatsUtils.h"
#include "PeakFinder.h"
#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;

PanTompkinResult pan_tompkin(const vector<double>& ecg, double fs, int gr) {
    PanTompkinResult result;
    result.delay = 0;
    if (ecg.empty()) return result;

    // 1. Pre-processing: Remove Mean
    double avg = mean(ecg);
    vector<double> x = ecg;
    for (auto& v : x) v -= avg;

    // 2. Bandpass Filter (5-15 Hz)
    vector<double> b_bp, a_bp;
    butter(3, { 5.0 * 2.0 / fs, 15.0 * 2.0 / fs }, b_bp, a_bp);
    vector<double> ecg_h = filtfilt(b_bp, a_bp, x);

    // Normalize ecg_h
    double max_h_val = 0.0;
    for (auto v : ecg_h) max_h_val = max(max_h_val, abs(v));
    if (max_h_val > 0) { for (auto& v : ecg_h) v /= max_h_val; }

    // 3. Derivative Filter (Delay = 2 samples)
    vector<double> b_d = { 1, 2, 0, -2, -1 };
    for (auto& v : b_d) v *= (fs / 8.0);
    vector<double> ecg_d = filter(b_d, { 1.0 }, ecg_h);

    // Normalize ecg_d
    double max_d = 0.0;
    for (auto v : ecg_d) max_d = max(max_d, abs(v));
    if (max_d > 0) { for (auto& v : ecg_d) v /= max_d; }

    // 4. Squaring
    vector<double> ecg_s(ecg_d.size());
    for (size_t i = 0; i < ecg_d.size(); ++i) ecg_s[i] = ecg_d[i] * ecg_d[i];

    // 5. Moving Average Integration (MATLAB 'same' emulation)
    int win = (int)round(0.150 * fs);
    vector<double> ecg_m(ecg_s.size(), 0.0);
    int half_win = win / 2;

    for (int i = 0; i < (int)ecg_s.size(); ++i) {
        double current_sum = 0;
        int count = 0;
        for (int j = i - half_win; j <= i + half_win; ++j) {
            if (j >= 0 && j < (int)ecg_s.size()) {
                current_sum += ecg_s[j];
                count++;
            }
        }
        ecg_m[i] = current_sum / win;
    }

    // 6. Initial Threshold Estimation (CRITICAL FIX)
    // We look at the first 2 seconds to estimate signal/noise levels
    int init_samples = min((int)ecg_m.size(), (int)(2.0 * fs));
    double max_init = 0.0;
    for (int i = 0; i < init_samples; ++i) max_init = max(max_init, ecg_m[i]);

    double SIG_LEV = max_init * 0.5;
    double NOISE_LEV = mean(vector<double>(ecg_m.begin(), ecg_m.begin() + init_samples)) * 0.5;
    double THR_SIG = NOISE_LEV + 0.25 * (SIG_LEV - NOISE_LEV);

    // 7. Find Peaks in Integrated Signal
    vector<double> pks;
    vector<size_t> locs;
    findpeaks(ecg_m, pks, locs, round(0.2 * fs));

    int total_delay = 2 + (win / 2);
    vector<size_t> qrs_i_raw;
    vector<double> qrs_amp_raw;
    int last_qrs_idx = -1;

    // 8. Detection Loop
    for (size_t i = 0; i < locs.size(); ++i) {
        size_t loc = locs[i];
        double peak_val = pks[i];

        // Search Back Logic: If no heart beat found for a while, lower threshold
        if (last_qrs_idx != -1 && (int(loc) - last_qrs_idx) > (1.5 * fs)) {
            THR_SIG *= 0.5; // Lower the bar to find missed peaks
        }

        // Check against adaptive threshold
        if (peak_val > THR_SIG) {
            // Find R-peak in Bandpass signal (ecg_h) within 150ms of integrated peak
            int compensated_loc = (int)loc - total_delay;
            int radius = (int)round(0.15 * fs);
            int start = max(0, compensated_loc - radius);
            int end = min((int)ecg_h.size() - 1, compensated_loc + radius);

            auto h_res = max_element_index(ecg_h, start, end + 1);

            qrs_i_raw.push_back(start + h_res.second);
            qrs_amp_raw.push_back(h_res.first);

            // Update Signal Level
            SIG_LEV = 0.125 * peak_val + 0.875 * SIG_LEV;
            last_qrs_idx = (int)loc;
        }
        else {
            // Update Noise Level
            NOISE_LEV = 0.125 * peak_val + 0.875 * NOISE_LEV;
        }

        // Update Threshold for next iteration
        THR_SIG = NOISE_LEV + 0.25 * (SIG_LEV - NOISE_LEV);
    }

    result.qrs_i_raw = qrs_i_raw;
    result.qrs_amp_raw = qrs_amp_raw;
    result.delay = total_delay;
    return result;
}
