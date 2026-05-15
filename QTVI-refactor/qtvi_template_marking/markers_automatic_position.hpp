#pragma once
//
// This module automatically detects the locations of the five markers:
//
// Beginning of Q wave
// Beginning of T wave
// End of T wave
// 
// Onset of PPG
// Peak of PPG
// 
// These are only approximations, and will need to be adjusted by the user (that is the whole point of this GUI)
//

#include <vector>
#include <iostream>
#include <fstream>
namespace ecg_markers {

    inline std::vector<double> first_derivative(const std::vector<double>& v) {
        const int N = static_cast<int>(v.size());
        std::vector<double> d1(N, 0.0);
        if (N < 2) return d1;

        for (int i = 1; i < N - 1; ++i) {
            d1[i] = 0.5 * (v[i + 1] - v[i - 1]);
        }
        d1[0] = v[1] - v[0];
        d1[N - 1] = v[N - 1] - v[N - 2];
        return d1;
    }

    inline std::pair<int, bool> r_peak(const std::vector<double>& v) {
        /*
            Get the R peak, and if its negative (which means lead inverted)
            you only select in the first half because sometimes the templates contain
            a second R peak - if this gets fixed upstream you can remove it
        */
        const int half = static_cast<int>(v.size()) / 2;
        if (half < 1) return { -1, true };
        // Use deviation from the mean — DC offset shouldn't decide polarity.
        double mean = 0.0;
        for (int i = 0; i < half; ++i) mean += v[i];
        mean /= half;

        auto it = std::max_element(v.begin(), v.begin() + half,
            [mean](double a, double b) {
                return std::abs(a - mean) < std::abs(b - mean);
            });

        int idx = static_cast<int>(it - v.begin());
        bool isPositive = (*it >= mean);   // deviates UP from baseline?
        return { idx, isPositive };
    }


    // ---- ECG Marker detectors -----------------------------------

    inline int detect_q_begin(const std::vector<double>& ecg_signal) {

        auto [r_idx, is_positive] = r_peak(ecg_signal);

        // If R-peak is negative, flip the signal so the standard upright-R
        // algorithm works unchanged. We work on a local copy either way.
        std::vector<double> upright_signal = ecg_signal;
        if (!is_positive) {
            for (auto& x : upright_signal) x = -x;
        }

        // Bound the search to just after the R peak
        int search_lim = std::min(r_idx-1, static_cast<int>(ecg_signal.size()));

        const auto d1 = first_derivative(upright_signal);

        // Step 1: Q-trough (lowest beat value up to and including r peak)
        int qTrough = 0;
        double qVal = upright_signal[0];
        for (int i = 1; i <= search_lim; ++i) {
            if (upright_signal[i] < qVal) { qVal = upright_signal[i]; qTrough = i; }
        }
        if (qTrough < 1) return -1;


        // Step 2: walk back; first sample whose slope is >= 0 is the
        // end of the flat PR-segment and the start of the Q-descent.
        for (int i = qTrough - 1; i >= 0; --i) {
            if (d1[i] >= 0.0) return i;
        }
        return -1;
    }


    inline int detect_t_begin(const std::vector<double>& ecg_signal) {

        auto [r_idx, is_positive] = r_peak(ecg_signal);
        if (r_idx < 0) return -1;

        // Flip to upright so the threshold/argmax logic works the same way
        // for inverted-R leads (T usually has the same polarity as R).
        std::vector<double> upright = ecg_signal;
        if (!is_positive) {
            for (auto& x : upright) x = -x;
        }

        const int N = static_cast<int>(upright.size());

        // ST-baseline: median of a flat window after the S-trough recovery.
        const int st_lo = r_idx + 50;
        const int st_hi = std::min(r_idx + 100, N);
        if (st_hi - st_lo < 5) return -1;
        std::vector<double> st_window(upright.begin() + st_lo,
            upright.begin() + st_hi);
        std::nth_element(st_window.begin(),
            st_window.begin() + st_window.size() / 2,
            st_window.end());
        const double st_baseline = st_window[st_window.size() / 2];

        // T-peak search range: avoid QRS-and-recovery on the front, and the
        // next beat creeping in on the back.
        const int t_lo = r_idx + 100;
        const int t_hi = std::min(static_cast<int>(0.55 * N), N);
        if (t_hi - t_lo < 10) return -1;

        int t_peak = t_lo;
        double t_peak_val = upright[t_lo];
        for (int i = t_lo + 1; i < t_hi; ++i) {
            if (upright[i] > t_peak_val) { t_peak_val = upright[i]; t_peak = i; }
        }

        const double t_height = t_peak_val - st_baseline;
        if (t_height <= 0.0) return -1;  // no upward T detected

        // Walk back from T-peak; first sample below 10% of T-height is T-begin.
        const double thresh = st_baseline + 0.10 * t_height;
        for (int i = t_peak - 1; i > r_idx; --i) {
            if (upright[i] < thresh) return i;
        }
        return -1;
    }
    inline int detect_t_end(const std::vector<double>& ecg_signal) {
        auto [r_idx, is_positive] = r_peak(ecg_signal);
        if (r_idx < 0) return -1;

        std::vector<double> upright = ecg_signal;
        if (!is_positive) for (auto& x : upright) x = -x;

        const int N = static_cast<int>(upright.size());
        const int st_lo = r_idx + 50;
        const int st_hi = std::min(r_idx + 100, N);
        if (st_hi - st_lo < 5) return -1;

        std::vector<double> w(upright.begin() + st_lo, upright.begin() + st_hi);
        std::nth_element(w.begin(), w.begin() + w.size() / 2, w.end());
        const double baseline = w[w.size() / 2];

        double st_noise = 0.0;
        for (int i = st_lo; i < st_hi; ++i) {
            st_noise = std::max(st_noise, std::abs(upright[i] - baseline));
        }
        const double tol = 2.0 * st_noise;

        // Find T-begin first (first sample above baseline+tol).
        int t_begin = -1;
        for (int i = st_hi; i < N; ++i) {
            if (upright[i] > baseline + tol) { t_begin = i; break; }
        }
        if (t_begin < 0) return -1;

        // T-end: first sample after T-begin that drops back below baseline+tol.
        for (int i = t_begin + 1; i < N; ++i) {
            if (upright[i] <= baseline + tol) return i;
        }
        return detect_t_begin(ecg_signal)*1.5; //if you can't find it, just do 20% added to the begin t
    }


    // ---- PPG detectors --------------------------------------------------
    // `pulse` is one PPG template waveform (one cardiac cycle, onset to
    // end). No R-peak hint -- PPG features are found directly from the
    // pulse shape.

    inline int detect_ppg_onset(const std::vector<double>& pulse) {
        return 1;
    }

    inline int detect_ppg_peak(const std::vector<double>& pulse) {
        auto it = std::max_element(pulse.begin(), pulse.end());
        return static_cast<int>(it - pulse.begin());
    }

} // namespace ecg_markers