#pragma once
//
// Automatic detection of marker positions on ECG and PPG templates.
//
// All detectors return a non-negative position. When the heuristic
// detection fails (no clear P-wave, no clear T-end, no dicrotic notch,
// etc.) the detector falls back to a sensible default position rather
// than returning -1. This ensures every marker is visible and draggable
// in the GUI even on poor-quality templates.
//
// Markers detected here:
//   ECG (per channel):  P-onset, Q-begin, T-begin, T-end
//   PPG:                Onset, Peak, Dicrotic notch, 50% upstroke, End
//

#include <vector>
#include <algorithm>
#include <cmath>

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
        if (half < 1) return { 0, true };
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

    // P-wave onset. Falls back to a position before the R-peak (or
    // near the start of the template) when the heuristic fails.
    inline int detect_p_begin(const std::vector<double>& ecg_signal) {
        const int N = static_cast<int>(ecg_signal.size());
        auto [r_idx, is_positive] = r_peak(ecg_signal);
        if (r_idx < 0) return std::min(5, std::max(0, N - 1));

        std::vector<double> upright = ecg_signal;
        if (!is_positive) for (auto& x : upright) x = -x;

        const int q_idx_approx = std::max(0, r_idx - 40);
        if (q_idx_approx < 20) return std::max(0, r_idx / 4);

        // PR-baseline window before Q, avoiding the P-wave itself.
        const int baseline_lo = std::max(0, q_idx_approx - 200);
        const int baseline_hi = std::max(0, q_idx_approx - 150);
        if (baseline_hi - baseline_lo < 5) return std::max(0, q_idx_approx - 60);

        std::vector<double> w(upright.begin() + baseline_lo,
            upright.begin() + baseline_hi);
        std::nth_element(w.begin(), w.begin() + w.size() / 2, w.end());
        const double baseline = w[w.size() / 2];

        double noise = 0.0;
        for (int i = baseline_lo; i < baseline_hi; ++i) {
            noise = std::max(noise, std::abs(upright[i] - baseline));
        }
        const double thresh = baseline + 2.0 * noise;

        const int p_lo = baseline_hi;
        const int p_hi = q_idx_approx;
        if (p_hi - p_lo < 5) return std::max(0, q_idx_approx - 60);

        int p_peak = p_lo;
        double p_peak_val = upright[p_lo];
        for (int i = p_lo + 1; i < p_hi; ++i) {
            if (upright[i] > p_peak_val) { p_peak_val = upright[i]; p_peak = i; }
        }
        if (p_peak_val <= thresh) return std::max(0, q_idx_approx - 60);

        for (int i = p_peak - 1; i >= baseline_lo; --i) {
            if (upright[i] <= thresh) return i;
        }
        return p_peak;
    }

    inline int detect_q_begin(const std::vector<double>& ecg_signal) {
        const int N = static_cast<int>(ecg_signal.size());
        auto [r_idx, is_positive] = r_peak(ecg_signal);

        std::vector<double> upright_signal = ecg_signal;
        if (!is_positive) for (auto& x : upright_signal) x = -x;

        int search_lim = std::min(r_idx - 1, N);
        const auto d1 = first_derivative(upright_signal);

        // Q-trough: lowest beat value up to and including the R peak.
        int qTrough = 0;
        double qVal = upright_signal[0];
        for (int i = 1; i <= search_lim; ++i) {
            if (upright_signal[i] < qVal) { qVal = upright_signal[i]; qTrough = i; }
        }
        if (qTrough < 1) return std::max(0, r_idx - 20);

        // Walk back; first sample whose slope is >= 0 is the start of
        // the Q-descent.
        for (int i = qTrough - 1; i >= 0; --i) {
            if (d1[i] >= 0.0) return i;
        }
        return qTrough;
    }

    inline int detect_t_begin(const std::vector<double>& ecg_signal) {
        const int N = static_cast<int>(ecg_signal.size());
        auto [r_idx, is_positive] = r_peak(ecg_signal);
        if (r_idx < 0) return std::max(0, N / 3);

        std::vector<double> upright = ecg_signal;
        if (!is_positive) for (auto& x : upright) x = -x;

        const int st_lo = r_idx + 50;
        const int st_hi = std::min(r_idx + 100, N);
        if (st_hi - st_lo < 5)
            return std::min(N - 1, r_idx + (N - r_idx) / 3);

        std::vector<double> st_window(upright.begin() + st_lo,
            upright.begin() + st_hi);
        std::nth_element(st_window.begin(),
            st_window.begin() + st_window.size() / 2,
            st_window.end());
        const double st_baseline = st_window[st_window.size() / 2];

        const int t_lo = r_idx + 100;
        const int t_hi = std::min(static_cast<int>(0.55 * N), N);
        if (t_hi - t_lo < 10) return std::min(N - 1, r_idx + 150);

        int t_peak = t_lo;
        double t_peak_val = upright[t_lo];
        for (int i = t_lo + 1; i < t_hi; ++i) {
            if (upright[i] > t_peak_val) { t_peak_val = upright[i]; t_peak = i; }
        }

        const double t_height = t_peak_val - st_baseline;
        if (t_height <= 0.0) return t_lo;

        const double thresh = st_baseline + 0.10 * t_height;
        for (int i = t_peak - 1; i > r_idx; --i) {
            if (upright[i] < thresh) return i;
        }
        return t_peak;
    }

    inline int detect_t_end(const std::vector<double>& ecg_signal) {
        const int N = static_cast<int>(ecg_signal.size());
        auto [r_idx, is_positive] = r_peak(ecg_signal);
        if (r_idx < 0) return std::max(0, (2 * N) / 3);

        std::vector<double> upright = ecg_signal;
        if (!is_positive) for (auto& x : upright) x = -x;

        const int st_lo = r_idx + 50;
        const int st_hi = std::min(r_idx + 100, N);
        if (st_hi - st_lo < 5) return std::max(0, N - 5);

        std::vector<double> w(upright.begin() + st_lo, upright.begin() + st_hi);
        std::nth_element(w.begin(), w.begin() + w.size() / 2, w.end());
        const double baseline = w[w.size() / 2];

        double st_noise = 0.0;
        for (int i = st_lo; i < st_hi; ++i) {
            st_noise = std::max(st_noise, std::abs(upright[i] - baseline));
        }
        const double tol = 2.0 * st_noise;

        int t_begin = -1;
        for (int i = st_hi; i < N; ++i) {
            if (upright[i] > baseline + tol) { t_begin = i; break; }
        }
        if (t_begin < 0) return std::max(0, N - 5);

        for (int i = t_begin + 1; i < N; ++i) {
            if (upright[i] <= baseline + tol) return i;
        }
        const int fallback = static_cast<int>(t_begin * 1.5);
        return std::clamp(fallback, t_begin + 1, N - 1);
    }


    // ---- PPG detectors --------------------------------------------------
    // `pulse` is one PPG template waveform (one cardiac cycle, onset to
    // end). No R-peak hint -- PPG features are found directly from the
    // pulse shape.

    inline int detect_ppg_onset(const std::vector<double>& pulse) {
        const int N = static_cast<int>(pulse.size());
        if (N < 2) return 0;
        const int half = std::max(2, N / 2);

        int idx = 0;
        double v = pulse[0];
        for (int i = 1; i < half; ++i) {
            if (pulse[i] < v) { v = pulse[i]; idx = i; }
        }
        // If the min is at sample 0, prefer a slightly inside position.
        if (idx == 0) return std::min(5, N - 1);
        return idx;
    }

    inline int detect_ppg_peak(const std::vector<double>& pulse) {
        if (pulse.empty()) return 0;
        auto it = std::max_element(pulse.begin(), pulse.end());
        return static_cast<int>(it - pulse.begin());
    }

    inline int detect_ppg_end(const std::vector<double>& pulse) {
        const int N = static_cast<int>(pulse.size());
        if (N < 4) return std::max(0, N - 1);
        const int peak = detect_ppg_peak(pulse);
        if (peak < 0 || peak >= N - 2) return std::max(0, N - 1);

        int end = peak + 1;
        double endVal = pulse[end];
        for (int i = peak + 2; i < N; ++i) {
            if (pulse[i] < endVal) { endVal = pulse[i]; end = i; }
        }
        return end;
    }

    // Dicrotic notch: local minimum on the descending limb between the
    // systolic peak and the end of the pulse.
    inline int detect_ppg_dicrotic(const std::vector<double>& pulse) {
        const int N = static_cast<int>(pulse.size());
        const int peak = detect_ppg_peak(pulse);
        const int end = detect_ppg_end(pulse);
        if (peak < 0 || end < 0 || end - peak < 10)
            return std::clamp((peak + end) / 2, 0, N - 1);

        const int margin = std::max(2, (end - peak) / 10);
        const int lo = peak + margin;
        const int hi = end - 1;
        if (hi - lo < 3)
            return std::clamp(peak + (end - peak) / 3, 0, N - 1);

        int best = -1;
        double bestVal = 1e300;
        for (int i = lo + 1; i < hi; ++i) {
            if (pulse[i] <= pulse[i - 1] && pulse[i] <= pulse[i + 1]) {
                if (pulse[i] < bestVal) { bestVal = pulse[i]; best = i; }
            }
        }
        if (best < 0) return std::clamp(peak + (end - peak) / 3, 0, N - 1);
        return best;
    }

    // 50% point on the upstroke: between PPG onset (foot) and peak.
    inline int detect_ppg_50(const std::vector<double>& pulse) {
        const int N = static_cast<int>(pulse.size());
        const int onset = detect_ppg_onset(pulse);
        const int peak = detect_ppg_peak(pulse);
        if (onset < 0 || peak < 0 || peak <= onset)
            return std::clamp((onset + peak) / 2, 0, N - 1);

        const double halfVal = 0.5 * (pulse[onset] + pulse[peak]);
        for (int i = onset + 1; i <= peak; ++i) {
            if (pulse[i] >= halfVal) return i;
        }
        return std::clamp((onset + peak) / 2, 0, N - 1);
    }

} // namespace ecg_markers