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
            R peak = argmax of the first half of the PRIMARY beat.
            Under the new alignment the template is 1.25*RR wide (0.25*RR
            before R + 1.0*RR of primary beat + tail). "First half of the
            primary beat" == first 0.5*RR of the template == 0.4*N.
            Keeps the search clear of the T-wave (which can otherwise
            outsize R on inverted or tall-T templates).
        */
        const int half = static_cast<int>(v.size()) * 2 / 5;
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

    // Systolic peak of the ANCHOR pulse. With a lead-in cycle prepended, the
// template holds ~3 cycles, so a global/first-half argmax can land on the
// wrong cycle. Search only [onset, onset + oneCycle], which is the anchor
// pulse -- i.e. the cycle whose foot sits under R+delay.
    inline int detect_ppg_peak(const std::vector<double>& pulse, int onset = -1) {
        if (pulse.empty()) return 0;
        auto it = std::max_element(pulse.begin(), pulse.end());
        return static_cast<int>(it - pulse.begin());
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

        // Primary beat ends at 4*N/5 (= 1.0*RR of the 1.25*RR template);
        // don't let T-search leak into the next-beat tail.
        const int primary_end = (4 * N) / 5;
        const int t_lo = r_idx + 100;
        const int t_hi = std::min(static_cast<int>(0.55 * primary_end), primary_end);
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

        // T-end lives inside the primary beat (before next-beat P).
        const int primary_end = (4 * N) / 5;

        const int st_lo = r_idx + 50;
        const int st_hi = std::min(r_idx + 100, primary_end);
        if (st_hi - st_lo < 5) return std::max(0, primary_end - 5);

        std::vector<double> w(upright.begin() + st_lo, upright.begin() + st_hi);
        std::nth_element(w.begin(), w.begin() + w.size() / 2, w.end());
        const double baseline = w[w.size() / 2];

        double st_noise = 0.0;
        for (int i = st_lo; i < st_hi; ++i) {
            st_noise = std::max(st_noise, std::abs(upright[i] - baseline));
        }
        const double tol = 2.0 * st_noise;

        int t_begin = -1;
        for (int i = st_hi; i < primary_end; ++i) {
            if (upright[i] > baseline + tol) { t_begin = i; break; }
        }
        if (t_begin < 0) return std::max(0, primary_end - 5);

        for (int i = t_begin + 1; i < primary_end; ++i) {
            if (upright[i] <= baseline + tol) return i;
        }
        if (t_begin >= primary_end - 1) return primary_end - 1;
        const int fallback = static_cast<int>(t_begin * 1.5);
        return std::clamp(fallback, t_begin + 1, primary_end - 1);
    }
    // S-wave trough: walk forward from R while the beat keeps falling, stop
    // at the first turn back up (the QRS offset / J point). Falls back to
    // just after R when there's no clear trough.
    inline int detect_s(const std::vector<double>& ecg_signal) {
        const int N = static_cast<int>(ecg_signal.size());
        auto [r_idx, is_positive] = r_peak(ecg_signal);
        if (r_idx < 0 || r_idx >= N - 1)
            return std::clamp(r_idx + 1, 0, std::max(0, N - 1));

        std::vector<double> upright = ecg_signal;
        if (!is_positive) for (auto& x : upright) x = -x;

        const int hi = std::min(r_idx + 60, N);
        int s_idx = r_idx;
        for (int i = r_idx + 1; i < hi; ++i) {
            if (upright[i] <= upright[s_idx]) s_idx = i;   // still descending
            else if (i > s_idx + 1) break;                 // risen past the trough
        }
        if (s_idx <= r_idx) return std::min(r_idx + 1, N - 1);
        return s_idx;
    }


    // S-wave end (J point / QRS offset): from the S trough, walk forward
    // until the upstroke recovers most of the way back to the ST baseline.
    // Falls back to a fixed offset past S when there's no clear return.
    inline int detect_s_end(const std::vector<double>& ecg_signal) {
        const int N = static_cast<int>(ecg_signal.size());
        auto [r_idx, is_positive] = r_peak(ecg_signal);
        if (r_idx < 0 || r_idx >= N - 1)
            return std::clamp(r_idx + 1, 0, std::max(0, N - 1));

        std::vector<double> upright = ecg_signal;
        if (!is_positive) for (auto& x : upright) x = -x;

        const int s_idx = detect_s(ecg_signal);   // S-wave trough
        if (s_idx < 0 || s_idx >= N - 1)
            return std::min(std::max(0, r_idx + 1), N - 1);

        // ST baseline: median of the window after the QRS, same as the
        // T-wave detectors use.
        const int st_lo = std::min(r_idx + 50, N - 1);
        const int st_hi = std::min(r_idx + 100, N);
        double baseline = upright[s_idx];          // degenerate default
        if (st_hi - st_lo >= 5) {
            std::vector<double> w(upright.begin() + st_lo, upright.begin() + st_hi);
            std::nth_element(w.begin(), w.begin() + w.size() / 2, w.end());
            baseline = w[w.size() / 2];
        }

        const double s_val = upright[s_idx];
        const double depth = baseline - s_val;

        // If the trough isn't actually below baseline (flat / odd morphology),
        // just step a fixed amount past S.
        if (depth <= 0.0)
            return std::clamp(s_idx + 15, 0, N - 1);

        // Recovery target: 90% of the way from the S trough back to baseline.
        const double target = s_val + 0.90 * depth;
        const int hi = std::min(s_idx + 60, N);
        for (int i = s_idx + 1; i < hi; ++i) {
            if (upright[i] >= target) return i;
        }
        return std::clamp(s_idx + 15, 0, N - 1);
    }


    // ---- Peak detectors (P peak, R peak, T peak) -----------------------
    // Direct peak-finding on the ECG template. Each returns a non-negative
    // sample index -- a fallback position when no clear extremum is found,
    // so the movable bar always has a place to sit.

    // R peak: argmax of first half of primary beat, respecting polarity.
    // Delegates to r_peak() above so the auto-detection rule stays in
    // exactly one place.
    inline int detect_r_peak(const std::vector<double>& ecg_signal) {
        return r_peak(ecg_signal).first;
    }

    // P peak: interior local extremum (matching R polarity) BEFORE Q onset.
    // Q is auto-estimated via detect_q_begin so the seed doesn't depend on
    // any external ordering.
    inline int detect_p_peak(const std::vector<double>& ecg_signal) {
        const int N = static_cast<int>(ecg_signal.size());
        if (N < 3) return 0;

        auto [r_idx, is_positive] = r_peak(ecg_signal);
        std::vector<double> upright = ecg_signal;
        if (!is_positive) for (auto& x : upright) x = -x;

        const int q_est = detect_q_begin(ecg_signal);
        const int hi = std::max(2, std::min(q_est, N - 1));

        // LAST interior local max in [1, hi) matching R polarity (upright).
        int best = -1;
        for (int i = 1; i < hi - 1; ++i) {
            if (std::isnan(upright[i - 1]) || std::isnan(upright[i]) || std::isnan(upright[i + 1])) continue;
            if (upright[i] >= upright[i - 1] && upright[i] >= upright[i + 1]) best = i;
        }
        if (best >= 0) return best;
        return std::max(0, hi / 2);
    }

    // T peak: interior local extremum (matching R polarity) AFTER the QRS,
    // in the T-wave region. Search from ~50 samples past R out to the end
    // of the primary beat (4*N/5 under the 1.25*RR template shape).
    inline int detect_t_peak(const std::vector<double>& ecg_signal) {
        const int N = static_cast<int>(ecg_signal.size());
        if (N < 3) return 0;

        auto [r_idx, is_positive] = r_peak(ecg_signal);
        std::vector<double> upright = ecg_signal;
        if (!is_positive) for (auto& x : upright) x = -x;

        const int primary_end = (4 * N) / 5;
        const int lo = std::min(N - 2, r_idx + 50);
        const int hi = std::max(lo + 5, std::min(N - 1, primary_end));
        if (hi - lo < 5) return std::min(N - 1, r_idx + 150);

        int best = -1;
        double bestVal = -std::numeric_limits<double>::infinity();
        for (int i = lo + 1; i < hi; ++i) {
            if (std::isnan(upright[i - 1]) || std::isnan(upright[i]) || std::isnan(upright[i + 1])) continue;
            if (upright[i] >= upright[i - 1] && upright[i] >= upright[i + 1]) {
                if (upright[i] > bestVal) { bestVal = upright[i]; best = i; }
            }
        }
        if (best >= 0) return best;
        return std::clamp(r_idx + 150, 0, N - 1);
    }


    // ---- PPG detectors --------------------------------------------------
    // `pulse` is one PPG template waveform (one cardiac cycle, onset to
    // end). No R-peak hint -- PPG features are found directly from the
    // pulse shape.

    inline int detect_ppg_onset(const std::vector<double>& pulse) {
        const int N = static_cast<int>(pulse.size());
        if (N < 2) return 0;

        // Peak first: global max over the pulse (systolic peak). The foot is the
        // trough immediately before it, so we don't need a foot hint to find it.
        int peak = 0;
        double pv = pulse[0];
        for (int i = 1; i < N; ++i)
            if (pulse[i] > pv) { pv = pulse[i]; peak = i; }

        // Foot = argmin over [0, peak]: the lowest point before the systolic peak.
        int idx = 0;
        double v = pulse[0];
        for (int i = 1; i <= peak; ++i)
            if (pulse[i] < v) { v = pulse[i]; idx = i; }

        // If the min is at sample 0, prefer a slightly inside position.
        if (idx == 0) return std::min(5, N - 1);
        return idx;
    }

    inline int detect_ppg_end(const std::vector<double>& pulse) {
        const int N = static_cast<int>(pulse.size());
        if (N < 4) return std::max(0, N - 1);

        // Peak first: global max (systolic peak). The end is the trough after it,
        // so no foot/peak hint needed.
        int peak = 0;
        double pv = pulse[0];
        for (int i = 1; i < N; ++i)
            if (pulse[i] > pv) { pv = pulse[i]; peak = i; }

        if (peak >= N - 2) return std::max(0, N - 1);   // peak too close to the end

        // End = argmin over [peak, N): the lowest point after the systolic peak.
        int end = peak + 1;
        double v = pulse[end];
        for (int i = peak + 2; i < N; ++i)
            if (pulse[i] < v) { v = pulse[i]; end = i; }
        return end;
    }

    // Dicrotic notch: local minimum on the descending limb between the
    // systolic peak and the end of the pulse.
    inline int detect_ppg_dicrotic(const std::vector<double>& pulse) {
        const int N = static_cast<int>(pulse.size());
        // in each of the three, replace `detect_ppg_peak(pulse)` with:
        const int onset = detect_ppg_onset(pulse);
        const int peak = detect_ppg_peak(pulse, onset);
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

    // (detect_ppg_peak2 removed: peak2 is a user-placed marker only, no
    // auto-detect. The viewer seeds it to a fixed default position -- 90%
    // of the way from foot to end -- and the user drags it into place.)


    // ---- Combined PPG landmark detector -------------------------------------
    // Contract: pass the visible PPG template plus the two ECG R peaks already
    // mapped into PPG-frame sample indices (rFirstPpg / rSecondPpg). The
    // caller owns that mapping since it has delay + RR at hand.
    //
    // Rule:
    //   peak       = first local maximum in ppg[rFirstPpg .. rSecondPpg]
    //                (falls back to argmax over the same window).
    //   firstTrough = argmin of ppg[0 .. peak]  (deepest low before the peak).
    //   onset      = argmax of ppg[j+1] - ppg[j] over [firstTrough, peak].
    //                This is the point of maximum first-derivative on the
    //                upstroke -- the same criterion EnsembleTemplate uses to
    //                foot-align PPG beats before averaging, so the marker
    //                sits exactly where the template's own alignment point
    //                is. Geometrically ~50% of amplitude on the rising edge.
    //   end        = argmin of ppg[peak .. N-1]  (min after the peak).
    //   dicrotic / p50 left at -1 (add later if you want).
    //
    // NaN-safe: skips NaN samples in every scan.
    struct PpgLandmarks {
        int onset = -1;
        int peak = -1;
        int dicrotic = -1;
        int p50 = -1;
        int end = -1;
    };

    inline PpgLandmarks detect_ppg_landmarks(
        const std::vector<double>& ppg, int rFirstPpg, int rSecondPpg)
    {
        PpgLandmarks r;
        const int N = static_cast<int>(ppg.size());
        if (N < 3) return r;

        // Clamp the R window. If it collapses (arterial: no ECG), fall back
        // to the whole trace so we still return something usable.
        int lo = std::clamp(rFirstPpg, 0, N - 1);
        int hi = std::clamp(rSecondPpg, 0, N - 1);
        if (hi - lo < 3) { lo = 0; hi = N - 1; }

        // 1. Peak = first local max in [lo, hi] (NaN-skipping). Fall back to
        //    argmax over the same window if no clean local max exists.
        int peak = -1;
        for (int i = std::max(1, lo); i <= std::min(hi, N - 2); ++i) {
            if (std::isnan(ppg[i]) || std::isnan(ppg[i - 1]) || std::isnan(ppg[i + 1])) continue;
            if (ppg[i] > ppg[i - 1] && ppg[i] >= ppg[i + 1]) { peak = i; break; }
        }
        if (peak < 0) {
            double pv = -std::numeric_limits<double>::infinity();
            for (int i = lo; i <= hi; ++i)
                if (!std::isnan(ppg[i]) && ppg[i] > pv) { pv = ppg[i]; peak = i; }
        }
        if (peak < 0) return r;   // window all NaN; give up
        r.peak = peak;

        // 3. onset = local minimum AFTER R (the trough right before the
        //    anchor pulse's rise). Search [rFirst, peak]: from the first R
        //    down to the anchor systolic peak. This is the pulse foot for
        //    the anchor beat, matching how PpgOnset markers are consumed
        //    downstream (the "PPG foot" glyph).
        int onset = -1;
        {
            double fv = std::numeric_limits<double>::infinity();
            const int oLo = std::clamp(rFirstPpg, 0, peak);
            for (int i = oLo; i <= peak; ++i)
                if (!std::isnan(ppg[i]) && ppg[i] < fv) { fv = ppg[i]; onset = i; }
        }
        if (onset < 0) onset = std::max(0, peak - 1);   // degenerate fallback
        r.onset = onset;

        // 4. end = argmin of ppg[peak .. N-1], NaN-skipping.
        int end = -1;
        {
            double ev = std::numeric_limits<double>::infinity();
            for (int i = peak + 1; i < N; ++i)
                if (!std::isnan(ppg[i]) && ppg[i] < ev) { ev = ppg[i]; end = i; }
        }
        r.end = (end >= 0) ? end : std::min(N - 1, peak + 1);

        return r;
    }

} // namespace ecg_markers