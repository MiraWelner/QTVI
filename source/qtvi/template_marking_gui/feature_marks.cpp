//
// feature_marks.cpp -- implementations for FeatureMarks.
// See feature_marks.hpp for the public interface.
//

#include "feature_marks.hpp"
#include "TemplateBinIO.hpp"    // for TemplateBin (used by seed_all)

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>

// =========================================================================
// Private helpers
// =========================================================================

std::vector<double> FeatureMarks::first_derivative(const std::vector<double>& v) {
    const int N = static_cast<int>(v.size());
    std::vector<double> d1(N, 0.0);
    if (N < 2) return d1;
    for (int i = 1; i < N - 1; ++i) d1[i] = 0.5 * (v[i + 1] - v[i - 1]);
    d1[0] = v[1] - v[0];
    d1[N - 1] = v[N - 1] - v[N - 2];
    return d1;
}

int FeatureMarks::detect_s(const std::vector<double>& ecg_signal) {
    const int N = static_cast<int>(ecg_signal.size());
    auto [r_idx, is_positive] = r_peak(ecg_signal);
    if (r_idx < 0 || r_idx >= N - 1)
        return std::clamp(r_idx + 1, 0, std::max(0, N - 1));

    std::vector<double> upright = ecg_signal;
    if (!is_positive) for (auto& x : upright) x = -x;

    const int hi = std::min(r_idx + 60, N);
    int s_idx = r_idx;
    for (int i = r_idx + 1; i < hi; ++i) {
        if (upright[i] <= upright[s_idx]) s_idx = i;
        else if (i > s_idx + 1) break;
    }
    if (s_idx <= r_idx) return std::min(r_idx + 1, N - 1);
    return s_idx;
}

// =========================================================================
// Fixed
// =========================================================================

std::pair<int, bool> FeatureMarks::r_peak(const std::vector<double>& v) {
    if (v.empty()) return { 0, true };
    // Search the first 3/4 so the next-beat onset in the tail doesn't
    // compete with the actual R spike.
    const auto hi = v.begin() + (v.size() * 3) / 4;
    const double baseline = std::accumulate(v.begin(), hi, 0.0)
        / static_cast<int>(hi - v.begin());
    auto it = std::max_element(v.begin(), hi,
        [baseline](double a, double b) {
            return std::abs(a - baseline) < std::abs(b - baseline);
        });
    return { static_cast<int>(it - v.begin()), (*it >= baseline) };
}

int FeatureMarks::detect_r_peak(const std::vector<double>& ecg_signal) {
    return r_peak(ecg_signal).first;
}

// =========================================================================
// Reactive
// =========================================================================

int FeatureMarks::compute_q_peak(const std::vector<double>& ecg,
    int q_begin, int r_peak_idx)
{
    const int N = static_cast<int>(ecg.size());
    if (N == 0 || q_begin < 0 || r_peak_idx < 0
        || q_begin >= N || r_peak_idx >= N || q_begin > r_peak_idx)
        return -1;

    auto [ridx, is_positive] = r_peak(ecg);
    (void)ridx;
    const bool up = is_positive;

    int q = q_begin;
    for (int i = q_begin; i <= r_peak_idx; ++i) {
        if (std::isnan(ecg[i]) || std::isnan(ecg[q])) continue;
        if (up ? ecg[i] < ecg[q] : ecg[i] > ecg[q]) q = i;
    }
    return q;
}

int FeatureMarks::compute_s_peak(const std::vector<double>& ecg,
    int r_peak_idx, int s_end)
{
    const int N = static_cast<int>(ecg.size());
    if (N == 0 || r_peak_idx < 0 || s_end < 0
        || r_peak_idx >= N || s_end >= N || r_peak_idx > s_end)
        return -1;

    auto [ridx, is_positive] = r_peak(ecg);
    (void)ridx;
    const bool up = is_positive;

    int s = r_peak_idx;
    for (int i = r_peak_idx; i <= s_end; ++i) {
        if (std::isnan(ecg[i]) || std::isnan(ecg[s])) continue;
        if (up ? ecg[i] < ecg[s] : ecg[i] > ecg[s]) s = i;
    }
    return s;
}

FeatureMarks::EcgGlyphs FeatureMarks::compute_ecg_glyphs(
    const std::vector<double>& ecg,
    int p_peak, int q_begin, int r_peak,
    int s_end, int t_peak, int t_end)
{
    EcgGlyphs g;
    g.p_peak = p_peak;
    g.q_begin = q_begin;
    g.r_peak = r_peak;
    g.s_end = s_end;
    g.t_peak = t_peak;
    g.t_end = t_end;
    g.q_peak = compute_q_peak(ecg, q_begin, r_peak);
    g.s_peak = compute_s_peak(ecg, r_peak, s_end);
    return g;
}

FeatureMarks::PpgGlyphs FeatureMarks::compute_ppg_glyphs(
    const std::vector<double>& v, int foot)
{
    PpgGlyphs r;
    r.foot = foot;
    const int N = static_cast<int>(v.size());
    if (foot < 0 || N < 3) return r;

    // P1: local max in first 75% after foot; O at window midpoint if none.
    int p1 = -1, p1O = -1;
    {
        const int p1hi = foot + (3 * (N - foot)) / 4;
        double best = -1e300;
        for (int i = foot + 1; i < p1hi && i < N - 1; ++i) {
            if (std::isnan(v[i]) || std::isnan(v[i - 1]) || std::isnan(v[i + 1])) continue;
            if (v[i] >= v[i - 1] && v[i] >= v[i + 1] && v[i] > best) {
                best = v[i]; p1 = i;
            }
        }
        if (p1 < 0) p1O = (foot + p1hi) / 2;
    }
    r.p1 = p1; r.p1_fallback = p1O;
    const int p1Ref = (p1 >= 0) ? p1 : p1O;

    // P50: sample nearest 50% amplitude on foot -> P1 upslope.
    if (foot >= 0 && p1Ref > foot &&
        !std::isnan(v[foot]) && !std::isnan(v[p1Ref])) {
        const double target = 0.5 * (v[foot] + v[p1Ref]);
        double bestDiff = 1e300;
        int p50 = -1;
        for (int i = foot; i <= p1Ref; ++i) {
            if (std::isnan(v[i])) continue;
            const double d = std::abs(v[i] - target);
            if (d < bestDiff) { bestDiff = d; p50 = i; }
        }
        const double band = 0.25 * std::abs(v[p1Ref] - v[foot]);
        if (p50 < 0 || bestDiff > band) {
            r.p50 = -1;
            r.p50_fallback = (foot + p1Ref) / 2;
        }
        else {
            r.p50 = p50;
        }
    }
    else if (foot >= 0 && p1Ref > foot) {
        r.p50_fallback = (foot + p1Ref) / 2;
    }

    // End: argmin after P1.
    int end = -1;
    if (p1Ref >= 0) {
        double best = 1e300;
        for (int i = p1Ref + 1; i < N; ++i)
            if (!std::isnan(v[i]) && v[i] < best) { best = v[i]; end = i; }
    }
    r.end = end;

    // Dicrotic: least-negative slope in middle 50% of downslope.
    if (p1Ref >= 0 && end > p1Ref + 4) {
        int firstMin = -1;
        for (int i = p1Ref + 1; i < end && i < N - 1; ++i) {
            if (std::isnan(v[i]) || std::isnan(v[i - 1]) || std::isnan(v[i + 1])) continue;
            if (v[i] <= v[i - 1] && v[i] <= v[i + 1]) { firstMin = i; break; }
        }
        const int slopeStart = (firstMin > 0) ? firstMin : p1Ref;
        const int span = end - slopeStart;
        if (span >= 4) {
            const int mid_lo = slopeStart + span / 4;
            const int mid_hi = slopeStart + (3 * span) / 4;
            int bestIdx = -1;
            double bestSlope = -1e300;
            for (int i = std::max(1, mid_lo); i <= std::min(N - 2, mid_hi); ++i) {
                if (std::isnan(v[i - 1]) || std::isnan(v[i + 1])) continue;
                const double slope = 0.5 * (v[i + 1] - v[i - 1]);
                if (slope > bestSlope) { bestSlope = slope; bestIdx = i; }
            }
            if (bestIdx >= 0 && bestSlope >= 0.0) {
                r.dic = bestIdx;
                r.notch_found = true;
            }
            else {
                r.dic_fallback = (slopeStart + end) / 2;
            }
        }
    }

    // P2: local max between notch (or fallback) and end.
    const int p2Lo = (r.dic >= 0) ? r.dic
        : (r.dic_fallback >= 0 ? r.dic_fallback : p1Ref);
    if (p2Lo >= 0 && end > p2Lo + 1) {
        int p2 = -1;
        double best = -1e300;
        for (int i = p2Lo + 1; i < end; ++i) {
            if (std::isnan(v[i]) || std::isnan(v[i - 1]) || std::isnan(v[i + 1])) continue;
            if (v[i] >= v[i - 1] && v[i] >= v[i + 1] && v[i] > best) {
                best = v[i]; p2 = i;
            }
        }
        if (p2 >= 0) r.p2 = p2;
        else         r.p2_fallback = (p2Lo + end) / 2;
    }

    return r;
}

// =========================================================================
// Movable (auto-detected seeds)
// =========================================================================

int FeatureMarks::detect_q_begin(const std::vector<double>& ecg_signal) {
    const int N = static_cast<int>(ecg_signal.size());
    auto [r_idx, is_positive] = r_peak(ecg_signal);

    std::vector<double> upright_signal = ecg_signal;
    if (!is_positive) for (auto& x : upright_signal) x = -x;

    int search_lim = std::min(r_idx - 1, N);
    const auto d1 = first_derivative(upright_signal);

    int qTrough = 0;
    double qVal = upright_signal[0];
    for (int i = 1; i <= search_lim; ++i) {
        if (upright_signal[i] < qVal) { qVal = upright_signal[i]; qTrough = i; }
    }
    if (qTrough < 1) return std::max(0, r_idx - 20);

    for (int i = qTrough - 1; i >= 0; --i) {
        if (d1[i] >= 0.0) return i;
    }
    return qTrough;
}

int FeatureMarks::detect_p_peak(const std::vector<double>& ecg_signal) {
    const int N = static_cast<int>(ecg_signal.size());
    if (N < 3) return 0;

    auto [r_idx, is_positive] = r_peak(ecg_signal);
    std::vector<double> upright = ecg_signal;
    if (!is_positive) for (auto& x : upright) x = -x;

    const int q_est = detect_q_begin(ecg_signal);
    const int hi = std::max(2, std::min(q_est, N - 1));

    int best = -1;
    for (int i = 1; i < hi - 1; ++i) {
        if (std::isnan(upright[i - 1]) || std::isnan(upright[i]) || std::isnan(upright[i + 1])) continue;
        if (upright[i] >= upright[i - 1] && upright[i] >= upright[i + 1]) best = i;
    }
    if (best >= 0) return best;
    return std::max(0, hi / 2);
}

int FeatureMarks::detect_s_end(const std::vector<double>& ecg_signal) {
    const int N = static_cast<int>(ecg_signal.size());
    auto [r_idx, is_positive] = r_peak(ecg_signal);
    if (r_idx < 0 || r_idx >= N - 1)
        return std::clamp(r_idx + 1, 0, std::max(0, N - 1));

    std::vector<double> upright = ecg_signal;
    if (!is_positive) for (auto& x : upright) x = -x;

    const int s_idx = detect_s(ecg_signal);
    if (s_idx < 0 || s_idx >= N - 1)
        return std::min(std::max(0, r_idx + 1), N - 1);

    const int st_lo = std::min(r_idx + 50, N - 1);
    const int st_hi = std::min(r_idx + 100, N);
    double baseline = upright[s_idx];
    if (st_hi - st_lo >= 5) {
        std::vector<double> w(upright.begin() + st_lo, upright.begin() + st_hi);
        std::nth_element(w.begin(), w.begin() + w.size() / 2, w.end());
        baseline = w[w.size() / 2];
    }

    const double s_val = upright[s_idx];
    const double depth = baseline - s_val;
    if (depth <= 0.0)
        return std::clamp(s_idx + 15, 0, N - 1);

    const double target = s_val + 0.90 * depth;
    const int hi = std::min(s_idx + 60, N);
    for (int i = s_idx + 1; i < hi; ++i) {
        if (upright[i] >= target) return i;
    }
    return std::clamp(s_idx + 15, 0, N - 1);
}

int FeatureMarks::detect_t_peak(const std::vector<double>& ecg_signal) {
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

int FeatureMarks::detect_t_end(const std::vector<double>& ecg_signal) {
    const int N = static_cast<int>(ecg_signal.size());
    auto [r_idx, is_positive] = r_peak(ecg_signal);
    if (r_idx < 0) return std::max(0, (2 * N) / 3);

    std::vector<double> upright = ecg_signal;
    if (!is_positive) for (auto& x : upright) x = -x;

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

// -------------------------------------------------------------------------
// PPG detectors
// -------------------------------------------------------------------------

int FeatureMarks::detect_ppg_onset(const std::vector<double>& pulse) {
    const int N = static_cast<int>(pulse.size());
    if (N < 2) return 0;

    int peak = 0;
    double pv = pulse[0];
    for (int i = 1; i < N; ++i)
        if (pulse[i] > pv) { pv = pulse[i]; peak = i; }

    int idx = 0;
    double v = pulse[0];
    for (int i = 1; i <= peak; ++i)
        if (pulse[i] < v) { v = pulse[i]; idx = i; }

    if (idx == 0) return std::min(5, N - 1);
    return idx;
}

int FeatureMarks::detect_ppg_peak(const std::vector<double>& pulse) {
    if (pulse.empty()) return 0;
    auto it = std::max_element(pulse.begin(), pulse.end());
    return static_cast<int>(it - pulse.begin());
}

int FeatureMarks::detect_ppg_end(const std::vector<double>& pulse) {
    const int N = static_cast<int>(pulse.size());
    if (N < 4) return std::max(0, N - 1);

    int peak = 0;
    double pv = pulse[0];
    for (int i = 1; i < N; ++i)
        if (pulse[i] > pv) { pv = pulse[i]; peak = i; }

    if (peak >= N - 2) return std::max(0, N - 1);

    int end = peak + 1;
    double v = pulse[end];
    for (int i = peak + 2; i < N; ++i)
        if (pulse[i] < v) { v = pulse[i]; end = i; }
    return end;
}

int FeatureMarks::detect_ppg_dicrotic(const std::vector<double>& pulse) {
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

// =========================================================================
// Bin-level seed
// =========================================================================
// One entry point for auto-seeding an entire TemplateBin. Was
// seedBinMarkers() in TemplateViewerWindow.cpp; moved here so all
// marker code lives in one class.

namespace {

    // Seed a pulse's markers from the window BETWEEN THE TWO R PEAKS.
    inline void seedPulse(const std::vector<double>& v, int rFirst, int rSecond,
        int& onset, int& peak, int& dicrotic, int& peak2, int& end)
    {
        const int n = static_cast<int>(v.size());
        if (n < 1) return;
        auto cl = [&](int x) { return std::clamp(x, 0, n - 1); };

        int lo = cl(rFirst), hi = cl(rSecond);
        if (hi - lo < 3) { lo = 0; hi = n - 1; }

        if (peak < 0) {
            int p = -1; double pv = -std::numeric_limits<double>::infinity();
            for (int i = lo; i <= hi; ++i)
                if (!std::isnan(v[i]) && v[i] > pv) { pv = v[i]; p = i; }
            if (p < 0) return;
            peak = p;
        }
        if (onset < 0) {
            int f = -1; double fv = std::numeric_limits<double>::infinity();
            for (int i = 0; i <= peak; ++i)
                if (!std::isnan(v[i]) && v[i] < fv) { fv = v[i]; f = i; }
            if (f < 0) f = std::max(0, peak - 1);
            onset = f;
        }
        if (end < 0) {
            int e = -1; double ev = std::numeric_limits<double>::infinity();
            for (int i = peak + 1; i < n; ++i)
                if (!std::isnan(v[i]) && v[i] < ev) { ev = v[i]; e = i; }
            if (e < 0) e = std::min(n - 1, peak + 1);
            end = e;
        }
        if (dicrotic < 0 && cl(end) > cl(onset) + 2) {
            const int base = cl(onset);
            std::vector<double> cyc(v.begin() + base, v.begin() + cl(end) + 1);
            dicrotic = cl(base + FeatureMarks::detect_ppg_dicrotic(cyc));
        }
        if (peak2 < 0 && end > onset) {
            peak2 = cl(onset + (9 * (end - onset)) / 10);
        }
    }

    inline int clampToVisible(int idx, int visN) {
        return std::clamp(idx, 0, visN - 1);
    }

} // anonymous

void FeatureMarks::seed_all(TemplateBin& b, double sampleRate) {

    // ---- PPG ------------------------------------------------------------
    if (b.ppgTemplate.empty()) {
        b.ppg_issue = 2;
        b.ppg_onset = b.ppg_p50 = b.ppg_peak = -1;
        b.ppg_dicrotic = b.ppg_peak2 = b.ppg_end = -1;
        b.ppg_onset_auto = b.ppg_p50_auto = b.ppg_peak_auto = -1;
        b.ppg_dicrotic_auto = b.ppg_peak2_auto = b.ppg_end_auto = -1;
    }
    else if (b.ppg_issue == 1) {
        b.ppg_onset = b.ppg_p50 = b.ppg_peak = -1;
        b.ppg_dicrotic = b.ppg_peak2 = b.ppg_end = -1;
        // Leave *_auto alone -- they're the original auto positions.
    }
    else {
        const std::vector<double>& v = b.ppgTemplate;
        const int N = static_cast<int>(v.size());

        const int rFirst = (sampleRate > 0.0)
            ? static_cast<int>(std::llround(0.3 * sampleRate)) : 0;
        const int footWin = (sampleRate > 0.0)
            ? static_cast<int>(std::llround(0.35 * sampleRate)) : N / 2;
        const int cycleGuess = std::max(1, N - rFirst);

        // foot: argmin in [rFirst, rFirst + footWin)
        int foot = -1;
        {
            double best = std::numeric_limits<double>::infinity();
            const int lo = std::clamp(rFirst, 0, N - 1);
            const int hi = std::min(N, lo + footWin);
            for (int i = lo; i < hi; ++i)
                if (!std::isnan(v[i]) && v[i] < best) { best = v[i]; foot = i; }
        }

        // peak: argmax in [foot+1, foot+cycleGuess)
        int peak = -1;
        if (foot >= 0) {
            double best = -std::numeric_limits<double>::infinity();
            const int hi = std::min(N, foot + cycleGuess);
            for (int i = foot + 1; i < hi; ++i)
                if (!std::isnan(v[i]) && v[i] > best) { best = v[i]; peak = i; }
        }

        // end: argmin in [peak+1, peak+cycleGuess)
        int end = -1;
        if (peak >= 0) {
            double best = std::numeric_limits<double>::infinity();
            const int hi = std::min(N, peak + cycleGuess);
            for (int i = peak + 1; i < hi; ++i)
                if (!std::isnan(v[i]) && v[i] < best) { best = v[i]; end = i; }
        }

        // dicrotic: least-negative slope in middle 50% of downslope
        int dic = -1;
        if (peak >= 0 && end > peak + 4) {
            int firstMin = -1;
            for (int i = peak + 1; i < end && i < N - 1; ++i) {
                if (std::isnan(v[i]) || std::isnan(v[i - 1]) || std::isnan(v[i + 1])) continue;
                if (v[i] <= v[i - 1] && v[i] <= v[i + 1]) { firstMin = i; break; }
            }
            const int slopeStart = (firstMin > 0) ? firstMin : peak;
            const int span = end - slopeStart;
            if (span >= 4) {
                const int mid_lo = slopeStart + span / 4;
                const int mid_hi = slopeStart + (3 * span) / 4;
                int bestIdx = -1;
                double bestSlope = -std::numeric_limits<double>::infinity();
                for (int i = std::max(1, mid_lo); i <= std::min(N - 2, mid_hi); ++i) {
                    if (std::isnan(v[i - 1]) || std::isnan(v[i + 1])) continue;
                    const double slope = 0.5 * (v[i + 1] - v[i - 1]);
                    if (slope > bestSlope) { bestSlope = slope; bestIdx = i; }
                }
                if (bestIdx >= 0 && bestSlope >= 0.0) dic = bestIdx;
                else                                  dic = (slopeStart + end) / 2;
            }
        }

        // p50: nearest 50% amplitude on foot->peak upslope
        int p50 = -1;
        if (foot >= 0 && peak > foot &&
            !std::isnan(v[foot]) && !std::isnan(v[peak])) {
            const double target = 0.5 * (v[foot] + v[peak]);
            double bestDiff = std::numeric_limits<double>::infinity();
            for (int i = foot; i <= peak; ++i) {
                if (std::isnan(v[i])) continue;
                const double d = std::abs(v[i] - target);
                if (d < bestDiff) { bestDiff = d; p50 = i; }
            }
            const double band = 0.25 * std::abs(v[peak] - v[foot]);
            if (bestDiff > band) p50 = (foot + peak) / 2;
        }
        else if (foot >= 0 && peak > foot) {
            p50 = (foot + peak) / 2;
        }

        // peak2: default 90% foot->end
        int peak2 = (foot >= 0 && end > foot)
            ? foot + (9 * (end - foot)) / 10
            : std::max(0, N - 1);

        // ALWAYS populate auto fields (fresh every loadSubject).
        b.ppg_onset_auto = foot;
        b.ppg_peak_auto = peak;
        b.ppg_end_auto = end;
        b.ppg_dicrotic_auto = dic;
        b.ppg_p50_auto = p50;
        b.ppg_peak2_auto = peak2;
        // User fields: only seed when unset (respect prior edits).
        if (b.ppg_onset < 0) b.ppg_onset = foot;
        if (b.ppg_peak < 0) b.ppg_peak = peak;
        if (b.ppg_end < 0) b.ppg_end = end;
        if (b.ppg_dicrotic < 0) b.ppg_dicrotic = dic;
        if (b.ppg_p50 < 0) b.ppg_p50 = p50;
        if (b.ppg_peak2 < 0) b.ppg_peak2 = peak2;
    }

    // ---- ECG (per channel) ---------------------------------------------
    ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };
    for (int c = 0; c < 3; ++c) {
        const auto& ecg = chs[c]->ecgTemplate_raw;
        if (ecg.empty()) {
            b.bad_r_ch[c] = true;
            b.p_peak_ch[c] = b.q_begin_ch[c] = b.r_peak_ch[c] = -1;
            b.s_end_ch[c] = b.t_peak_ch[c] = b.t_end_ch[c] = -1;
            b.p_peak_auto_ch[c] = b.q_begin_auto_ch[c] = b.r_peak_auto_ch[c] = -1;
            b.s_end_auto_ch[c] = b.t_peak_auto_ch[c] = b.t_end_auto_ch[c] = -1;
            continue;
        }

        const int visN = std::max(static_cast<int>(ecg.size()), 2);

        const int p_auto = clampToVisible(detect_p_peak(ecg), visN);
        const int q_auto = clampToVisible(detect_q_begin(ecg), visN);
        const int r_auto = clampToVisible(detect_r_peak(ecg), visN);
        const int s_auto = clampToVisible(detect_s_end(ecg), visN);
        const int tp_auto = clampToVisible(detect_t_peak(ecg), visN);
        const int te_auto = clampToVisible(detect_t_end(ecg), visN);

        // Auto fields always updated.
        b.p_peak_auto_ch[c] = p_auto;
        b.q_begin_auto_ch[c] = q_auto;
        b.r_peak_auto_ch[c] = r_auto;
        b.s_end_auto_ch[c] = s_auto;
        b.t_peak_auto_ch[c] = tp_auto;
        b.t_end_auto_ch[c] = te_auto;

        // User fields: only seed when unset. R peak is auto-only so
        // its user field is always overwritten with the fresh auto.
        if (b.p_peak_ch[c] < 0) b.p_peak_ch[c] = p_auto;
        if (b.q_begin_ch[c] < 0) b.q_begin_ch[c] = q_auto;
        b.r_peak_ch[c] = r_auto;
        if (b.s_end_ch[c] < 0) b.s_end_ch[c] = s_auto;
        if (b.t_peak_ch[c] < 0) b.t_peak_ch[c] = tp_auto;
        if (b.t_end_ch[c] < 0) b.t_end_ch[c] = te_auto;
    }

    // ---- Arterial (ABP / ART / ART_PULM) --------------------------------
    const int rFirstArt = (sampleRate > 0.0)
        ? static_cast<int>(std::llround(0.3 * sampleRate)) : 0;
    const int rSecondArt = rFirstArt + ((sampleRate > 0.0)
        ? static_cast<int>(std::llround(0.75 * sampleRate)) : 0);

    auto seedArterial = [&](const std::vector<double>& trace, uint8_t& issue,
        int& onset, int& peak, int& dicrotic, int& peak2, int& end,
        int& onset_auto, int& peak_auto, int& dic_auto, int& p2_auto, int& end_auto)
        {
            if (trace.empty()) {
                issue = 2;
                onset = peak = dicrotic = peak2 = end = -1;
                onset_auto = peak_auto = dic_auto = p2_auto = end_auto = -1;
                return;
            }
            if (issue == 1) {
                onset = peak = dicrotic = peak2 = end = -1;
                return;
            }
            int aOn = -1, aPk = -1, aDic = -1, aP2 = -1, aEnd = -1;
            seedPulse(trace, rFirstArt, rSecondArt, aOn, aPk, aDic, aP2, aEnd);
            onset_auto = aOn; peak_auto = aPk; dic_auto = aDic;
            p2_auto = aP2; end_auto = aEnd;
            if (onset < 0) onset = aOn;
            if (peak < 0) peak = aPk;
            if (dicrotic < 0) dicrotic = aDic;
            if (peak2 < 0) peak2 = aP2;
            if (end < 0) end = aEnd;
        };

    seedArterial(b.abpTemplate, b.abp_issue,
        b.abp_onset, b.abp_peak, b.abp_dicrotic, b.abp_peak2, b.abp_end,
        b.abp_onset_auto, b.abp_peak_auto, b.abp_dicrotic_auto,
        b.abp_peak2_auto, b.abp_end_auto);
    seedArterial(b.artTemplate, b.art_issue,
        b.art_onset, b.art_peak, b.art_dicrotic, b.art_peak2, b.art_end,
        b.art_onset_auto, b.art_peak_auto, b.art_dicrotic_auto,
        b.art_peak2_auto, b.art_end_auto);
    seedArterial(b.artPulmTemplate, b.art_pulm_issue,
        b.art_pulm_onset, b.art_pulm_peak, b.art_pulm_dicrotic,
        b.art_pulm_peak2, b.art_pulm_end,
        b.art_pulm_onset_auto, b.art_pulm_peak_auto, b.art_pulm_dicrotic_auto,
        b.art_pulm_peak2_auto, b.art_pulm_end_auto);
}