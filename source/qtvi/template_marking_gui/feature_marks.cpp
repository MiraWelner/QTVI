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
#include <vector>

// =========================================================================
// File-scope helpers for the reactive ECG glyphs.
// =========================================================================
namespace {

    // Half-window in samples for a +/-0.05 s search around a user marker.
    inline int win_005s(double fs) {
        return std::max(1, static_cast<int>(std::lround(0.05 * fs)));
    }

    // Least-squares polynomial fit (degree d) over samples [lo,hi] of v, with
    // x measured from lo. Returns coeffs c[0..d] (c[0] + c[1]x + c[2]x^2 ...).
    // Small normal-equation solve via Gaussian elimination with partial pivot.
    inline std::vector<double> polyfit(const std::vector<double>& v,
        int lo, int hi, int d)
    {
        const int n = hi - lo + 1;
        const int m = d + 1;
        std::vector<double> A(m * m, 0.0), b(m, 0.0);
        for (int k = 0; k < n; ++k) {
            const double y = v[lo + k];
            if (std::isnan(y)) continue;
            const double x = k;
            std::vector<double> px(2 * m - 1);
            double xp = 1.0;
            for (int j = 0; j < 2 * m - 1; ++j) { px[j] = xp; xp *= x; }
            for (int i = 0; i < m; ++i) {
                b[i] += px[i] * y;
                for (int j = 0; j < m; ++j) A[i * m + j] += px[i + j];
            }
        }
        for (int i = 0; i < m; ++i) {
            int piv = i;
            for (int r = i + 1; r < m; ++r)
                if (std::abs(A[r * m + i]) > std::abs(A[piv * m + i])) piv = r;
            for (int j = 0; j < m; ++j) std::swap(A[i * m + j], A[piv * m + j]);
            std::swap(b[i], b[piv]);
            const double diag = A[i * m + i];
            if (std::abs(diag) < 1e-12) continue;
            for (int r = 0; r < m; ++r) {
                if (r == i) continue;
                const double f = A[r * m + i] / diag;
                for (int j = 0; j < m; ++j) A[r * m + j] -= f * A[i * m + j];
                b[r] -= f * b[i];
            }
        }
        std::vector<double> c(m, 0.0);
        for (int i = 0; i < m; ++i) {
            const double diag = A[i * m + i];
            c[i] = (std::abs(diag) > 1e-12) ? b[i] / diag : 0.0;
        }
        return c;
    }

    // Evaluate a fitted polynomial (coeffs c[0..d]) at offset x.
    inline double poly_eval(const std::vector<double>& c, double x) {
        double s = 0.0, xp = 1.0;
        for (double coef : c) { s += coef * xp; xp *= x; }
        return s;
    }

    // "Elbow" of a degree-d polynomial fit over [lo,hi]: the sample whose
    // fitted value deviates most from the straight chord joining the fitted
    // endpoints. (Max-curvature of a cubic is useless here -- its 2nd
    // derivative is linear, so the extremum is always at a window edge.)
    inline int fit_elbow(const std::vector<double>& v, int lo, int hi, int d) {
        if (hi - lo < d + 1) return (lo + hi) / 2;
        const auto c = polyfit(v, lo, hi, d);
        const int n = hi - lo;
        const double y0 = poly_eval(c, 0.0);
        const double y1 = poly_eval(c, static_cast<double>(n));
        int best = lo; double bd = -1.0;
        for (int k = 0; k <= n; ++k) {
            const double chord = y0 + (y1 - y0) * (static_cast<double>(k) / n);
            const double dev = std::abs(poly_eval(c, static_cast<double>(k)) - chord);
            if (dev > bd) { bd = dev; best = lo + k; }
        }
        return best;
    }

} // anonymous


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
    const int N = static_cast<int>(v.size());
    if (N == 0) return { 0, true };

    // R is the first deflection after 0.176 of the window (its expected column
    // by snip geometry): scan forward and return the first local extremum.
    // A local max => upright R (positive); a local min => inverted R (negative).
    const int lo = std::clamp(static_cast<int>(std::lround(0.176 * N)), 1, N - 1);
    const int hi = std::clamp(static_cast<int>(std::lround(0.30 * N)), lo, N - 1);
    for (int i = lo; i < hi; ++i) {
        if (std::isnan(v[i - 1]) || std::isnan(v[i]) || std::isnan(v[i + 1])) continue;
        const bool isMax = v[i] >= v[i - 1] && v[i] >= v[i + 1];
        const bool isMin = v[i] <= v[i - 1] && v[i] <= v[i + 1];
        if (isMax && !isMin) return { i, true };    // first upward deflection
        if (isMin && !isMax) return { i, false };   // first downward deflection
    }
    return { lo, true };
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

// -------------------------------------------------------------------------
// Reactive ECG X-glyphs: each is auto-computed but tracks the user's movable
// markers live. Windows are +/-0.05 s around the relevant user marker.
// -------------------------------------------------------------------------

// P wave = max value within +/-0.05 s of the user's P-peak marker.
int FeatureMarks::compute_p_wave(const std::vector<double>& v, int pUser, double fs) {
    const int N = static_cast<int>(v.size());
    if (pUser < 0 || N == 0) return -1;
    const int w = win_005s(fs);
    const int lo = std::max(0, pUser - w), hi = std::min(N - 1, pUser + w);
    int best = lo; double bv = -std::numeric_limits<double>::infinity();
    for (int i = lo; i <= hi; ++i)
        if (!std::isnan(v[i]) && v[i] > bv) { bv = v[i]; best = i; }
    return best;
}

// T peak = max value between the user's T-begin and T-end markers.
int FeatureMarks::compute_t_wave(const std::vector<double>& v, int tBegin, int tEnd) {
    const int N = static_cast<int>(v.size());
    if (tBegin < 0 || tEnd < 0 || tBegin >= N || tEnd >= N) return -1;
    const int lo = std::min(tBegin, tEnd), hi = std::max(tBegin, tEnd);
    int best = lo; double bv = -std::numeric_limits<double>::infinity();
    for (int i = lo; i <= hi; ++i)
        if (!std::isnan(v[i]) && v[i] > bv) { bv = v[i]; best = i; }
    return best;
}

// R wave = argmax |v - baseline| within [q_begin, s_end] (the user markers),
// baseline = mean over that same window. Same rule as r_peak(), new window.
int FeatureMarks::compute_r_wave(const std::vector<double>& v, int qBegin, int sEnd) {
    const int N = static_cast<int>(v.size());
    if (qBegin < 0 || sEnd < 0 || qBegin >= N || sEnd >= N || qBegin >= sEnd)
        return -1;
    double base = 0.0; int n = 0;
    for (int i = 0; i < N; ++i)
        if (!std::isnan(v[i])) { base += v[i]; ++n; }
    if (n == 0) return -1;
    base /= n;
    int best = qBegin; double bd = -1.0;
    for (int i = qBegin; i <= sEnd; ++i) {
        if (std::isnan(v[i])) continue;
        const double d = std::abs(v[i] - base);
        if (d > bd) { bd = d; best = i; }
    }
    return best;
}

// S end = recovery dropoff/knee (steepest upstroke back toward baseline, then
// the point where that slope flattens) within +/-0.05 s of the user's S-end
// marker.
int FeatureMarks::compute_s_end(const std::vector<double>& v, int sUser, double fs) {
    const int N = static_cast<int>(v.size());
    if (sUser < 0 || N < 3) return sUser;
    const int w = win_005s(fs);
    const int lo = std::max(1, sUser - w), hi = std::min(N - 2, sUser + w);
    if (hi <= lo) return std::clamp(sUser, 0, N - 1);

    double mx = 0.0; int mi = lo;
    for (int i = lo; i <= hi; ++i) {
        if (std::isnan(v[i - 1]) || std::isnan(v[i + 1])) continue;
        const double s = 0.5 * (v[i + 1] - v[i - 1]);
        if (s > mx) { mx = s; mi = i; }
    }
    if (mx <= 0.0) return std::clamp(sUser, 0, N - 1);   // no clear recovery
    for (int i = mi; i <= hi; ++i) {
        if (std::isnan(v[i - 1]) || std::isnan(v[i + 1])) continue;
        const double s = 0.5 * (v[i + 1] - v[i - 1]);
        if (s <= 0.2 * mx) return i;                      // slope has flattened
    }
    return hi;
}

// Q onset = knee (max |curvature|) of a cubic fit within +/-0.05 s of the
// user's Q-begin marker.
int FeatureMarks::compute_q_onset(const std::vector<double>& v, int qUser, double fs) {
    const int N = static_cast<int>(v.size());
    if (qUser < 0 || N < 4) return qUser;
    const int w = win_005s(fs);

    // Orient so Q is a downward dip.
    auto [r_idx, is_positive] = r_peak(v); (void)r_idx;
    std::vector<double> u = v;
    if (!is_positive) for (auto& x : u) x = -x;

    // Q trough nearest the user's mark.
    const int a = std::max(0, qUser - w), b = std::min(N - 1, qUser + w);
    int qmin = a;
    for (int i = a; i <= b; ++i)
        if (!std::isnan(u[i]) && u[i] < u[qmin]) qmin = i;

    auto d1 = [&](int i) {
        if (i <= 0)     return u[1] - u[0];
        if (i >= N - 1) return u[N - 1] - u[N - 2];
        return 0.5 * (u[i + 1] - u[i - 1]);
        };

    // Onset = walk left from the trough until the descent flattens. "Flattened"
    // is judged relative to the steepest descent slope, so a gentle P-wave
    // downslope before Q is not mistaken for the descent into Q (a plain
    // slope >= 0 test overshoots to the P peak when the PR segment isn't flat).
    double steep = 0.0;
    for (int i = qmin; i > std::max(0, qmin - 2 * w); --i)
        steep = std::min(steep, d1(i));
    const double thr = 0.10 * steep;                 // steep <= 0, so thr <= 0
    const int stop = std::max(0, qmin - 4 * w);
    for (int i = qmin - 1; i >= stop; --i)
        if (d1(i) >= thr) return i;
    return std::max(0, qmin - 1);
}

// T end = knee (max |curvature|) of a cubic fit within +/-0.05 s of the
// user's T-end marker.
int FeatureMarks::compute_t_end(const std::vector<double>& v, int tEndUser, double fs) {
    const int N = static_cast<int>(v.size());
    if (tEndUser < 0 || N < 4) return tEndUser;
    const int w = win_005s(fs);
    const int lo = std::max(0, tEndUser - w), hi = std::min(N - 1, tEndUser + w);
    return fit_elbow(v, lo, hi, 3);
}

FeatureMarks::EcgGlyphs FeatureMarks::compute_ecg_glyphs(
    const std::vector<double>& ecg,
    int p_peak, int q_begin, int s_end, int t_begin, int t_end, double fs)
{
    EcgGlyphs g;
    g.p_wave = compute_p_wave(ecg, p_peak, fs);
    g.q_onset = compute_q_onset(ecg, q_begin, fs);
    g.r_wave = compute_r_wave(ecg, q_begin, s_end);   // window = user q_begin..s_end
    const int N = static_cast<int>(ecg.size());
    g.s_end = (s_end >= 0 && s_end < N) ? s_end : -1;  // S = the user's marker (no math)
    g.t_peak = compute_t_wave(ecg, t_begin, t_end);   // max between user T begin/end
    g.t_end = (t_end >= 0 && t_end < N) ? t_end : -1;  // T end = the user's marker
    return g;
}

FeatureMarks::PpgGlyphs FeatureMarks::compute_ppg_glyphs(
    const std::vector<double>& v, int foot, int dic, int peak2)
{
    PpgGlyphs g;
    const int N = static_cast<int>(v.size());
    if (N < 3) return g;
    auto cl = [&](int x) { return std::clamp(x, 0, N - 1); };

    // Foot = minimum near the foot marker (the marker is seeded in the trough,
    // so this small-window search refines/tracks it).
    if (foot >= 0) {
        const int w = std::max(3, N / 30);
        const int lo = cl(foot - w), hi = cl(foot + w);
        int fmin = lo; double best = std::numeric_limits<double>::infinity();
        for (int i = lo; i <= hi; ++i)
            if (!std::isnan(v[i]) && v[i] < best) { best = v[i]; fmin = i; }
        g.foot = fmin;
    }

    // Peak #1 (systolic) = max between the foot and the dicrotic-notch marker.
    if (g.foot >= 0 && dic > g.foot) {
        int p1 = g.foot; double best = -std::numeric_limits<double>::infinity();
        for (int i = g.foot; i <= cl(dic); ++i)
            if (!std::isnan(v[i]) && v[i] > best) { best = v[i]; p1 = i; }
        g.p1 = p1;
    }

    // Dicrotic notch and Peak #2 are placed manually -> glyph tracks the marker.
    if (dic >= 0) g.dic = cl(dic);
    if (peak2 >= 0) g.p2 = cl(peak2);
    g.notch_found = (dic >= 0);

    // End = first local minimum after the Peak #2 marker.
    if (peak2 >= 0) {
        for (int i = cl(peak2) + 1; i < N - 1; ++i) {
            if (std::isnan(v[i - 1]) || std::isnan(v[i]) || std::isnan(v[i + 1])) continue;
            if (v[i] <= v[i - 1] && v[i] <= v[i + 1]) { g.end = i; break; }
        }
    }

    // Not used in this scheme.
    g.p50 = -1;
    return g;
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

    // Max value before Q begin.
    const int q_est = detect_q_begin(ecg_signal);
    const int hi = std::max(1, std::min(q_est, N));
    int best = -1; double bv = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < hi; ++i)
        if (!std::isnan(upright[i]) && upright[i] > bv) { bv = upright[i]; best = i; }
    return (best >= 0) ? best : 0;
}

int FeatureMarks::detect_t_begin(const std::vector<double>& ecg_signal) {
    // First local maximum right of S end (the T wave), then walk left to its
    // foot. begin/end bracket that peak.
    const int N = static_cast<int>(ecg_signal.size());
    if (N < 3) return 0;
    auto [r_idx, is_positive] = r_peak(ecg_signal);
    (void)r_idx;
    std::vector<double> upright = ecg_signal;
    if (!is_positive) for (auto& x : upright) x = -x;

    const int s_end = detect_s_end(ecg_signal);
    const int lo = std::clamp(s_end + 1, 1, N - 2);
    const int hi = std::min(N - 1, (2 * N) / 3);       // first 2/3 only

    int tPeak = -1;
    double tBest = -std::numeric_limits<double>::infinity();
    for (int i = lo; i < hi; ++i) {
        if (std::isnan(upright[i - 1]) || std::isnan(upright[i]) || std::isnan(upright[i + 1])) continue;
        if (upright[i] > upright[i - 1] && upright[i] >= upright[i + 1]) {
            if (upright[i] > tBest) { tBest = upright[i]; tPeak = i; }   // tallest = the T wave, not the first ST bump
        }
    }
    if (tPeak < 0) return std::clamp(s_end + (int)std::lround(0.1 * N), 0, N - 1);

    int i = tPeak;   // left foot
    while (i - 1 > s_end && !std::isnan(upright[i - 1]) && upright[i - 1] < upright[i]) --i;
    return i;
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

int FeatureMarks::detect_t_end(const std::vector<double>& ecg_signal) {
    // First local maximum right of S end, then walk right to its foot.
    const int N = static_cast<int>(ecg_signal.size());
    if (N < 3) return std::max(0, N - 1);
    auto [r_idx, is_positive] = r_peak(ecg_signal);
    (void)r_idx;
    std::vector<double> upright = ecg_signal;
    if (!is_positive) for (auto& x : upright) x = -x;

    const int s_end = detect_s_end(ecg_signal);
    const int lo = std::clamp(s_end + 1, 1, N - 2);
    const int hi = std::min(N - 1, (2 * N) / 3);       // first 2/3 only

    int tPeak = -1;
    double tBest = -std::numeric_limits<double>::infinity();
    for (int i = lo; i < hi; ++i) {
        if (std::isnan(upright[i - 1]) || std::isnan(upright[i]) || std::isnan(upright[i + 1])) continue;
        if (upright[i] > upright[i - 1] && upright[i] >= upright[i + 1]) {
            if (upright[i] > tBest) { tBest = upright[i]; tPeak = i; }   // tallest = the T wave, not the first ST bump
        }
    }
    if (tPeak < 0) return std::clamp(s_end + (int)std::lround(0.1 * N), 0, N - 1);

    int i = tPeak;   // right foot
    while (i + 1 < hi && !std::isnan(upright[i + 1]) && upright[i + 1] < upright[i]) ++i;
    return std::clamp(std::max(i, s_end + 1), 0, N - 1);   // T end must sit right of S-end
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

int FeatureMarks::detect_ppg_t80(const std::vector<double>& pulse) {
    const int N = static_cast<int>(pulse.size());
    const int foot = detect_ppg_onset(pulse);
    const int peak = detect_ppg_peak(pulse);
    if (foot < 0 || peak <= foot
        || std::isnan(pulse[foot]) || std::isnan(pulse[peak]))
        return std::clamp((foot + peak) / 2, 0, std::max(0, N - 1));

    const double target = pulse[foot] + 0.80 * (pulse[peak] - pulse[foot]);
    int best = foot; double bestDiff = std::numeric_limits<double>::infinity();
    for (int i = foot; i <= peak; ++i) {
        if (std::isnan(pulse[i])) continue;
        const double d = std::abs(pulse[i] - target);
        if (d < bestDiff) { bestDiff = d; best = i; }
    }
    return best;
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
        b.ppg_onset = b.ppg_p50 = b.ppg_t80 = b.ppg_peak = -1;
        b.ppg_dicrotic = b.ppg_peak2 = b.ppg_end = -1;
        b.ppg_onset_auto = b.ppg_p50_auto = b.ppg_t80_auto = b.ppg_peak_auto = -1;
        b.ppg_dicrotic_auto = b.ppg_peak2_auto = b.ppg_end_auto = -1;
    }
    else if (b.ppg_issue == 1) {
        b.ppg_onset = b.ppg_p50 = b.ppg_t80 = b.ppg_peak = -1;
        b.ppg_dicrotic = b.ppg_peak2 = b.ppg_end = -1;
        // Leave *_auto alone -- they're the original auto positions.
    }
    else {
        const std::vector<double>& v = b.ppgTemplate;
        const int N = static_cast<int>(v.size());
        // Visible PPG window: the display clips the PPG to the ECG time axis
        // (ppgStartSample()==0), so only min(ppg_len, ecg_len) samples show.
        // Seed everything within W so no marker lands past the screen edge.
        int W = N;
        {
            const ChannelTemplateData* cc[3] = { &b.ch1, &b.ch2, &b.ch3 };
            int mn = N; bool any = false;
            for (const auto* ch : cc) {
                const int l = static_cast<int>(ch->ecgTemplate_raw.size());
                if (l > 0) { mn = any ? std::min(mn, l) : l; any = true; }
            }
            if (any) W = std::min(N, mn);
        }
        W = std::max(2, W);
        auto pct = [&](double f) {
            return std::clamp(static_cast<int>(std::lround(f * W)), 0, W - 1);
            };

        // Movable-marker spawn points.
        // Seeds only; the user drags from here.
        const int dic = pct(0.60);   // Dicrotic notch  60%
        const int peak2 = pct(0.65);   // Peak #2         65%

        // Systolic peak = argmax over the visible window; Foot = the minimum
        // BEFORE that peak (seeds the foot bar in the trough); End = the
        // minimum between the peak and the end of the visible window; T80 =
        // 80% of the way DOWN from peak toward end; p50 = temporal midpoint
        // foot..peak.
        int peak = 0;
        {
            double best = -std::numeric_limits<double>::infinity();
            for (int i = 0; i < W; ++i)
                if (!std::isnan(v[i]) && v[i] > best) { best = v[i]; peak = i; }
        }
        int foot = 0;
        {
            const int lo = std::min(pct(0.10), peak);   // ignore the first 10%
            double best = std::numeric_limits<double>::infinity();
            for (int i = lo; i <= peak; ++i)
                if (!std::isnan(v[i]) && v[i] < best) { best = v[i]; foot = i; }
        }
        int end = std::min(W - 1, peak + 1);
        {
            double best = std::numeric_limits<double>::infinity();
            for (int i = peak + 1; i < W; ++i)   // min between peak and end of screen
                if (!std::isnan(v[i]) && v[i] < best) { best = v[i]; end = i; }
        }

        int t80 = std::clamp((peak + end) / 2, 0, W - 1);
        if (end > peak && !std::isnan(v[peak]) && !std::isnan(v[end])) {
            const double target = v[peak] + 0.80 * (v[end] - v[peak]);  // 80% down
            double bestDiff = std::numeric_limits<double>::infinity();
            for (int i = peak; i <= end; ++i) {
                if (std::isnan(v[i])) continue;
                const double d = std::abs(v[i] - target);
                if (d < bestDiff) { bestDiff = d; t80 = i; }
            }
        }

        // P50: temporal midpoint between onset (foot) and systolic peak.
        int p50 = (foot + peak) / 2;

        // Auto fields always updated.
        b.ppg_onset_auto = foot;
        b.ppg_peak_auto = peak;
        b.ppg_end_auto = end;
        b.ppg_dicrotic_auto = dic;
        b.ppg_p50_auto = p50;
        b.ppg_t80_auto = t80;
        b.ppg_peak2_auto = peak2;

        // Movable markers: seed only when unset (respect prior edits).
        if (b.ppg_onset < 0)    b.ppg_onset = foot;
        if (b.ppg_dicrotic < 0) b.ppg_dicrotic = dic;
        if (b.ppg_peak2 < 0)    b.ppg_peak2 = peak2;
        if (b.ppg_t80 < 0)    b.ppg_t80 = t80;
        if (b.ppg_end < 0)      b.ppg_end = end;
        // Auto-only markers: always refresh.
        b.ppg_peak = peak;
        b.ppg_p50 = p50;
    }

    // ---- ECG (per channel) ---------------------------------------------
    // Initial user-marker spawn points as fractions of the plot length.
    // Seeds only; the user drags from here. R stays auto-detected.
    ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };
    for (int c = 0; c < 3; ++c) {
        const auto& ecg = chs[c]->ecgTemplate_raw;
        if (ecg.empty()) {
            b.bad_r_ch[c] = true;
            b.p_peak_ch[c] = b.q_begin_ch[c] = b.r_peak_ch[c] = -1;
            b.s_end_ch[c] = b.t_begin_ch[c] = b.t_end_ch[c] = -1;
            b.p_peak_auto_ch[c] = b.q_begin_auto_ch[c] = b.r_peak_auto_ch[c] = -1;
            b.s_end_auto_ch[c] = b.t_begin_auto_ch[c] = b.t_end_auto_ch[c] = -1;
            continue;
        }

        const int visN = std::max(static_cast<int>(ecg.size()), 2);
        auto cl = [&](int x) { return clampToVisible(x, visN); };
        auto pct = [&](double f) { return static_cast<int>(std::lround(f * ecg.size())); };

        const int p_auto = cl(pct(0.05));            // P       5% in
        const int q_auto = cl(pct(0.10));            // Q begin 10% in
        const int r_auto = cl(detect_r_peak(ecg));   // R stays auto-detected
        const int s_auto = cl(pct(0.25));            // S end   25% in
        const int tp_auto = cl(detect_t_begin(ecg));  // T begin bar: left foot of the T wave (auto)
        const int te_auto = cl(detect_t_end(ecg));    // T end:   right foot of the T wave (auto)

        // Auto fields always updated.
        b.p_peak_auto_ch[c] = p_auto;
        b.q_begin_auto_ch[c] = q_auto;
        b.r_peak_auto_ch[c] = r_auto;
        b.s_end_auto_ch[c] = s_auto;
        b.t_begin_auto_ch[c] = tp_auto;
        b.t_end_auto_ch[c] = te_auto;

        // User fields: only seed when unset. R peak is auto-only so its user
        // field is always overwritten with the fresh auto.
        if (b.p_peak_ch[c] < 0) b.p_peak_ch[c] = p_auto;
        if (b.q_begin_ch[c] < 0) b.q_begin_ch[c] = q_auto;
        b.r_peak_ch[c] = r_auto;
        if (b.s_end_ch[c] < 0) b.s_end_ch[c] = s_auto;
        if (b.t_begin_ch[c] < 0) b.t_begin_ch[c] = tp_auto;
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