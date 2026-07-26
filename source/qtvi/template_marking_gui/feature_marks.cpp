//
// feature_marks.cpp -- implementations for FeatureMarks.
// See feature_marks.hpp for the public interface.
//

#include "feature_marks.hpp"
#include "TemplateBinIO.hpp"    // for TemplateBin (used by seed_all)

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>   // TEMP: diagnostics
#include <limits>
#include <numeric>
#include <vector>
#include <functional>
#include "anchor_fit.hpp"

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

// QRS polarity from the KNOWN R column: positive if the R sample sits above
// the trace baseline (median), negative otherwise. Replaces the old r_peak()
// geometric guesser -- callers now pass the true R column (r_col).
static bool qrs_positive_at(const std::vector<double>& v, int r_idx) {
    const int N = static_cast<int>(v.size());
    if (N == 0 || r_idx < 0 || r_idx >= N || std::isnan(v[r_idx])) return true;
    std::vector<double> f;
    f.reserve(N);
    for (double x : v) if (!std::isnan(x)) f.push_back(x);
    if (f.empty()) return true;
    std::nth_element(f.begin(), f.begin() + f.size() / 2, f.end());
    const double baseline = f[f.size() / 2];
    return (v[r_idx] - baseline) >= 0.0;
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

    const bool up = qrs_positive_at(ecg, r_peak_idx);

    int q = q_begin;
    for (int i = q_begin; i <= r_peak_idx; ++i) {
        if (std::isnan(ecg[i]) || std::isnan(ecg[q])) continue;
        if (up ? ecg[i] < ecg[q] : ecg[i] > ecg[q]) q = i;
    }
    return q;
}

// S = the first opposite-polarity trough after R: walk right from R tracking
// the opposite-polarity extreme (minimum if the QRS is upright, maximum if
// inverted) and stop once the trace has clearly turned back. Rate-aware search
// window (~0.12 s), so it needs no s_end bound -- the single S-trough source,
// used for the s_end detection and for |R|+|S| normalization.
int FeatureMarks::compute_s_peak(const std::vector<double>& ecg, int r_idx, double fs)
{
    const int N = static_cast<int>(ecg.size());
    if (r_idx < 0 || r_idx >= N - 1)
        return std::clamp(r_idx + 1, 0, std::max(0, N - 1));

    const bool up = qrs_positive_at(ecg, r_idx);
    const int w = (fs > 0.0)
        ? std::max(4, static_cast<int>(std::lround(0.12 * fs)))
        : 60;
    const int hi = std::min(r_idx + w, N);

    int s = r_idx;
    for (int i = r_idx + 1; i < hi; ++i) {
        if (std::isnan(ecg[i])) continue;
        if (std::isnan(ecg[s])) { s = i; continue; }
        const bool moreExtreme = up ? (ecg[i] < ecg[s]) : (ecg[i] > ecg[s]);
        if (moreExtreme) s = i;
        else if (i > s + 1) break;   // turned back after the trough -> done
    }
    if (s <= r_idx) return std::min(r_idx + 1, N - 1);
    return s;
}

// -------------------------------------------------------------------------
// Reactive ECG X-glyphs: each is auto-computed but tracks the user's movable
// markers live. Windows are +/-0.05 s around the relevant user marker.
// -------------------------------------------------------------------------

int FeatureMarks::compute_p_peak(const std::vector<double>& ecg,
    int p_onset, int q_begin, int r_peak_idx)
{
    const int N = static_cast<int>(ecg.size());
    if (N == 0 || p_onset < 0 || q_begin < 0
        || p_onset >= N || q_begin >= N || p_onset > q_begin)
        return std::clamp(std::max(p_onset, 0), 0, std::max(0, N - 1));

    const bool is_positive = qrs_positive_at(ecg, r_peak_idx);
    std::vector<double> u = ecg;
    if (!is_positive) for (auto& x : u) x = -x;

    int best = p_onset;
    for (int i = p_onset; i <= q_begin; ++i)
        if (!std::isnan(u[i]) && u[i] > u[best]) best = i;
    return best;
}

// T peak = max value between the user's T-begin and T-end markers.
int FeatureMarks::compute_t_peak(const std::vector<double>& v, int tBegin, int tEnd) {
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
int FeatureMarks::compute_r_peak(const std::vector<double>& v, int qBegin, int sEnd) {
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

int FeatureMarks::compute_s_end(const std::vector<double>& v, int sUser, double fs) {
    const int N = static_cast<int>(v.size());
    if (sUser < 0 || N < 3) return sUser;
    const int w = win_005s(fs);
    const int lo = std::max(0, sUser - w);
    const int hi = std::min(N - 1, sUser + w / 4);
    if (hi <= lo) return std::clamp(sUser, 0, N - 1);

    // B = ST baseline: right edge of window (flat region after recovery).
    double B = v[hi];

    // E = extremum (sample furthest from B in the window).
    double E = B;
    double bestDist = 0.0;
    for (int i = lo; i <= hi; ++i) {
        if (std::isnan(v[i])) continue;
        const double d = std::abs(v[i] - B);
        if (d > bestDist) { bestDist = d; E = v[i]; }
    }

    auto fit = anchor_fit::selectAnchorModel(v, lo, hi);
    const double anchor = anchor_fit::anchorAtFraction(fit, lo, hi, B, E, 0.90);
    return std::clamp(static_cast<int>(std::round(anchor)), 0, N - 1);
}

int FeatureMarks::compute_q_onset(const std::vector<double>& v, int qUser, double fs, int r_idx) {
    const int N = static_cast<int>(v.size());
    if (qUser < 0 || qUser >= N || N < 4) return std::clamp(qUser, 0, std::max(0, N - 1));
    const int w = win_005s(fs);

    const bool is_positive = qrs_positive_at(v, r_idx);
    std::vector<double> u = v;
    if (!is_positive) for (auto& x : u) x = -x;

    // Window extends left into PQ baseline, barely right past seed.
    const int lo = std::max(0, qUser - w);
    const int hi = std::min(N - 1, qUser + w / 4);

    // B = PQ baseline: left edge of window.
    double B = u[lo];

    // E = Q trough (min in the window).
    double E = B;
    for (int i = lo; i <= hi; ++i)
        if (!std::isnan(u[i]) && u[i] < E) E = u[i];

    auto fit = anchor_fit::selectAnchorModel(u, lo, hi);
    const double anchor = anchor_fit::anchorAtFraction(fit, lo, hi, B, E, 0.10);
    return std::clamp(static_cast<int>(std::round(anchor)), 0, N - 1);
}

int FeatureMarks::compute_t_end(const std::vector<double>& v, int tEndUser, double fs) {
    const int N = static_cast<int>(v.size());
    if (tEndUser < 0 || N < 4) return tEndUser;
    const int w = win_005s(fs);
    const int lo = std::max(0, tEndUser - 2 * w);
    const int hi = std::min(N - 1, tEndUser + w / 4);

    // B = post-T baseline: right edge of window.
    double B = v[hi];

    // E = extremum (sample furthest from B in the window).
    double E = B;
    double bestDist = 0.0;
    for (int i = lo; i <= hi; ++i) {
        if (std::isnan(v[i])) continue;
        const double d = std::abs(v[i] - B);
        if (d > bestDist) { bestDist = d; E = v[i]; }
    }

    auto fit = anchor_fit::selectAnchorModel(v, lo, hi);
    const double anchor = anchor_fit::anchorAtFraction(fit, lo, hi, B, E, 0.90);
    return std::clamp(static_cast<int>(std::round(anchor)), 0, N - 1);
}

// P onset: anchor-fit on the ascending onset of the P wave. No glyph;
// used only to bound the P-peak search (P peak = max between p_onset
// and q_begin).
int FeatureMarks::compute_p_onset(const std::vector<double>& v, int pUser, double fs, int r_idx) {
    const int N = static_cast<int>(v.size());
    if (pUser < 0 || pUser >= N || N < 4) return std::clamp(pUser, 0, std::max(0, N - 1));
    const int w = win_005s(fs);

    const bool is_positive = qrs_positive_at(v, r_idx);
    std::vector<double> u = v;
    if (!is_positive) for (auto& x : u) x = -x;

    // Window extends left into pre-P baseline, barely right past seed.
    const int lo = std::max(0, pUser - w);
    const int hi = std::min(N - 1, pUser + w / 4);

    // B = pre-P baseline: left edge of window.
    double B = u[lo];

    // E = P peak (max in the window).
    double E = -std::numeric_limits<double>::infinity();
    for (int i = lo; i <= hi; ++i)
        if (!std::isnan(u[i]) && u[i] > E) E = u[i];

    auto fit = anchor_fit::selectAnchorModel(u, lo, hi);
    const double anchor = anchor_fit::anchorAtFraction(fit, lo, hi, B, E, 0.10);
    return std::clamp(static_cast<int>(std::round(anchor)), 0, N - 1);
}

// T begin: anchor-fit on the ascending onset of the T wave.
int FeatureMarks::compute_t_begin(const std::vector<double>& v, int tBeginUser, double fs, int r_idx) {
    const int N = static_cast<int>(v.size());
    if (tBeginUser < 0 || tBeginUser >= N || N < 4) return std::clamp(tBeginUser, 0, std::max(0, N - 1));
    const int w = win_005s(fs);

    const bool is_positive = qrs_positive_at(v, r_idx);
    std::vector<double> u = v;
    if (!is_positive) for (auto& x : u) x = -x;

    // Window extends left into ST baseline, short right past marker.
    const int lo = std::max(0, tBeginUser - w);
    const int hi = std::min(N - 1, tBeginUser + w / 4);

    // B = ST baseline: left edge of window (isoelectric after S-end).
    double B = u[lo];

    // E = T peak (max in window — transition ascends toward it).
    double E = -std::numeric_limits<double>::infinity();
    for (int i = lo; i <= hi; ++i)
        if (!std::isnan(u[i]) && u[i] > E) E = u[i];

    auto fit2 = anchor_fit::selectAnchorModel(u, lo, hi);
    const double anchor2 = anchor_fit::anchorAtFraction(fit2, lo, hi, B, E, 0.10);
    return std::clamp(static_cast<int>(std::round(anchor2)), 0, N - 1);
}

FeatureMarks::EcgGlyphs FeatureMarks::compute_ecg_glyphs(
    const std::vector<double>& ecg,
    int p_peak, int q_begin, int s_end, int t_begin, int t_end, double fs)
{
    EcgGlyphs g;
    const int N = static_cast<int>(ecg.size());
    g.p_peak_glyph = (p_peak >= 0 && p_peak < N) ? p_peak : -1;
    g.q_begin_glyph = (q_begin >= 0 && q_begin < N) ? q_begin : -1;
    g.r_peak_glyph = compute_r_peak(ecg, q_begin, s_end);
    g.s_end_glyph = (s_end >= 0 && s_end < N) ? s_end : -1;
    g.t_peak_glyph = compute_t_peak(ecg, t_begin, t_end);
    g.t_end_glyph = (t_end >= 0 && t_end < N) ? t_end : -1;
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

int FeatureMarks::detect_q_begin(const std::vector<double>& ecg_signal, int r_idx) {
    const int N = static_cast<int>(ecg_signal.size());
    const bool is_positive = qrs_positive_at(ecg_signal, r_idx);

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

int FeatureMarks::detect_p_peak(const std::vector<double>& ecg_signal, int r_idx) {
    const int N = static_cast<int>(ecg_signal.size());
    if (N < 3) return 0;

    const bool is_positive = qrs_positive_at(ecg_signal, r_idx);
    std::vector<double> upright = ecg_signal;
    if (!is_positive) for (auto& x : upright) x = -x;

    // Max value before Q begin.
    const int q_est = detect_q_begin(ecg_signal, r_idx);
    const int hi = std::max(1, std::min(q_est, N));
    int best = -1; double bv = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < hi; ++i)
        if (!std::isnan(upright[i]) && upright[i] > bv) { bv = upright[i]; best = i; }
    return (best >= 0) ? best : 0;
}

int FeatureMarks::detect_t_begin(const std::vector<double>& ecg_signal, int r_idx, double fs) {
    // First local maximum right of S end (the T wave), then walk left to its
    // foot. begin/end bracket that peak.
    const int N = static_cast<int>(ecg_signal.size());
    if (N < 3) return 0;
    const bool is_positive = qrs_positive_at(ecg_signal, r_idx);
    std::vector<double> upright = ecg_signal;
    if (!is_positive) for (auto& x : upright) x = -x;

    const int s_end = detect_s_end(ecg_signal, r_idx, fs);
    const int lo = std::clamp(s_end + 1, 1, N - 2);
    const int hi = std::min(N - 1, N / 2);              // first half only

    int tPeak = -1;
    double tBest = -std::numeric_limits<double>::infinity();
    for (int i = lo; i < hi; ++i) {
        if (std::isnan(upright[i - 1]) || std::isnan(upright[i]) || std::isnan(upright[i + 1])) continue;
        if (upright[i] > upright[i - 1] && upright[i] >= upright[i + 1]) {
            if (upright[i] > tBest) { tBest = upright[i]; tPeak = i; }   // tallest = the T wave, not the first ST bump
        }
    }
    if (tPeak < 0) return std::clamp(s_end + (int)std::lround(0.1 * N), 0, hi);

    int i = tPeak;   // left foot
    while (i - 1 > s_end && !std::isnan(upright[i - 1]) && upright[i - 1] < upright[i]) --i;
    return i;
}

int FeatureMarks::detect_s_end(const std::vector<double>& ecg_signal, int r_idx, double fs) {
    const int N = static_cast<int>(ecg_signal.size());
    const bool is_positive = qrs_positive_at(ecg_signal, r_idx);
    if (r_idx < 0 || r_idx >= N - 1)
        return std::clamp(r_idx + 1, 0, std::max(0, N - 1));

    std::vector<double> upright = ecg_signal;
    if (!is_positive) for (auto& x : upright) x = -x;

    const int s_idx = compute_s_peak(ecg_signal, r_idx, fs);   // single trough source
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

int FeatureMarks::detect_t_end(const std::vector<double>& ecg_signal, int r_idx, double fs) {
    // First local maximum right of S end, then walk right to its foot.
    const int N = static_cast<int>(ecg_signal.size());
    if (N < 3) return std::max(0, N - 1);
    const bool is_positive = qrs_positive_at(ecg_signal, r_idx);
    std::vector<double> upright = ecg_signal;
    if (!is_positive) for (auto& x : upright) x = -x;

    const int s_end = detect_s_end(ecg_signal, r_idx, fs);
    const int lo = std::clamp(s_end + 1, 1, N - 2);
    const int hi = std::min(N - 1, N / 2);

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

        // ---- PPG autodetect positions (all set here) -----------------------
        // Systolic peak + foot: prefer the construction-time fiducials (exact,
        // one pulse); else a first-pulse-bounded search.
        if (b.ppg_peak_construct >= 0 && b.ppg_peak_construct < W) {
            b.ppg_peak_auto = b.ppg_peak_construct;
            b.ppg_onset_auto = (b.ppg_onset_construct >= 0 &&
                b.ppg_onset_construct <= b.ppg_peak_auto)
                ? b.ppg_onset_construct : 0;
        }
        else {
            int pk = 0; double best = -std::numeric_limits<double>::infinity();
            const int firstPulseEnd = std::max(2, pct(0.40));
            for (int i = 0; i < firstPulseEnd && i < W; ++i)
                if (!std::isnan(v[i]) && v[i] > best) { best = v[i]; pk = i; }
            b.ppg_peak_auto = pk;
            int ft = 0; best = std::numeric_limits<double>::infinity();
            for (int i = std::min(pct(0.10), pk); i <= pk; ++i)
                if (!std::isnan(v[i]) && v[i] < best) { best = v[i]; ft = i; }
            b.ppg_onset_auto = ft;
        }

        // End = minimum between the peak and the end of the visible window.
        b.ppg_end_auto = std::min(W - 1, b.ppg_peak_auto + 1);
        {
            double best = std::numeric_limits<double>::infinity();
            for (int i = b.ppg_peak_auto + 1; i < W; ++i)
                if (!std::isnan(v[i]) && v[i] < best) { best = v[i]; b.ppg_end_auto = i; }
        }

        // T80 = 80% of the way down from peak to end.
        b.ppg_t80_auto = std::clamp((b.ppg_peak_auto + b.ppg_end_auto) / 2, 0, W - 1);
        if (b.ppg_end_auto > b.ppg_peak_auto &&
            !std::isnan(v[b.ppg_peak_auto]) && !std::isnan(v[b.ppg_end_auto])) {
            const double target = v[b.ppg_peak_auto]
                + 0.80 * (v[b.ppg_end_auto] - v[b.ppg_peak_auto]);
            double bestDiff = std::numeric_limits<double>::infinity();
            for (int i = b.ppg_peak_auto; i <= b.ppg_end_auto; ++i) {
                if (std::isnan(v[i])) continue;
                const double d = std::abs(v[i] - target);
                if (d < bestDiff) { bestDiff = d; b.ppg_t80_auto = i; }
            }
        }

        // Dicrotic notch = 1/3 from systolic peak to end. P50 = onset..peak
        // midpoint. Peak #2 = 65% of the window.
        b.ppg_dicrotic_auto = (b.ppg_end_auto > b.ppg_peak_auto)
            ? b.ppg_peak_auto + (b.ppg_end_auto - b.ppg_peak_auto) / 3
            : pct(0.60);
        b.ppg_p50_auto = (b.ppg_onset_auto + b.ppg_peak_auto) / 2;
        b.ppg_peak2_auto = pct(0.65);

        // ---- seed the movable bars once (only when unset) ------------------
        if (b.ppg_onset < 0) b.ppg_onset = b.ppg_onset_auto;
        if (b.ppg_dicrotic < 0) b.ppg_dicrotic = b.ppg_dicrotic_auto;
        if (b.ppg_peak2 < 0) b.ppg_peak2 = b.ppg_peak2_auto;
        if (b.ppg_t80 < 0) b.ppg_t80 = b.ppg_t80_auto;
        if (b.ppg_end < 0) b.ppg_end = b.ppg_end_auto;
        // Auto-only bars: always refreshed.
        b.ppg_peak = b.ppg_peak_auto;
        b.ppg_p50 = b.ppg_p50_auto;
    }

    // ---- ECG (per channel) ---------------------------------------------
    // Ranges defined from landmarks only. Anchor-fit (per spec) finds
    // the precise location. Marker = glyph; user drags from here.
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
        auto ms = [&](double sec) { return static_cast<int>(std::lround(sec * sampleRate)); };

        // R is the deterministic anchor column from alignment.
        const int r_auto = cl(chs[c]->r_col_raw);

        // ---- Isoelectric vertical reference (PQ zero) --------------------
        // Median of the signal in a fixed PQ window before R (40-100 ms).
        double B_iso;
        {
            const int a = cl(r_auto - ms(0.100));
            const int b = cl(r_auto - ms(0.040));
            std::vector<double> s;
            for (int i = std::min(a, b); i <= std::max(a, b); ++i)
                if (!std::isnan(ecg[i])) s.push_back(ecg[i]);
            if (s.empty()) B_iso = ecg[cl(r_auto)];
            else { std::nth_element(s.begin(), s.begin() + s.size() / 2, s.end()); B_iso = s[s.size() / 2]; }
        }

        // ---- Landmark search bands, in SECONDS relative to R -------------
        //  *** TUNABLE: replace with your landmark-range table. ***
        //  {start, end} seconds from the R peak; negative = before R.
        //  Anchor-fit runs inside the band; B = B_iso, E = band extremum,
        //  f is the literal spec fraction (onset 0.10, offset 0.90).
        struct Band { double a, b; };
        const Band Q_ONSET_BAND = { -0.060, -0.010 };  // QRS onset
        const Band S_END_BAND = { 0.010,  0.090 };  // J-point / S recovery
        const Band T_BEGIN_BAND = { 0.080,  0.180 };  // T-wave foot
        const Band T_END_BAND = { 0.200,  0.400 };  // T-wave offset
        const Band P_PEAK_BAND = { -0.220, -0.120 };  // P-wave peak (argmax)

        // Anchor-fit inside a band: B = isoelectric, E = band extremum, literal f.
        auto fitInBand = [&](Band bnd, double f, const char* tag) -> int {
            int lo = cl(r_auto + ms(bnd.a));
            int hi = cl(r_auto + ms(bnd.b));
            if (hi - lo < 4) hi = cl(lo + 4);
            double E = B_iso; double bd = 0.0;
            for (int i = lo; i <= hi; ++i)
                if (!std::isnan(ecg[i]) && std::abs(ecg[i] - B_iso) > bd) { bd = std::abs(ecg[i] - B_iso); E = ecg[i]; }
            auto fit = anchor_fit::selectAnchorModel(ecg, lo, hi);
            // ---- TEMP DIAGNOSTIC (remove when done) ----
            const double L = B_iso + f * (E - B_iso);
            bool crossed = false;
            for (int t = lo; t < hi; ++t) {
                const double y0 = fit.f((double)t), y1 = fit.f((double)(t + 1));
                if (std::abs(y0 - L) < 1e-12 || (y0 - L) * (y1 - L) < 0.0) { crossed = true; break; }
            }
            const double anchor = anchor_fit::anchorAtFraction(fit, lo, hi, B_iso, E, f);
            const int ai = cl(static_cast<int>(std::round(anchor)));
            std::fprintf(stderr,
                "[feat ch%d] %-8s r=%d win=[%d,%d] B_iso=%.4g E=%.4g L=%.4g f(lo)=%.4g -> %s land=%.1f raw=%.4g np=%d\n",
                c, tag, r_auto, lo, hi, B_iso, E, L, fit.f((double)lo),
                crossed ? "crossing" : "MIDPOINT-FALLBACK", anchor, ecg[ai], fit.nparams);
            // ---- end diagnostic ----
            return cl(static_cast<int>(std::round(anchor)));
            };

        // Amplitude landmark inside a band: argmax |ecg - isoelectric|.
        auto peakInBand = [&](Band bnd, const char* tag) -> int {
            int lo = cl(r_auto + ms(bnd.a));
            int hi = cl(r_auto + ms(bnd.b));
            if (hi < lo) std::swap(lo, hi);
            int p = lo; double bd = -1.0;
            for (int i = lo; i <= hi; ++i)
                if (!std::isnan(ecg[i]) && std::abs(ecg[i] - B_iso) > bd) { bd = std::abs(ecg[i] - B_iso); p = i; }
            std::fprintf(stderr, "[feat ch%d] %-8s r=%d win=[%d,%d] B_iso=%.4g -> peak@%d (ecg=%.4g)\n",
                c, tag, r_auto, lo, hi, B_iso, p, ecg[p]);   // TEMP DIAGNOSTIC
            return p;
            };

        // Q-onset: right bound is capped at the Q peak (the initial QRS
        // deflection's extremum = first turning point walking left from R), so
        // the fit window is [PQ baseline, Q peak] and E can never be the R.
        // Then follow the spec crossing (f = 0.10) inside that window.
        int q_auto;
        {
            // Q peak = first local extremum in the ~80 ms left of R.
            const int qs = cl(r_auto - ms(0.080));
            int qpk = -1;
            for (int i = r_auto - 1; i > qs; --i) {
                if (std::isnan(ecg[i - 1]) || std::isnan(ecg[i]) || std::isnan(ecg[i + 1])) continue;
                const double dL = ecg[i] - ecg[i - 1];
                const double dR = ecg[i + 1] - ecg[i];
                if (dL != 0.0 && dL * dR < 0.0) { qpk = i; break; }   // local extremum
            }
            if (qpk < 0) qpk = cl(r_auto - ms(0.010));   // monophasic R: no Q, stop just before R

            int lo = cl(r_auto + ms(Q_ONSET_BAND.a));    // PQ-baseline side
            int hi = qpk;                                 // capped at the Q peak
            if (hi - lo < 4) lo = cl(hi - 4);

            double E = B_iso; double bd = 0.0;
            for (int i = lo; i <= hi; ++i)
                if (!std::isnan(ecg[i]) && std::abs(ecg[i] - B_iso) > bd) { bd = std::abs(ecg[i] - B_iso); E = ecg[i]; }
            auto fit = anchor_fit::selectAnchorModel(ecg, lo, hi);
            const double anchor = anchor_fit::anchorAtFraction(fit, lo, hi, B_iso, E, 0.10);
            q_auto = cl(static_cast<int>(std::round(anchor)));
            std::fprintf(stderr, "[feat ch%d] Q-onset  r=%d Qpk=%d win=[%d,%d] B_iso=%.4g E=%.4g land=%.1f raw=%.4g\n",
                c, r_auto, qpk, lo, hi, B_iso, E, anchor, ecg[q_auto]);   // TEMP DIAGNOSTIC
        }
        // Onsets use f = 0.10, offsets f = 0.90 (literal spec fractions).
        const int s_auto = fitInBand(S_END_BAND, 0.90, "S-end");
        const int tp_auto = fitInBand(T_BEGIN_BAND, 0.10, "T-begin");
        const int te_auto = fitInBand(T_END_BAND, 0.90, "T-end");
        const int p_auto = peakInBand(P_PEAK_BAND, "P-peak");

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