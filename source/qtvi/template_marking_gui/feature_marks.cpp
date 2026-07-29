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

// P begin: anchor-fit elbow on the ascending onset of the P wave. The
// human-editable P-onset marker. Mirrors compute_q_onset (f = 0.10).
int FeatureMarks::compute_p_begin(const std::vector<double>& v, int pUser, double fs, int r_idx) {
    const int N = static_cast<int>(v.size());
    if (pUser < 0 || pUser >= N || N < 4) return std::clamp(pUser, 0, std::max(0, N - 1));
    const int w = win_005s(fs);
    const bool is_positive = qrs_positive_at(v, r_idx);
    std::vector<double> u = v;
    if (!is_positive) for (auto& x : u) x = -x;
    const int lo = std::max(0, pUser - w);
    const int hi = std::min(N - 1, pUser + w / 4);
    double B = u[lo];
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

AnchorLocator make_anchor_locator(AnchorType type, int r_col, double fs) {
    switch (type) {
    case AnchorType::R_PEAK:  return [r_col](const std::vector<double>&) { return r_col; };
    case AnchorType::Q_ONSET:
        return [r_col, fs](const std::vector<double>& b) {
            const int N = static_cast<int>(b.size());
            if (r_col < 2 || r_col >= N) return r_col;

            // Primary: real Q onset.
            const int q = FeatureMarks::detect_q_begin(b, r_col);
            // detect_q_begin returns r_col-20 (fixed) when it found no Q
            // trough. Treat that exact value as "Q not found" and fall back
            // to the R-upstroke onset instead of the rate-dependent offset.
            if (q != std::max(0, r_col - 20)) return q;

            // R-upstroke onset: walk left from R down the steep rise to
            // where the slope flattens to <10% of the peak upstroke slope.
            // Scan window is 50 ms before R (rate-independent).
            const bool pos = qrs_positive_at(b, r_col);
            auto up = [&](int i) { const double v = b[i]; return pos ? v : -v; };

            const int win = std::max(2, static_cast<int>(std::lround(0.050 * fs)));
            const int scanLo = std::max(1, r_col - win);

            double maxSlope = 0.0;
            for (int i = r_col; i > scanLo; --i) {
                const double s = up(i) - up(i - 1);
                if (s > maxSlope) maxSlope = s;
            }
            if (maxSlope <= 0.0) return r_col;   // no rise: anchor on R

            const double thresh = 0.10 * maxSlope;
            for (int i = r_col; i > scanLo; --i) {
                const double s = up(i) - up(i - 1);
                if (s < thresh) return i;        // slope flattened -> onset
            }
            return scanLo;                        // never flattened in window
            };    case AnchorType::J_POINT: return [r_col, fs](const std::vector<double>& b) { return FeatureMarks::detect_s_end(b, r_col, fs); };
    case AnchorType::P_PEAK:  return [r_col](const std::vector<double>& b) { return FeatureMarks::detect_p_peak(b, r_col); };
    case AnchorType::P_ONSET: // no pure detector: seed at P-peak, refine with compute_p_begin
        return [r_col, fs](const std::vector<double>& b) {
            const int pk = FeatureMarks::detect_p_peak(b, r_col);
            return FeatureMarks::compute_p_begin(b, pk, fs, r_col);
            };
    case AnchorType::T_PEAK:  // no detect_t_peak: bracket via t_begin/t_end
        return [r_col, fs](const std::vector<double>& b) {
            return FeatureMarks::compute_t_peak(b,
                FeatureMarks::detect_t_begin(b, r_col, fs),
                FeatureMarks::detect_t_end(b, r_col, fs));
            };
    }
    return [](const std::vector<double>&) { return -1; };
}

FeatureMarks::EcgGlyphs FeatureMarks::compute_ecg_glyphs(const std::vector<double>& ecg, int p_peak, int q_begin, int s_end, int t_begin, int t_end, double fs)
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

int FeatureMarks::amplitude_crossing(const std::vector<double>& v, int a, int b, double frac) {
    const int N = static_cast<int>(v.size());
    if (a < 0 || b < 0 || b <= a || b >= N) return -1;
    const double va = v[a], vb = v[b];
    if (std::isnan(va) || std::isnan(vb)) return -1;
    const double target = va + frac * (vb - va);
    int best = a; double bestDiff = std::numeric_limits<double>::infinity();
    for (int i = a; i <= b; ++i) {
        if (std::isnan(v[i])) continue;
        const double d = std::abs(v[i] - target);
        if (d < bestDiff) { bestDiff = d; best = i; }
    }
    return best;
}

FeatureMarks::PpgFiducials FeatureMarks::detect_ppg_fiducials(
    const std::vector<double>& v, int constructPeak, int constructOnset,
    int firstR, int W, double fs)
{
    PpgFiducials g;
    const int N = static_cast<int>(v.size());
    if (N < 3) return g;
    const int Wc = std::clamp(W, 2, N);   // visible window; nothing below is ever placed past Wc-1
    auto cl = [&](int x) { return std::clamp(x, 0, Wc - 1); };
    auto pct = [&](double f) { return std::clamp(static_cast<int>(std::lround(f * Wc)), 0, Wc - 1); };
    const int marginSamples = (fs > 0.0)
        ? std::max(1, static_cast<int>(std::lround(0.010 * fs)))
        : 1;

    // ---- PEAK + ONSET: original detection method, INDEPENDENT of end.
    // Kept independent because the display normalization uses onset as the
    // per-sample perfusion-index baseline (v[onset] is the divisor in
    // pulse_norm) -- if onset lands wrong, the entire displayed trace goes
    // off-scale. Decoupling from end means end detection can't corrupt the
    // display no matter what happens downstream.
    if (constructPeak >= 0 && constructPeak < Wc
        && (firstR < 0 || constructPeak >= firstR)) {
        g.peak = constructPeak;
        g.onset = (constructOnset >= 0 && constructOnset <= g.peak)
            ? constructOnset : 0;
    }
    else {
        // Peak = argmax over [firstR, pct(0.40)). Search starts at the
        // ACTUAL first R sample (passed in from the ECG channel's
        // r_col_raw), not a proportional guess -- PPG templates are
        // R-anchored, so peak by definition must sit at or after R1.
        // If firstR wasn't provided (< 0), skip the constraint and search
        // from 0 as a last resort.
        const int firstPulseStart = (firstR >= 0 && firstR < Wc) ? firstR : 0;
        const int firstPulseEnd = std::max(firstPulseStart + 2, pct(0.40));
        int pk = firstPulseStart; double best = -std::numeric_limits<double>::infinity();
        for (int i = firstPulseStart; i < firstPulseEnd; ++i)
            if (!std::isnan(v[i]) && v[i] > best) { best = v[i]; pk = i; }
        g.peak = pk;
        // Onset = argmin over [pct(0.10), peak] -- the trough just before
        // the systolic upstroke.
        int ft = std::min(pct(0.10), pk); best = std::numeric_limits<double>::infinity();
        for (int i = std::min(pct(0.10), pk); i <= pk; ++i)
            if (!std::isnan(v[i]) && v[i] < best) { best = v[i]; ft = i; }
        g.onset = ft;
    }

    // ---- END: curve-fit detection, bounded to start AFTER peak so it
    // can't degenerately land on the pulse's own foot/onset trough (which
    // is very often deeper than the true diastolic end-of-cycle trough).
    //
    // Stage 1 (coarse locate): plain raw-signal argmin over [peak+1, Wc-1]
    // -- O(n), no model fitting. Full-window BIC model selection here was
    // the earlier performance regression (piecewise-linear breakpoint
    // search is O(n^2), running once per bin).
    // Stage 2 (refine): a small +/-10ms window centered on the coarse
    // point -- THAT small window is where the (cheap, since it's tiny)
    // curve fit runs, for sub-sample precision without paying for it over
    // the whole trace. Bounded to [0, Wc-1] throughout, so end can never
    // land in the off-screen tail-overlap region past the visible window.
    {
        const int lo = std::min(g.peak + 1, Wc - 1);
        const int hi = Wc - 1;
        int coarse = hi;
        {
            double best = std::numeric_limits<double>::infinity();
            for (int i = lo; i <= hi; ++i)
                if (!std::isnan(v[i]) && v[i] < best) { best = v[i]; coarse = i; }
        }

        const int rlo = cl(coarse - marginSamples);
        const int rhi = cl(coarse + marginSamples);
        int refined = coarse;
        if (rhi - rlo >= 3) {
            const auto rfit = anchor_fit::selectAnchorModel(v, rlo, rhi);
            const double rt = anchor_fit::anchorAtPeak(rfit, rlo, rhi, /*findMin=*/true);
            refined = cl(static_cast<int>(std::lround(rt)));
        }
        else {
            double best = std::numeric_limits<double>::infinity();
            for (int i = rlo; i <= rhi; ++i)
                if (!std::isnan(v[i]) && v[i] < best) { best = v[i]; refined = i; }
        }

        g.end = refined;
        g.end_found = true;   // the refinement always yields SOME position
    }

    // Scale threshold for the 2nd-derivative extrema tests below.
    // Requiring d2 > kThresh * pulse_amplitude keeps the tests
    // scale-invariant (works whether the raw signal is in mV, ADC
    // counts, or normalized units) and filters out noise-driven
    // near-zero d2 values that were previously being called notches.
    const double pulseAmp = (!std::isnan(v[g.peak]) && !std::isnan(v[g.end]))
        ? std::abs(v[g.peak] - v[g.end])
        : 0.0;
    const double kThresh = 0.001;    // 1% of pulse amplitude per sample^2
    const double curvThresh = kThresh * pulseAmp;

    // Amplitude landmarks: samples where the trace has descended a given
    // fraction of the peak-to-end amplitude range.
    int amp30 = amplitude_crossing(v, g.peak, g.end, 0.30);
    int amp60 = amplitude_crossing(v, g.peak, g.end, 0.60);
    int amp80 = amplitude_crossing(v, g.peak, g.end, 0.80);
    if (amp30 < 0) amp30 = g.peak;
    if (amp60 < 0) amp60 = g.end;
    if (amp80 < 0) amp80 = g.end;
    amp30 = std::clamp(amp30, 0, Wc - 1);
    amp60 = std::clamp(amp60, amp30, Wc - 1);
    amp80 = std::clamp(amp80, amp60, Wc - 1);

    // ---- DICROTIC NOTCH: sample with the steepest (most positive) 2nd
    // derivative in [30%-height-down, 60%-height-down]. The 2nd derivative
    // d2[i] = v[i+1] - 2v[i] + v[i-1] is positive where the curve is
    // concave-up (a "cup" -- physically, a dip). The threshold requires
    // d2 to exceed a small fraction of pulse amplitude, so noise-level
    // fluctuations don't get called a notch. If no sample clears the
    // threshold, draw the O in the middle of the range.
    {
        const int lo = amp30;
        const int hi = std::max(lo, amp60);

        int notch = -1;
        double bestCurv = curvThresh;
        for (int i = lo + 1; i < hi; ++i) {
            if (std::isnan(v[i - 1]) || std::isnan(v[i]) || std::isnan(v[i + 1])) continue;
            const double d2 = v[i + 1] - 2.0 * v[i] + v[i - 1];
            if (d2 > bestCurv) { bestCurv = d2; notch = i; }
        }

        if (notch >= 0) {
            g.dicrotic = std::clamp(notch, lo, hi);
            g.notch_found = true;
        }
        else {
            g.dicrotic = cl((lo + hi) / 2);
            g.notch_found = false;
        }
    }

    // ---- PEAK2 (diastolic peak): sample with the least steep (most
    // negative) 2nd derivative in [DN + 5ms, T80]. Lower bound in x is
    // 5ms past the notch (skip the notch's own curvature-recovery
    // region); upper bound is the 80%-height-down marker (T80). If no
    // sample clears the threshold in magnitude, draw the O in the middle.
    {
        const int fiveMsSamples = (fs > 0.0)
            ? std::max(1, static_cast<int>(std::lround(0.005 * fs)))
            : 1;
        const int lo = std::clamp(g.dicrotic + fiveMsSamples, 0, Wc - 1);
        const int hi = std::max(lo, amp80);

        int p2 = -1;
        double bestCurv = -curvThresh;
        for (int i = lo + 1; i < hi; ++i) {
            if (std::isnan(v[i - 1]) || std::isnan(v[i]) || std::isnan(v[i + 1])) continue;
            const double d2 = v[i + 1] - 2.0 * v[i] + v[i - 1];
            if (d2 < bestCurv) { bestCurv = d2; p2 = i; }
        }

        if (p2 >= 0) {
            g.peak2 = std::clamp(p2, lo, hi);
            g.peak2_found = true;
        }
        else {
            g.peak2 = cl((lo + hi) / 2);
            g.peak2_found = false;
        }
    }

    // ---- T80 / P50: amplitude crossings (shared helper -- same formula
    // the GUI's reactive T80/P50 glyphs use). ------------------------------
    g.t80 = amplitude_crossing(v, g.peak, g.end, 0.80);
    if (g.t80 < 0) g.t80 = cl((g.peak + g.end) / 2);
    g.p50 = amplitude_crossing(v, g.onset, g.peak, 0.50);
    if (g.p50 < 0) g.p50 = cl((g.onset + g.peak) / 2);

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

        // Movable-marker spawn points.
        // Seeds only; the user drags from here.

        // ---- PPG autodetect positions: ONE call, single source of truth
        // for every fiducial. Nothing here recomputes anything independently
        // -- the GUI's frozen glyphs read straight from these same *_auto
        // fields (see BinPlotWidget::captureGlyphSnapshot), so "peak" (etc.)
        // can never mean two different things in two different places.
        {
            // firstR = the actual R position, taken from any available ECG
            // channel's r_col_raw (all three should agree, since ECG and
            // PPG templates are R-anchored to the same sample). -1 if no
            // ECG channel is present, which lets detect_ppg_fiducials fall
            // back to searching from 0.
            int firstR = -1;
            const ChannelTemplateData* cc[3] = { &b.ch1, &b.ch2, &b.ch3 };
            for (const auto* ch : cc) {
                if (!ch->ecgTemplate_raw.empty() && ch->r_col_raw >= 0) {
                    firstR = ch->r_col_raw;
                    break;
                }
            }

            const auto pf = FeatureMarks::detect_ppg_fiducials(
                v, b.ppg_peak_construct, b.ppg_onset_construct, firstR, W, sampleRate);
            b.ppg_peak_auto = pf.peak;
            b.ppg_onset_auto = pf.onset;
            b.ppg_peak2_auto = pf.peak2;           b.ppg_peak2_found_auto = pf.peak2_found;
            b.ppg_end_auto = pf.end;               b.ppg_end_found_auto = pf.end_found;
            b.ppg_dicrotic_auto = pf.dicrotic;     b.ppg_dicrotic_found_auto = pf.notch_found;
            b.ppg_t80_auto = pf.t80;
            b.ppg_p50_auto = pf.p50;
        }

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
        const Band T_BEGIN_BAND = { 0.090,  0.180 };  // T-begin: .a = ST-baseline start (.b unused; hi = T peak)
        const Band P_PEAK_BAND = { -0.220, -0.120 };  // P-wave peak (argmax)
        const Band P_BEGIN_BAND = { -0.320, -0.250 };  // P-wave onset (left of P peak)

        // Amplitude landmark inside a band: argmax |ecg - isoelectric|.
        auto peakInBand = [&](Band bnd, const char* tag) -> int {
            int lo = cl(r_auto + ms(bnd.a));
            int hi = cl(r_auto + ms(bnd.b));
            if (hi < lo) std::swap(lo, hi);
            int p = lo; double bd = -1.0;
            for (int i = lo; i <= hi; ++i)
                if (!std::isnan(ecg[i]) && std::abs(ecg[i] - B_iso) > bd) { bd = std::abs(ecg[i] - B_iso); p = i; }
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
        }
        // S-end: mirror of Q-onset. Cap the LEFT bound at the S peak (first
        // turning point walking right from R = the S extremum) so E is the S,
        // not the R tail. Then follow the spec crossing over [S peak, ST
        // baseline]. To land at the END of the S wave -- the J-point, where the
        // wave has recovered to baseline (the spec's "80-100% of an offset") --
        // the target sits near baseline. In the literal formula L = B + f*(E-B)
        // that is f = 0.10 (== 90% recovered from the extremum), NOT 0.90 which
        // would sit at the trough. Recovery is monotonic S -> baseline, so the
        // forward crossing lands at the J-point.
        int s_auto;
        {
            const int ss = cl(r_auto + ms(0.080));
            int spk = -1;
            for (int i = r_auto + 1; i < ss; ++i) {
                if (std::isnan(ecg[i - 1]) || std::isnan(ecg[i]) || std::isnan(ecg[i + 1])) continue;
                const double dL = ecg[i] - ecg[i - 1];
                const double dR = ecg[i + 1] - ecg[i];
                if (dL != 0.0 && dL * dR < 0.0) { spk = i; break; }   // local extremum (S)
            }
            if (spk < 0) spk = cl(r_auto + ms(0.010));   // monophasic: no S, start just after R

            int lo = spk;                                 // capped at the S peak
            int hi = cl(r_auto + ms(S_END_BAND.b));       // ST-baseline side
            if (hi - lo < 4) hi = cl(lo + 4);

            double E = B_iso; double bd = 0.0;
            for (int i = lo; i <= hi; ++i)
                if (!std::isnan(ecg[i]) && std::abs(ecg[i] - B_iso) > bd) { bd = std::abs(ecg[i] - B_iso); E = ecg[i]; }
            auto fit = anchor_fit::selectAnchorModel(ecg, lo, hi);
            const double anchor = anchor_fit::anchorAtFraction(fit, lo, hi, B_iso, E, 0.10);
            s_auto = cl(static_cast<int>(std::round(anchor)));
        }
        // T peak: largest deflection from baseline in the T-wave region. It
        // splits the T wave into its onset (T-begin) and offset (T-end) sides.
        int t_peak;
        {
            const int tlo = cl(r_auto + ms(0.150));
            const int thi = cl(std::min(visN - 1, r_auto + ms(0.400)));
            t_peak = tlo; double bd = 0.0;
            for (int i = tlo; i <= thi; ++i)
                if (!std::isnan(ecg[i]) && std::abs(ecg[i] - B_iso) > bd) { bd = std::abs(ecg[i] - B_iso); t_peak = i; }
        }

        // T-begin: T-wave foot. Mirror of Q-onset -- cap the RIGHT bound at the
        // T peak so E is the T (the old fixed window ended before the T peak, so
        // E was a tiny ST value and the crossing fell back to the midpoint).
        // Window [ST baseline, T peak], spec crossing with near-baseline target
        // (f=0.10) -> the foot where the T departs baseline.
        int tp_auto;
        {
            int lo = cl(r_auto + ms(T_BEGIN_BAND.a));   // ST-segment baseline
            int hi = t_peak;                             // capped at the T peak
            if (hi - lo < 4) lo = cl(hi - 4);

            double E = B_iso; double bd = 0.0;
            for (int i = lo; i <= hi; ++i)
                if (!std::isnan(ecg[i]) && std::abs(ecg[i] - B_iso) > bd) { bd = std::abs(ecg[i] - B_iso); E = ecg[i]; }
            auto fit = anchor_fit::selectAnchorModel(ecg, lo, hi);
            const double anchor = anchor_fit::anchorAtFraction(fit, lo, hi, B_iso, E, 0.10);
            tp_auto = cl(static_cast<int>(std::round(anchor)));
        }
        // T-end: mirror of S-end. Use the shared T peak as the LEFT bound so E
        // is the T, then the spec crossing over [T peak, post-T baseline] with a
        // near-baseline target so the forward crossing lands where the T returns
        // to baseline = the actual T offset (not up at the T peak). Literal
        // L = B + f*(E-B) with f=0.10 == 90% recovered = "80-100% of an offset".
        int te_auto;
        {
            int lo = t_peak;                                        // shared T peak
            int hi = cl(std::min(visN - 1, r_auto + ms(0.520)));     // post-T baseline
            if (hi - lo < 4) hi = cl(lo + 4);

            // E = the T-peak value. Constrain search to a tight window
            // around t_peak (not the whole [lo, hi]) so that any later
            // noise / U wave that happens to be more distant from B_iso
            // than the T peak itself can't flip E to the wrong extremum
            // -- that flip puts L on the wrong side of baseline, kills
            // the anchor-fit crossing, and drops T-end at the midpoint
            // fallback.
            const int eWin = ms(0.020);   // +/-20 ms around t_peak
            const int eLo = cl(t_peak - eWin);
            const int eHi = cl(std::min(hi, t_peak + eWin));
            double E = B_iso; double bd = 0.0;
            for (int i = eLo; i <= eHi; ++i)
                if (!std::isnan(ecg[i]) && std::abs(ecg[i] - B_iso) > bd) { bd = std::abs(ecg[i] - B_iso); E = ecg[i]; }
            auto fit = anchor_fit::selectAnchorModel(ecg, lo, hi);
            const double anchor = anchor_fit::anchorAtFraction(fit, lo, hi, B_iso, E, 0.10);
            te_auto = cl(static_cast<int>(std::round(anchor)));
        }
        const int p_auto = peakInBand(P_PEAK_BAND, "P-peak");
        // P-begin: P-wave foot. Window [pre-P baseline, P peak], anchor-fit
        // crossing near baseline (f=0.10). Right bound capped at the P peak.
        int pb_auto;
        {
            int lo = cl(r_auto + ms(P_BEGIN_BAND.a));   // pre-P baseline
            int hi = p_auto;                             // capped at the P peak
            if (hi - lo < 4) lo = cl(hi - 4);
            double E = B_iso; double bd = 0.0;
            for (int i = lo; i <= hi; ++i)
                if (!std::isnan(ecg[i]) && std::abs(ecg[i] - B_iso) > bd) { bd = std::abs(ecg[i] - B_iso); E = ecg[i]; }
            auto fit = anchor_fit::selectAnchorModel(ecg, lo, hi);
            const double anchor = anchor_fit::anchorAtFraction(fit, lo, hi, B_iso, E, 0.10);
            pb_auto = cl(static_cast<int>(std::round(anchor)));
        }

        // Auto fields always updated.
        b.p_peak_auto_ch[c] = p_auto;
        b.q_begin_auto_ch[c] = q_auto;
        b.r_peak_auto_ch[c] = r_auto;
        b.s_end_auto_ch[c] = s_auto;
        b.t_begin_auto_ch[c] = tp_auto;
        b.t_end_auto_ch[c] = te_auto;
        b.p_begin_auto_ch[c] = pb_auto;

        // User fields: only seed when unset. R peak is auto-only so its user
        // field is always overwritten with the fresh auto.
        if (b.p_peak_ch[c] < 0) b.p_peak_ch[c] = p_auto;
        if (b.q_begin_ch[c] < 0) b.q_begin_ch[c] = q_auto;
        b.r_peak_ch[c] = r_auto;
        if (b.s_end_ch[c] < 0) b.s_end_ch[c] = s_auto;
        if (b.t_begin_ch[c] < 0) b.t_begin_ch[c] = tp_auto;
        if (b.t_end_ch[c] < 0) b.t_end_ch[c] = te_auto;
        if (b.p_begin_ch[c] < 0) b.p_begin_ch[c] = pb_auto;
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