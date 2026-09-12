/*feature_marks.cpp -- implementations for FeatureMarks.
See feature_marks.hpp for the public interface*/

#include "feature_marks.hpp"
#include "template_marking_gui\template_marking_bin_io.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstdio>
#include <cstdint>
#include <limits>
#include <numeric>
#include <vector>
#include <functional>
#include "template_anchoring\anchor_fit.hpp"
#include "template_anchoring\landmark_admissibility.hpp"
#include "subsample_refine.hpp"
#include "ppg_derivative.hpp"
#include "ppg_dicrotic.hpp"

namespace {

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
} // anonymous


// =========================================================================
// Private helpers
// =========================================================================

// QRS polarity from the KNOWN R column: positive if the R sample sits above
// the trace baseline (median), negative otherwise. Replaces the old r_peak()
// geometric guesser -- callers now pass the true R column (r_col).
bool FeatureMarks::qrs_positive_at(const std::vector<double>& v, int r_idx) {
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

// Search spans, in seconds. Q-peak within 50 ms before R (mirror of
// S_PEAK_WIN_S); Q-onset within 50 ms before the Q peak.
namespace { constexpr double Q_PEAK_WIN_S = 0.040; constexpr double Q_ONSET_WIN_S = 0.050; }
namespace { constexpr double S_PEAK_WIN_S = 0.100; constexpr double J_POINT_WIN_S = 0.100; }
namespace { constexpr double Q_MIN_DEPTH = 0.0075; }

double FeatureMarks::sample_at(const std::vector<double>& v, double p) {
    const int n = static_cast<int>(v.size());
    const double NaND = std::numeric_limits<double>::quiet_NaN();
    if (n == 0 || !std::isfinite(p) || p < 0.0 || p > n - 1) return NaND;
    const int i = static_cast<int>(std::floor(p));
    const double f = p - static_cast<double>(i);
    if (f == 0.0) return v[i];
    const double a = v[i], b = v[std::min(n - 1, i + 1)];
    if (std::isnan(a) || std::isnan(b)) return NaND;   // a gap stays a gap
    return a + f * (b - a);
}

// SUB-SAMPLE POSITIONS ARE THE ONLY POSITIONS.
//
// compute_q_peak, compute_s_peak and compute_t_peak each ran a subsample_refine
// fit and then lround'ed the result back to an integer column, so the
// refinement cost its compute and returned the seed the detector had already
// found. All three now return the refined double and NOTHING here rounds it.
// Where a caller needs the trace's value at a landmark it interpolates
// (sample_at below) rather than indexing a rounded column: a landmark at
// 104.37 has an amplitude, and it is not ecg[104].
double FeatureMarks::compute_q_peak(const std::vector<double>& ecg, int r_idx, double fs)
{
    const int N = static_cast<int>(ecg.size());
    if (r_idx <= 0 || r_idx >= N) return -1;
    auto cl = [&](int i) { return std::clamp(i, 0, N - 1); };

    const bool is_positive = FeatureMarks::qrs_positive_at(ecg, r_idx);
    std::vector<double> u = ecg;
    if (!is_positive) for (auto& x : u) x = -x;

    // Search range: [R - 50 ms, R]. Q is the trough just before R, so on the
    // upright copy it is the minimum.
    const int qp_lo = cl(r_idx - static_cast<int>(std::lround(Q_PEAK_WIN_S * fs)));
    const int qp_hi = cl(r_idx);
    int qSeed = qp_hi;
    double qv = std::numeric_limits<double>::infinity();
    for (int i = qp_lo; i <= qp_hi; ++i)
        if (!std::isnan(u[i]) && u[i] < qv) { qv = u[i]; qSeed = i; }

    // A REAL Q is a strict interior minimum -- the signal turns back up before R
    // -- and deep enough to be a wave rather than the meeting point of the P
    // tail and the R upstroke.
    if (qSeed <= qp_lo || qSeed >= qp_hi) return -1;

    // PQ isoelectric reference: median of [R - 100 ms, R - 40 ms] on the same
    // upright copy (the window seed_all uses for B_iso).
    double b_iso;
    {
        const int a = cl(r_idx - static_cast<int>(std::lround(0.100 * fs)));
        const int b = cl(r_idx - static_cast<int>(std::lround(0.040 * fs)));
        std::vector<double> s;
        for (int i = std::min(a, b); i <= std::max(a, b); ++i)
            if (!std::isnan(u[i])) s.push_back(u[i]);
        if (s.empty()) return -1;
        std::nth_element(s.begin(), s.begin() + s.size() / 2, s.end());
        b_iso = s[s.size() / 2];
    }
    const double rAmp = u[r_idx] - b_iso;
    if (rAmp > 0.0 && (b_iso - u[qSeed]) < Q_MIN_DEPTH * rAmp) return -1;   // too shallow to be a Q

    // Symmetric extremum: Gaussian-weighted quadratic, sigma = 4.
    return std::clamp(subsample_refine::symmetricExtremum(u, qSeed, 4.0),
        0.0, static_cast<double>(N - 1));
}

// S = the first opposite-polarity trough after R: walk right from R tracking
// the opposite-polarity extreme (minimum if the QRS is upright, maximum if
// inverted) and stop once the trace has clearly turned back. Rate-aware search
// window (~0.12 s), so it needs no s_end bound -- the single S-trough source,
// used for the s_end detection and for |R|+|S| normalization.
double FeatureMarks::compute_s_peak(const std::vector<double>& ecg, int r_idx, double fs)
{
    const int N = static_cast<int>(ecg.size());
    if (r_idx < 0 || r_idx >= N - 1)
        return static_cast<double>(std::clamp(r_idx + 1, 0, std::max(0, N - 1)));
    auto cl = [&](int i) { return std::clamp(i, 0, N - 1); };

    const bool is_positive = FeatureMarks::qrs_positive_at(ecg, r_idx);
    std::vector<double> u = ecg;
    if (!is_positive) for (auto& x : u) x = -x;

    // Search range: [R, R + 50 ms]. S is the trough just after R, so on the
    // upright copy it is the minimum. Mirror of compute_q_peak's range.
    const int sp_lo = cl(r_idx);
    const int sp_hi = cl(r_idx + static_cast<int>(std::lround(S_PEAK_WIN_S * fs)));

    int sSeed = sp_lo;
    double sv = std::numeric_limits<double>::infinity();
    for (int i = sp_lo; i <= sp_hi; ++i)
        if (!std::isnan(u[i]) && u[i] < sv) { sv = u[i]; sSeed = i; }

    // Symmetric extremum: Gaussian-weighted quadratic, sigma = 4.
    return std::clamp(subsample_refine::symmetricExtremum(u, sSeed, 4.0),
        0.0, static_cast<double>(N - 1));
}

// -------------------------------------------------------------------------
// Reactive ECG X-glyphs: each is auto-computed but tracks the user's movable
// markers live. Windows are +/-0.05 s around the relevant user marker.
// -------------------------------------------------------------------------

// T peak = max value between the user's T-begin and T-end markers.
// T-peak: single canonical finder. T-peak is an AMPLITUDE landmark (a peak),
// so per spec it is placed at the fitted peak and refined as an ASYMMETRIC
// EXTREMUM -- cubic fit on Gaussian-weighted samples with the derivative
// solved analytically, sigma = 15 -- NOT via the transition-upsample path used
// for onsets/offsets. Bracketed by the S-end and T-end bars (the T-wave's position varies
// with heart rate, so it can't be derived from R alone). Both the movable bar
// (T_PEAK locator) and the auto glyph call this; the sigma=15 refinement is
// folded in here so neither caller re-refines.
// Parameters named for what they ARE -- the two bars bracketing the search --
// not for a t_begin landmark that no longer exists and never held a value. The
// order does not matter; they are min/max'd below.
double FeatureMarks::compute_t_peak(const std::vector<double>& v, double bracketSEnd, double bracketTEnd) {
    // T peak = the sample furthest from the bracket-endpoint baseline, over the
    // window between the S-end and T-end bars, clamped into the trace's FINITE
    // region. Because the window always contains signal, the extremum always
    // resolves: no failure branch, no fallback, and the only -1 is a trace with
    // no finite sample at all (an empty trace, not a T-detection case).
    const int N = static_cast<int>(v.size());
    if (N < 1) return -1.0;

    int fFin = -1, lFin = -1;
    for (int i = 0; i < N; ++i)      if (!std::isnan(v[i])) { fFin = i; break; }
    for (int i = N - 1; i >= 0; --i) if (!std::isnan(v[i])) { lFin = i; break; }
    if (fFin < 0) return -1.0;   // entirely NaN: no data

    // Bracket ordered and clamped to [fFin, lFin]. A -1 bracket (absent bar)
    // clamps to the finite edge rather than aborting; a collapsed bracket
    // becomes one sample.
    const double loD = std::min(bracketSEnd, bracketTEnd);
    const double hiD = std::max(bracketSEnd, bracketTEnd);
    int lo = std::clamp(static_cast<int>(std::ceil(loD)), fFin, lFin);
    int hi = std::clamp(static_cast<int>(std::floor(hiD)), fFin, lFin);
    if (hi < lo) hi = lo;

    int a = lo, b = hi;
    while (a <= b && std::isnan(v[a])) ++a;
    while (b >= a && std::isnan(v[b])) --b;
    if (b < a) return static_cast<double>(lo);   // window was a NaN gap

    // Baseline = mean of the (finite) window ends; both sit at the T's feet, so
    // their mean is a local isoelectric reference. The peak is the sample
    // furthest from it -- polarity-independent, so an inverted T works too.
    const double B = 0.5 * (v[a] + v[b]);
    int best = a; double bd = -1.0;
    for (int i = a; i <= b; ++i)
        if (!std::isnan(v[i]) && std::abs(v[i] - B) > bd) { bd = std::abs(v[i] - B); best = i; }

    // Refine on a copy oriented so the peak is a maximum (sigma = 15), then keep
    // it inside the finite window. Coarse column stands if refinement is
    // non-finite -- the detected position, not a default.
    std::vector<double> u = v;
    if (v[best] < B) for (auto& x : u) x = -x;
    const double p = subsample_refine::asymmetricExtremum(u, best, 15.0);
    return std::isfinite(p)
        ? std::clamp(p, static_cast<double>(lo), static_cast<double>(hi))
        : static_cast<double>(best);
}

double FeatureMarks::compute_p_peak(const std::vector<double>& v, double loIn, double hiIn, double fs)
{
    // P peak = the largest deviation from the endpoint chord over the pre-QRS
    // window. The window is clamped into the trace's FINITE region, so it
    // always contains signal and the argmax always resolves. There is no
    // failure branch and no fallback: the only -1 is a trace with no finite
    // sample at all, which is not a P-detection case but an empty trace.
    const int N = static_cast<int>(v.size());
    if (N < 1 || fs <= 0.0) return -1.0;

    int fFin = -1, lFin = -1;
    for (int i = 0; i < N; ++i)      if (!std::isnan(v[i])) { fFin = i; break; }
    for (int i = N - 1; i >= 0; --i) if (!std::isnan(v[i])) { lFin = i; break; }
    if (fFin < 0) return -1.0;   // entirely NaN: no data

    // Bracket ordered, 20 ms trimmed off the top so the Q upstroke cannot win,
    // then clamped to [fFin, lFin]. A -1 bracket (absent bar) clamps to the
    // finite edge rather than aborting; a collapsed bracket becomes one sample.
    const double loD = std::min(loIn, hiIn);
    const double hiD = std::max(loIn, hiIn) - 0.020 * fs;
    int lo = std::clamp(static_cast<int>(std::ceil(loD)), fFin, lFin);
    int hi = std::clamp(static_cast<int>(std::floor(hiD)), fFin, lFin);
    if (hi < lo) hi = lo;

    int a = lo, b = hi;
    while (a <= b && std::isnan(v[a])) ++a;
    while (b >= a && std::isnan(v[b])) --b;
    if (b < a) return static_cast<double>(fFin);   // window was a NaN gap

    const double yL = v[a], yR = v[b];
    const double invSpan = (b > a) ? 1.0 / static_cast<double>(b - a) : 0.0;

    int best = a; double bd = -std::numeric_limits<double>::infinity();
    for (int i = a; i <= b; ++i) {
        if (std::isnan(v[i])) continue;
        const double dev = v[i] - (yL + (yR - yL) * (i - a) * invSpan);
        if (dev > bd) { bd = dev; best = i; }
    }

    // sigma = 12, the P-peak sigma. Coarse argmax stands if refinement is
    // non-finite -- that is still the detected column, not a fallback default.
    const double p = subsample_refine::asymmetricExtremum(v, best, 12.0);
    return std::isfinite(p)
        ? std::clamp(p, static_cast<double>(fFin), static_cast<double>(lFin))
        : static_cast<double>(best);
}

// SINGLE canonical S-peak + J-point finder (spec steps 1-5). Used by BOTH the
// movable bar (J_POINT locator) and the auto glyph. Steps:
//   1) S-peak search range = [R, R + 50 ms].
//   2) S-peak = Gaussian-weighted quadratic (symmetricExtremum), sigma = 4.
//   3) The offset found here is the J-point.
//   4) J-point search range = [S-peak, S-peak + 50 ms].
//   5) 40-sample window, 4x cubic upsample, fit-and-select (piecewise-linear /
//      sigmoid / fractional-poly by BIC) via transitionAnchor; anchor placed
//      at the offset fractional level.
double FeatureMarks::compute_j_point(const std::vector<double>& v, double fs, int r_col) {
    const int N = static_cast<int>(v.size());
    if (r_col < 0 || r_col >= N - 1 || N < 4) return -1.0;
    auto cl = [&](int i) { return std::clamp(i, 0, N - 1); };
    auto cld = [&](double d) { return std::clamp(d, 0.0, static_cast<double>(N - 1)); };
    auto ms = [&](double s) { return static_cast<int>(std::lround(s * fs)); };

    const bool is_positive = FeatureMarks::qrs_positive_at(v, r_col);
    std::vector<double> u = v;
    if (!is_positive) for (auto& x : u) x = -x;

    // Steps 1-2: the S peak, from the single canonical finder -- argmin over
    // [R, R + 50 ms] refined by the sigma=4 Gaussian-weighted quadratic. Same
    // value computeEcgFeatures uses for |R|+|S|.
    const double sPeakD = FeatureMarks::compute_s_peak(v, r_col, fs);

    // Step 4: J-point search range = [S-peak, S-peak + 50 ms]. transitionAnchor
    // needs an integer seed and integer bounds, so the SEARCH is on columns --
    // but sPeakD is carried through to the return so a window that is too short
    // still reports the refined S position rather than a rounded one.
    const int sPeak = cl(static_cast<int>(std::floor(sPeakD)));
    const int lo = cl(sPeak);
    const int hi = cl(sPeak + ms(J_POINT_WIN_S));
    if (hi - lo < 4) return cld(sPeakD);

    // Baseline reference = right edge of the J-point window (recovered ST
    // level); extremum = the S trough at lo.
    const double baseline = u[hi];

    // 4x cubic-upsample transition fit-and-select. Offset anchor at the
    // recovered-baseline end: fraction 0.10, matching every other onset and
    // offset in this file (was 0.08, which had no stated reason).
    return cld(subsample_refine::transitionAnchor(u, sPeak, 0.10, 40, baseline, lo, hi));
}

// Q-onset finder, structured exactly like compute_j_point (spec steps 1-5,
// mirrored for an onset):
//   1) Q-peak search range = [R - 50 ms, R], via compute_q_peak.
//   2) Q-peak = Gaussian-weighted quadratic (symmetricExtremum), sigma = 4.
//   4) Q-onset search range = [Q-peak - 50 ms, Q-peak].
//   5) 40-sample window, 4x cubic upsample, fit-and-select via transitionAnchor;
//      onset anchor placed at 0-20% of the onset (fraction 0.10, baseline side).
//
// qPeakIn is an OPTIONAL SEED, not a requirement: pass a Q peak already in hand
// and this skips step 1, so detect_template_landmarks does not search for the
// same landmark twice. -1 (the default) means "find it yourself".
//
// `measured` (optional) reports whether the returned onset came from a real Q
// wave. FALSE means one of two things, and neither is comparable to a normal
// Q-onset:
//   - the monophasic-R fallback ran (no Q trough at all), which is a
//     slope-if walk down the R upstroke -- a different measurement, not
//     a worse version of the same one; or
//   - a Q peak was found but the fit window was too short, so the onset fell
//     back to the peak itself.
// The glyph layer draws a hollow circle for these rather than an X. When it is
// false because there was no Q trough, compute_q_peak has already returned -1
// for the same reason, so the Q-peak mark is simply absent -- the hollow onset
// and the missing peak are two readings of one return value and cannot disagree.
double FeatureMarks::compute_q_onset(const std::vector<double>& v, double fs, int r_idx, double qPeakIn, bool* measured) {
    if (measured) *measured = false;   // set true only on the fit path below

    const int N = static_cast<int>(v.size());
    if (r_idx <= 0 || r_idx >= N || N < 4) return -1.0;
    auto cl = [&](int i) { return std::clamp(i, 0, N - 1); };
    auto cld = [&](double d) { return std::clamp(d, 0.0, static_cast<double>(N - 1)); };
    auto ms = [&](double s) { return static_cast<int>(std::lround(s * fs)); };

    const bool is_positive = FeatureMarks::qrs_positive_at(v, r_idx);
    std::vector<double> u = v;
    if (!is_positive) for (auto& x : u) x = -x;

    // Steps 1-2: the Q peak, from the single canonical finder unless the caller
    // supplied one. -1 means no strict interior trough deep enough to be a Q
    // wave (monophasic R) -> R-upstroke fallback below.
    const double qPeakD = (qPeakIn >= 0.0) ? qPeakIn
        : FeatureMarks::compute_q_peak(v, r_idx, fs);

    if (qPeakD >= 0.0) {
        // transitionAnchor needs an integer seed and integer bounds, so the
        // SEARCH runs on columns. floor, not lround: hi is this same column, so
        // rounding up would put the seed outside its own window when the
        // sub-sample peak sits just below an integer.
        const int qPeak = cl(static_cast<int>(std::floor(qPeakD)));

        // Step 4: Q-onset search range = [Q-peak - 50 ms, Q-peak].
        const int lo = cl(qPeak - ms(Q_ONSET_WIN_S));
        const int hi = cl(qPeak);
        if (hi - lo >= 4) {
            // Baseline = left edge (PQ baseline, the level the onset rises FROM).
            const double baseline = u[lo];
            // Step 5: 40-sample 4x cubic-upsample transition fit-and-select.
            // Onset anchor at 0-20% (fraction 0.10, baseline side).
            if (measured) *measured = true;
            return cld(subsample_refine::transitionAnchor(u, qPeak, 0.10, 40, baseline, lo, hi));
        }
        // Window too short to fit: the onset is placed at the peak. A position,
        // but not an onset measurement -- hence measured stays false.
        return cld(qPeakD);
    }

    // No Q trough (monophasic R): R-upstroke onset. Walk left from R down the
    // steep rise to where the slope flattens to < 10% of the peak upstroke
    // slope. Scan window is 50 ms before R. measured stays false throughout.
    const int win = std::max(2, ms(0.050));
    const int scanLo = std::max(1, r_idx - win);
    // Record WHERE the steepest upstroke is, not just how steep. The old code
    // started the threshold walk at r_idx itself, but r_idx is the apex, where
    // the slope is ~0 -- so `s < thresh` fired on the first iteration and the
    // function returned r_idx every time. Walk left from the steepest point.
    int maxSlopePos = r_idx;
    double maxSlope = 0.0;
    for (int i = r_idx; i > scanLo; --i) {
        const double s = u[i] - u[i - 1];
        if (s > maxSlope) { maxSlope = s; maxSlopePos = i; }
    }
    if (maxSlope <= 0.0) return cl(r_idx);   // no rise: anchor on R
    const double thresh = 0.10 * maxSlope;
    for (int i = maxSlopePos; i > scanLo; --i) {
        const double s = u[i] - u[i - 1];
        if (s < thresh) return cl(i);         // slope flattened: onset
    }
    return cl(scanLo);
}

// T-offset: where the T wave returns to baseline. Window is bounded by T-begin
// -- [T-begin + 100 ms, T-begin + 200 ms] -- so the landmark chain runs
// J-point -> T-begin -> T-end, each bounding the next. Previously the window was
// seed +- 100 ms around a downhill-walk T-end estimate: the search for
// T-end was bounded by a guess at T-end.
double FeatureMarks::compute_t_end(const std::vector<double>& v, double fs, int r_col, double j_point) {
    const int N = static_cast<int>(v.size());
    if (N < 4 || r_col < 0 || r_col >= N) return -1.0;
    auto cl = [&](int i) { return std::clamp(i, 0, N - 1); };
    auto cld = [&](double d) { return std::clamp(d, 0.0, static_cast<double>(N - 1)); };

    // Window: [T-begin + 100 ms, T-begin + 350 ms].
    const int lo0 = cl(static_cast<int>(std::lround(j_point + 0.100 * fs)));
    const int hi = cl(static_cast<int>(std::lround(j_point + 0.350 * fs)));
    if (hi <= lo0 + 3) return cld(j_point);

    // B = post-T baseline at the right edge. E = the extremum in the window, by
    // |distance| so an inverted T behaves the same.
    const double B = v[hi];
    int ePos = lo0; double bestDist = 0.0;
    for (int i = lo0; i <= hi; ++i) {
        if (std::isnan(v[i])) continue;
        const double d = std::abs(v[i] - B);
        if (d > bestDist) { bestDist = d; ePos = i; }
    }
    // Start the fit AT the extremum: the T crosses a near-baseline level twice,
    // once rising and once recovering, and anchorAtFraction takes the FIRST
    // crossing. Excluding the upslope leaves only the recovery.
    const int lo = std::min(ePos, hi - 3);
    const double E = v[ePos];

    // Near-baseline target (f = 0.10 == 90% recovered), matching every other
    // onset and offset here.
    if (std::isnan(B) || std::isnan(E)) return -1.0;
    auto fit = anchor_fit::selectAnchorModel(v, lo, hi);
    const double af = anchor_fit::anchorAtFraction(fit, lo, hi, B, E, 0.02);
    if (!std::isfinite(af)) return -1.0;
    const int seed = std::clamp(static_cast<int>(std::round(af)), 0, N - 1);
    // NaN MUST NOT REACH cld. transitionAnchor returns NaN when its window sits
    // in a NaN run -- its own guard is `if (isnan(v)) v = local.front()`, which
    // does nothing when front is also NaN -- and std::clamp(NaN, lo, hi) returns
    // NaN, because neither of its comparisons fires. The NaN then travelled to
    // detect_template_landmarks, where keep() folded it to -1, and the T-end bar
    // and glyph both silently vanished with nothing anywhere reporting a
    // failure. Absent is reported explicitly instead.
    const double te = subsample_refine::transitionAnchor(v, seed, 0.02, 40, B, lo, hi);
    if (!std::isfinite(te)) return -1.0;
    return cld(te);
}

// P begin: anchor-fit elbow on the ascending onset of the P wave. The
// human-editable P-onset marker. Mirrors compute_q_onset (f = 0.10).
// P begin: anchor-fit elbow on the ascending onset of the P wave. The
// human-editable P-onset marker. Range is [start of ECG, P-peak] -- the onset
// precedes the peak, and cannot sit left of the first ECG sample. No other
// clamps.
double FeatureMarks::compute_p_begin(const std::vector<double>& v, double fs, int r_idx, double pPeakIn) {
    const int N = static_cast<int>(v.size());
    if (N < 1 || fs <= 0.0) return -1.0;

    int fFin = -1;
    for (int i = 0; i < N; ++i) if (!std::isnan(v[i])) { fFin = i; break; }
    if (fFin < 0) return -1.0;   // entirely NaN: no data

    // Upper bound is the P peak. compute_p_peak never returns -1 for a finite
    // trace, so neither does the onset. Bracket: from the finite start to
    // Q-onset (or R when Q is unresolved -- the QRS anchor, not an arbitrary
    // offset). No fallback path.
    double pPeak = pPeakIn;
    if (!(pPeak >= 0.0)) {
        const double q = FeatureMarks::compute_q_onset(v, fs, r_idx);
        const double hi = (q >= 0.0) ? q : static_cast<double>(r_idx);
        pPeak = FeatureMarks::compute_p_peak(v, static_cast<double>(fFin), hi, fs);
    }

    // Onset = the elbow on the upstroke leading into the peak, over
    // [start, peak]. The transition fit's own column stands if refinement is
    // non-finite -- the detected position, not a default.
    const int pUser = std::clamp(static_cast<int>(std::lround(pPeak)), fFin, N - 1);
    // 100 ms window back from the peak (was 50 ms): a broad P's onset can sit
    // well before the peak, and a short window took its baseline B already on
    // the rising P, which pushed the 10% crossing up the upstroke.
    const int w = std::max(1, static_cast<int>(std::lround(0.100 * fs)));
    const bool is_positive = FeatureMarks::qrs_positive_at(v, r_idx);
    std::vector<double> u = v;
    if (!is_positive) for (auto& x : u) x = -x;
    const int lo = std::max(fFin, pUser - w);
    const int hi = std::min(pUser, N - 1);
    const double B = u[std::clamp(lo, 0, N - 1)];
    const double pb = subsample_refine::transitionAnchor(u, pUser, 0.10, 40, B, lo, hi);
    return std::clamp(std::isfinite(pb) ? pb : static_cast<double>(pUser),
        static_cast<double>(fFin), pPeak);
}

AnchorLocator make_anchor_locator(AnchorType type, int r_col, double fs) {
    switch (type) {
    case AnchorType::R_PEAK:  return [r_col](const std::vector<double>&) {
        return static_cast<double>(r_col); };
    case AnchorType::Q_ONSET:
        // Single canonical finder: detect_q_onset supplies the seed;
        // compute_q_onset does the Q-peak search + sigma=4 refine + I-3
        // transitionAnchor refinement, and folds in the R-upstroke fallback
        // for monophasic-R (no-Q) beats. Same shape as J_POINT below.
        return [r_col, fs](const std::vector<double>& b) {
            return FeatureMarks::compute_q_onset(b, fs, r_col);
            };
    case AnchorType::J_POINT:
        // One call: compute_j_point re-derives the S peak itself and places the
        // anchor via transitionAnchor. No longer rounded: AnchorLocator is the
        // sub-sample form now.
        return [r_col, fs](const std::vector<double>& b) {
            return FeatureMarks::compute_j_point(b, fs, r_col);
            };
    case AnchorType::P_ONSET:
        // compute_p_begin detects the P peak itself, so this is one call.
        return [r_col, fs](const std::vector<double>& b) {
            return FeatureMarks::compute_p_begin(b, fs, r_col);
            };
    }
    return [](const std::vector<double>&) { return -1.0; };
}


// ============================================================================
// Reactive glyph bundles -- the single definition of every bracketed glyph.
// Pure functions of a trace + bracketing bar positions; nothing is cached, so
// the GUI (per repaint) and the writers (per row) cannot disagree. The caller
// decides WHICH bars to bracket with: the *_auto fields for an autodetect
// column, the user MarkerSet for a user column.
// ============================================================================

FeatureMarks::ReactiveEcg FeatureMarks::reactive_ecg(const std::vector<double>& ecg, double p_begin, double q_onset, double s_end, double t_end, double sampleRate)
{
    ReactiveEcg r;
    if (static_cast<int>(ecg.size()) < 3) return r;
    r.p_peak = compute_p_peak(ecg, p_begin, q_onset, sampleRate);
    r.t_peak = compute_t_peak(ecg, s_end, t_end);
    return r;
}

FeatureMarks::ReactivePpg FeatureMarks::reactive_ppg(const std::vector<double>& ppg, double onset, double peak, double dicrotic, double end)
{
    ReactivePpg r;
    if (static_cast<int>(ppg.size()) < 3) return r;

    // THE BRACKETS ROUND, THE RESULTS DO NOT. The bars arrive sub-sample, but
    // amplitude_crossing and crossing_at_level take an integer search grid and
    // return an interpolated position inside it, so a bracket half a sample
    // either way does not move the answer. Rounding here, once, keeps those
    // bodies integer-indexed instead of pushing fractional indices into every
    // v[i] and loop bound inside them.
    auto win = [](double x) {
        return (x < 0.0) ? -1 : static_cast<int>(std::lround(x));
        };
    const int iOnset = win(onset), iPeak = win(peak), iEnd = win(end);

    if (iOnset >= 0 && iPeak > iOnset)
        r.t50 = amplitude_crossing(ppg, iOnset, iPeak, 0.50);
    if (iPeak >= 0 && iEnd > iPeak) {
        r.t80 = amplitude_crossing(ppg, iPeak, iEnd, 0.80);
        // T80_rise / PW80 at t80's OWN absolute level (see detect_ppg_
        // fiducials for the rationale): upslope crossing of the same value.
        // The two amplitudes come from the UNROUNDED bar positions -- a level
        // is a measurement, not a window, and sample_at interpolates.
        if (iOnset >= 0 && iPeak > iOnset) {
            const double vp = sample_at(ppg, peak);
            const double ve = sample_at(ppg, end);
            if (std::isfinite(vp) && std::isfinite(ve)) {
                const double target = vp + 0.80 * (ve - vp);
                const double xr = crossing_at_level(ppg, iOnset, iPeak, target);
                if (xr >= 0.0) {
                    r.t80_rise = xr;
                    if (r.t80 >= 0.0 && r.t80 > r.t80_rise) r.pw80 = r.t80 - r.t80_rise;
                }
            }
        }
    }
    //diastolic peak - highest first dir between dn foot - 20ms
    // `dicrotic` IS NOT READ by this function, and never was: every output
    // above is bracketed by onset, peak and end. Dragging the dicrotic bar
    // therefore changes no reactive value. Left in the signature so the call
    // matches BinPlotWidget::reactiveGlyphs; worth revisiting separately.
    if (iPeak >= 0 && iEnd > iPeak) {
        const double t80 = amplitude_crossing(ppg, iPeak, iEnd, 0.80);
        r.peak2 = detect_ppg_peak2(ppg, iPeak, t80, iEnd);
    }
    return r;

}

// Position at which the trace reaches `frac` of the a-to-b amplitude.
//
// Returns a SUB-SAMPLE position. This used to return the nearest column, which
// quantised T80 and P50 to whole samples -- at 256 Hz that is a 3.9 ms floor on
// an interval whose whole clinical value is that small differences in it
// separate groups (the T80 entropy result in Section 6.3). The bracketing
// columns are found as before, then the crossing is interpolated between them.
double FeatureMarks::amplitude_crossing(const std::vector<double>& v, int a, int b, double frac) {
    const int N = static_cast<int>(v.size());
    if (a < 0 || b < 0 || b <= a || b >= N) return -1.0;
    const double va = v[a], vb = v[b];
    if (std::isnan(va) || std::isnan(vb)) return -1.0;
    const double target = va + frac * (vb - va);

    // First pair of adjacent finite samples that straddles the target.
    const bool rising = (vb >= va);
    for (int i = a + 1; i <= b; ++i) {
        if (std::isnan(v[i]) || std::isnan(v[i - 1])) continue;
        const bool crossed = rising ? (v[i] >= target && v[i - 1] < target)
            : (v[i] <= target && v[i - 1] > target);
        if (!crossed) continue;
        const double den = v[i] - v[i - 1];
        const double f = (den != 0.0) ? (target - v[i - 1]) / den : 0.0;
        return (i - 1) + std::clamp(f, 0.0, 1.0);
    }
    // No straddle (monotone miss, or a plateau at the target): fall back to the
    // closest column, as before, so the caller still gets a position.
    int best = a; double bestDiff = std::numeric_limits<double>::infinity();
    for (int i = a; i <= b; ++i) {
        if (std::isnan(v[i])) continue;
        const double d = std::abs(v[i] - target);
        if (d < bestDiff) { bestDiff = d; best = i; }
    }
    return static_cast<double>(best);
}

double FeatureMarks::crossing_at_level(const std::vector<double>& v, int a, int b, double target) {
    const int N = static_cast<int>(v.size());
    if (a < 0 || b < 0 || b <= a || b >= N) return -1.0;
    const double va = v[a], vb = v[b];
    if (std::isnan(va) || std::isnan(vb) || std::isnan(target)) return -1.0;
    const bool rising = (vb >= va);
    for (int i = a + 1; i <= b; ++i) {
        if (std::isnan(v[i]) || std::isnan(v[i - 1])) continue;
        const bool crossed = rising ? (v[i] >= target && v[i - 1] < target)
            : (v[i] <= target && v[i - 1] > target);
        if (!crossed) continue;
        const double den = v[i] - v[i - 1];
        const double f = (den != 0.0) ? (target - v[i - 1]) / den : 0.0;
        return (i - 1) + std::clamp(f, 0.0, 1.0);
    }
    return -1.0;
}

double FeatureMarks::first_crossing(const std::vector<double>& v, int a, int b, double frac) {
    const int N = static_cast<int>(v.size());
    if (a < 0 || b <= a || b >= N) return -1.0;
    const double va = v[a], vb = v[b];
    if (std::isnan(va) || std::isnan(vb)) return -1;
    const double target = va + frac * (vb - va);
    const bool rising = (vb >= va);
    for (int i = a + 1; i <= b; ++i) {
        if (std::isnan(v[i]) || std::isnan(v[i - 1])) continue;
        const bool crossed = rising ? (v[i] >= target && v[i - 1] < target)
            : (v[i] <= target && v[i - 1] > target);
        if (!crossed) continue;
        const double den = v[i] - v[i - 1];
        const double f = (den != 0.0) ? (target - v[i - 1]) / den : 0.0;
        // No clamp on f: the original returned lround((i-1)+f) with f
        // unclamped, and clamping here changed the result on beats where the
        // straddle produced f slightly outside [0,1] -- enough to move a
        // handful of up50 columns and shift the PPG shared width by one.
        return (i - 1) + f;   // interpolated, not rounded
    }
    return -1;
}

int FeatureMarks::trough_in(const std::vector<double>& v, int lo, int hi) {
    const int N = static_cast<int>(v.size());
    lo = std::max(0, lo);
    hi = std::min(hi, N - 1);
    int best = -1;
    double bestV = std::numeric_limits<double>::infinity();
    for (int i = lo; i <= hi; ++i)
        if (!std::isnan(v[i]) && v[i] < bestV) { bestV = v[i]; best = i; }
    return best;
}

// Interior maximum of the first derivative: d[i] greater than both neighbours.
// NOT the global maximum -- on a monotonic decay that is always the last step,
// so peak2 landed a few samples off the foot and drew underneath the end glyph.
// A diastolic shoulder is a LOCAL flattening; the end of the decay is not one.
// Returns -1 when diastole has no interior flattening, which on a featureless
// pulse is the honest answer.
double FeatureMarks::steepest_slope_in(const std::vector<double>& v, int lo, int hi) {
    const int N = static_cast<int>(v.size());
    lo = std::max(0, lo);
    hi = std::min(hi, N - 1);
    if (hi - lo < 4) return -1.0;

    auto d = [&](int i) { return v[i + 1] - v[i]; };   // step i, valid i in [lo, hi-1]

    int best = -1; double bestD = -std::numeric_limits<double>::infinity();
    for (int i = lo + 1; i < hi - 1; ++i) {
        if (std::isnan(v[i - 1]) || std::isnan(v[i]) ||
            std::isnan(v[i + 1]) || std::isnan(v[i + 2])) continue;
        const double dm = d(i - 1), d0 = d(i), dp = d(i + 1);
        if (d0 > dm && d0 >= dp && d0 > bestD) { bestD = d0; best = i; }
    }
    if (best < 0) return -1.0;

    // Sub-sample: parabola through the three derivative samples around best.
    const double dm = d(best - 1), d0 = d(best), dp = d(best + 1);
    const double den = dm - 2.0 * d0 + dp;
    if (std::abs(den) > 1e-12) {
        const double off = 0.5 * (dm - dp) / den;
        if (std::abs(off) <= 1.0) return static_cast<double>(best) + off;
    }
    return static_cast<double>(best);
}

FeatureMarks::PpgFiducials FeatureMarks::detect_ppg_fiducials(const std::vector<double>& v, int W, double ppgRate, double heightMeters)
{
    PpgFiducials g;
    const int N = static_cast<int>(v.size());
    if (N < 3) { fprintf(stderr, "[ppg] bail: N=%d\n", N); return g; }
    const int Wc = std::clamp(W, 2, N);   // visible window; nothing is ever placed past Wc-1
    auto cl = [&](int x) { return std::clamp(x, 0, Wc - 1); };
    // Fractional clamp, and column floor/ceil for the helpers that still need
    // an integer search grid. Positions themselves stay fractional.
    auto cld = [&](double x) { return std::clamp(x, 0.0, static_cast<double>(Wc - 1)); };
    auto iFloor = [&](double x) {
        return std::clamp(static_cast<int>(std::floor(x)), 0, Wc - 1);
        };
    auto iCeil = [&](double x) {
        return std::clamp(static_cast<int>(std::ceil(x)), 0, Wc - 1);
        };

    // ---- Systolic peak: the FIRST upstroke in the visible window, then
 // symmetric extremum (sigma = 8).
 //
 // No R-derived gate. [...comment block stays as-is...]

 // Skip any leading NaN run: the template's first samples sit before the
 // first R (construction pads by `pad` seconds), so a partial pulse there
 // can win the slope gate and drag every landmark onto the left edge.
    int lo0 = 0;
    while (lo0 < Wc && std::isnan(v[lo0])) ++lo0;

    int pkSeed = FeatureMarks::detect_ppg_upstroke_peak(v, lo0, Wc);
    if (pkSeed < 0) {
        // No detectable upstroke - use argmax for ppg peak
        int best = -1; double bv = -std::numeric_limits<double>::infinity();
        for (int i = lo0; i < Wc; ++i)
            if (!std::isnan(v[i]) && v[i] > bv) { bv = v[i]; best = i; }
        pkSeed = best;
    }
    if (pkSeed < 0) return g;              // all-NaN window: nothing to mark
    g.peak = cld(subsample_refine::symmetricExtremum(v, pkSeed, 8.0));

    // A systolic peak with no room for a foot before it is a head fragment,
    // not a pulse. Bail rather than pile every landmark at sample 0.
    if (g.peak < 3) { return g; }

    //systolic foot: asymmetric extremum (sigma = 8) on the rising shoulder
    auto refine_foot = [&](int seed) {
        return cld(subsample_refine::asymmetricExtremum(v, seed, 8.0));
        };

    //systolic foot end of cycle: asymmetric extremum (sigma = 8) after the peak
    auto refine_end = [&](int seed) {
        return cld(subsample_refine::asymmetricExtremum(v, seed, 8.0));
        };

    // Systolic foot: the trough before the peak.
    {
        const int seed = trough_in(v, 0, iFloor(g.peak) - 1);
        g.onset = refine_foot(seed >= 0 ? seed : 0);
    }

    // Pulse end: the trough after the peak.
    {
        const int seed = trough_in(v, std::min(iCeil(g.peak) + 1, Wc - 1), Wc - 1);
        g.end = refine_end(seed >= 0 ? seed : Wc - 1);
    }

    // Dicrotic notch (placeholder tier for now).
    g.dicrotic = cld(g.peak + 0.12 * ppgRate);
    g.notch_found = false;
    g.dn_tier = 0;
    g.dn_confidence = 0.0;

    // ---- T80 / T50: amplitude crossings (the same helper the GUI's reactive
    // T80/T50 glyphs call, so the two can't disagree). ----------------------
    g.t80 = amplitude_crossing(v, iFloor(g.peak), iCeil(g.end), 0.80);
    if (g.t80 < 0) g.t80 = cld(0.5 * (g.peak + g.end));

    // Diastolic peak, bounded on the right by t80 -- so it must come AFTER t80
    // is known. In the block above it read g.t80 while g.t80 was still default,
    // so it always took the fraction-of-span fallback and never saw the bound.
    g.peak2 = detect_ppg_peak2(v, iFloor(g.dicrotic), g.t80, iFloor(g.end));
    // T80_rise: the UPSLOPE (onset->peak) point at the SAME absolute
    // amplitude t80 sits at -- i.e. the 80%-downslope level, measured on the
    // way up, NOT an 80% of onset->peak crossing (which would be a different
    // level). Width pw80 = t80 (downslope) - t80_rise (upslope) at that one
    // shared level. Uses the peak/end anchors amplitude_crossing used for
    // t80, so the level is identical by construction.
    {
        const double vp = sample_at(v, g.peak);
        const double ve = sample_at(v, g.end);
        if (std::isfinite(vp) && std::isfinite(ve)) {
            g.t80_rise_y = vp + 0.80 * (ve - vp);   // == t80's own amplitude
            const double xr = crossing_at_level(v, iFloor(g.onset), iCeil(g.peak), g.t80_rise_y);
            if (xr >= 0.0) {
                g.t80_rise = xr;
                if (g.t80 >= 0.0 && g.t80 > g.t80_rise) g.pw80 = g.t80 - g.t80_rise;
            }
        }
    }
    g.t50 = amplitude_crossing(v, iFloor(g.onset), iCeil(g.peak), 0.50);
    if (g.t50 < 0) g.t50 = cld(0.5 * (g.onset + g.peak));

    // ---- Derivative fiducials: VPG u/v/w, APG a-f, JPG p1/p2. One call;
    // the definitions cross-reference each other, so detect() owns the
    // resolution order.
    {
        const auto d = ppg_deriv::buildDerivatives(v, ppgRate);
        const auto fd = ppg_deriv::detect(d, iFloor(g.onset), iFloor(g.peak),
            iFloor(g.dicrotic), iFloor(g.peak2), Wc);
        g.u = fd.u;  g.v = fd.v;  g.w = fd.w;
        g.a = fd.a;  g.b = fd.b;  g.c = fd.c;
        g.d = fd.d;  g.e = fd.e;  g.f = fd.f;
        g.p1 = fd.p1;  g.p2 = fd.p2;

        // Derived indices from the same derivatives + resolved points. RI reads
        // the original pulse amplitude at p1/p2; SI needs height (NaN => NaN).
        const auto ix = ppg_deriv::computeIndices(v, d, fd, ppgRate, heightMeters);
        g.ba = ix.ba;  g.ca = ix.ca;  g.da = ix.da;  g.ea = ix.ea;  g.fa = ix.fa;
        g.agi = ix.agi;  g.ri = ix.ri;  g.si = ix.si;
        g.foundMask = ix.foundMask;
    }

    return g;
}



// =========================================================================
// Movable (auto-detected seeds)
// =========================================================================

// P-end: walk forward from the P peak until the signal recovers to within
// 10% of a post-P baseline estimate (mirrors detect_s_end's recovery-walk,
// just anchored on P instead of S). Needed for the PQ segment (spec: "end
// of P to immediately before Q onset"), which is a genuinely different
// landmark from P-onset -- no prior detector existed for it.
int FeatureMarks::detect_p_end(const std::vector<double>& ecg_signal, int r_idx, double fs,
    double pPeakIn) {
    const int N = static_cast<int>(ecg_signal.size());
    const bool is_positive = FeatureMarks::qrs_positive_at(ecg_signal, r_idx);
    if (r_idx < 0 || r_idx >= N)
        return std::clamp(r_idx, 0, std::max(0, N - 1));

    std::vector<double> upright = ecg_signal;
    if (!is_positive) for (auto& x : upright) x = -x;

    const double p_idx = (pPeakIn >= 0.0) ? pPeakIn
        : compute_p_peak(ecg_signal, 0.0, compute_q_onset(ecg_signal, fs, r_idx), fs);
    if (p_idx < 0.0 || p_idx >= static_cast<double>(N - 1))
        return p_idx + 1.0;

    // Post-P baseline: a short window just after the peak (P is much
    // shorter than T, so this window is smaller than detect_s_end's).
    const int pb_lo = std::min(p_idx + 20, static_cast<double>(N - 1));
    const int pb_hi = std::min(p_idx + 50, static_cast<double>(N));
    double baseline = upright[p_idx];
    if (pb_hi - pb_lo >= 5) {
        std::vector<double> w(upright.begin() + pb_lo, upright.begin() + pb_hi);
        std::nth_element(w.begin(), w.begin() + w.size() / 2, w.end());
        baseline = w[w.size() / 2];
    }

    const double p_val = upright[p_idx];
    const double depth = baseline - p_val;
    if (depth <= 0.0)
        return std::clamp(p_idx + 15, 0.0, static_cast<double>(N - 1));

    const double target = p_val + 0.90 * depth;
    const int hi = std::min(p_idx + 60, static_cast<double>(N));
    for (int i = p_idx + 1; i < hi; ++i) {
        if (std::isnan(upright[i])) continue;
        if (upright[i] >= target) return i;
    }
    return std::clamp(hi - 1, 0, N - 1);
}




// -------------------------------------------------------------------------
// PPG detectors
// -------------------------------------------------------------------------

// Systolic peak from the upstroke. See the note at the declaration for why
// argmax is wrong here.
//
// Steps, none of which compare amplitudes of competing peaks:
//   1. FIRST significant upstroke. The largest slope only sets a scale; the
//      anchor is the first slope run reaching kSlopeFrac of it. Anchoring on
//      the STEEPEST rise instead fails the same way argmax does -- with a
//      strong reflected wave its upstroke is the steeper one.
//   2. Peak = first point after the anchor where the slope stops being
//      positive, i.e. the first local maximum. A taller later peak is
//      unreachable by construction: the walk stops at the first apex.
//   3. The run must PERSIST (>= h samples). Noise clears any slope gate for a
//      sample or two; a systolic upstroke holds for tens of ms.
//
// The derivative window h scales with the range: a FIXED window does not work
// across heart rates, because at 50 bpm the rise is spread over ~3x the samples
// of a 110 bpm rise so its per-sample slope is ~3x smaller while fixed-window
// noise is unchanged.
int FeatureMarks::detect_ppg_upstroke_peak(const std::vector<double>& v,
    int lo, int hi)
{
    const int n = static_cast<int>(v.size());
    if (hi <= 0 || hi > n) hi = n;
    lo = std::max(0, lo);
    if (hi - lo < 5) return -1;

    const int h = std::max(3, (hi - lo) / 50);

    // Smoothed central-difference derivative; NaN-safe.
    std::vector<double> d(n, std::numeric_limits<double>::quiet_NaN());
    for (int i = std::max(lo, h); i + h < hi; ++i) {
        if (std::isnan(v[i - h]) || std::isnan(v[i + h])) continue;
        d[i] = (v[i + h] - v[i - h]) / (2.0 * h);
    }

    double maxSlope = 0.0;
    for (int i = lo; i < hi; ++i)
        if (!std::isnan(d[i]) && d[i] > maxSlope) maxSlope = d[i];
    if (maxSlope <= 0.0) return -1;                 // flat / no rise
    const double gate = 0.25 * maxSlope;
    const int minRun = std::max(2, h);

    int anchor = -1;
    for (int i = lo; i < hi; ++i) {
        if (std::isnan(d[i]) || d[i] < gate) continue;
        int j = i, bestJ = i, held = 0;
        double bestD = d[i];
        while (j < hi && (std::isnan(d[j]) || d[j] >= gate)) {
            if (!std::isnan(d[j])) {
                ++held;
                if (d[j] > bestD) { bestD = d[j]; bestJ = j; }
            }
            ++j;
        }
        if (held >= minRun) { anchor = bestJ; break; }
        i = j;                                      // too brief: keep looking
    }
    if (anchor < 0) return -1;

    int pk = -1;
    for (int i = anchor; i + 1 < hi; ++i) {
        if (std::isnan(d[i])) continue;
        if (d[i] <= 0.0) { pk = i; break; }
    }
    if (pk < 0) pk = hi - 1;                        // apex at/past the boundary

    // The smoothed derivative crosses zero slightly past the true apex (it
    // averages over +/-h), so settle onto the local maximum sample.
    {
        const int wlo = std::max(lo, pk - h - 1);
        const int whi = std::min(hi, pk + h + 2);
        int bestI = pk; double bestV = -std::numeric_limits<double>::infinity();
        for (int i = wlo; i < whi; ++i)
            if (!std::isnan(v[i]) && v[i] > bestV) { bestV = v[i]; bestI = i; }
        pk = bestI;
    }
    return pk;
}

int FeatureMarks::detect_ppg_onset(const std::vector<double>& pulse) {
    const int N = static_cast<int>(pulse.size());
    if (N < 2) return 0;

    // Both steps delegate: upstroke peak, then the shared trough primitive.
    // No local loops, so this cannot drift from detect_ppg_fiducials' onset.
    const int peak = detect_ppg_upstroke_peak(pulse);
    if (peak < 0) return std::min(5, N - 1);
    const int idx = trough_in(pulse, 0, peak);
    if (idx <= 0) return std::min(5, N - 1);
    return idx;
}

// PPG systolic peak, sub-sample refined. Coarse seed = FIRST peak via the
// upstroke (was argmax, i.e. the tallest peak); refined with a Gaussian-
// weighted quadratic (symmetric extremum, sigma = 8, the per-landmark PPG-peak
// sigma). Returns a FLOAT position so downstream fiducials that key off the
// peak (onset/t80/dicrotic/end brackets) inherit the sub-sample peak.
double FeatureMarks::detect_ppg_peak(const std::vector<double>& pulse) {
    if (pulse.empty()) return 0.0;
    const int seed = detect_ppg_upstroke_peak(pulse);
    if (seed < 0) return 0.0;
    return subsample_refine::symmetricExtremum(pulse, seed, 8.0);
}

int FeatureMarks::detect_ppg_end(const std::vector<double>& pulse) {
    const int N = static_cast<int>(pulse.size());
    if (N < 4) return std::max(0, N - 1);

    // Same two shared primitives as detect_ppg_fiducials' end.
    const int peak = detect_ppg_upstroke_peak(pulse);
    if (peak < 0 || peak >= N - 2) return std::max(0, N - 1);
    const int end = trough_in(pulse, peak + 1, N - 1);
    return (end >= 0) ? end : std::max(0, N - 1);
}

int FeatureMarks::detect_ppg_dicrotic(const std::vector<double>& pulse, int peak) {
    const int N = static_cast<int>(pulse.size());
    peak = std::clamp(peak, 0, std::max(0, N - 1));
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

double FeatureMarks::detect_ppg_peak2(const std::vector<double>& v, int sysPeak, double t80, int end)
{
    //highest first derivative between systolic peak and t80
    const int N = static_cast<int>(v.size());
    if (sysPeak < 0 || end <= sysPeak) return -1.0;
    const int lo = sysPeak + std::max(1, (end - sysPeak) / 20);
    int hi;
    if (t80 >= 0.0)
        hi = std::clamp(static_cast<int>(std::floor(t80)), lo + 1, std::min(end - 1, N - 1));
    else
        hi = std::clamp(lo + static_cast<int>(std::lround(0.4 * (end - sysPeak))),
            lo + 1, std::min(end - 1, N - 1));
    if (hi - lo < 3) return -1.0;
    return steepest_slope_in(v, lo, hi);
}

void FeatureMarks::seed_all(TemplateBin& b, double sampleRate, double ppgRate, AnchorType anchor,
    double heightMeters) {
    // Per-anchor ECG user markers are seeded into this anchor's set.
    // (No bin-wide marker handle. Landmarks are per (lead, slot, anchor) now --
    // see TemplateBin::slotMarks -- so the set is fetched inside the per-channel
    // loop below, where the lead is known.)

    // ---- PPG ------------------------------------------------------------
    if (b.ppgTemplate.empty()) {
        b.bad_ppg = 2;
        b.ppg_onset = b.ppg_t50 = b.ppg_t80 = b.ppg_peak = -1;
        b.ppg_dicrotic = b.ppg_peak2 = b.ppg_end = -1;
        b.ppg_onset_auto = b.ppg_t50_auto = b.ppg_t80_auto = b.ppg_peak_auto = -1;
        b.ppg_dicrotic_auto = b.ppg_peak2_auto = b.ppg_end_auto = -1;
        b.ppg_u_auto = b.ppg_v_auto = b.ppg_w_auto = -1;
        b.ppg_a_auto = b.ppg_b_auto = b.ppg_c_auto = -1;
        b.ppg_d_auto = b.ppg_e_auto = b.ppg_f_auto = -1;
        b.ppg_p1_auto = b.ppg_p2_auto = -1;
        b.ppg_ba_auto = b.ppg_ca_auto = b.ppg_da_auto = NAN;
        b.ppg_ea_auto = b.ppg_fa_auto = NAN;
        b.ppg_agi_auto = b.ppg_ri_auto = b.ppg_si_auto = NAN;
        b.ppg_found_mask_auto = 0;
        b.ppg_dn_tier_auto = 3;  b.ppg_dn_confidence_auto = 0.0;
    }
    else if (b.bad_ppg == 1) {
        b.ppg_onset = b.ppg_t50 = b.ppg_t80 = b.ppg_peak = -1;
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
            const auto pf = FeatureMarks::detect_ppg_fiducials(v, W, ppgRate, heightMeters);

            // THE ONE ROUNDING BOUNDARY.
            //
            // NO ROUNDING ANY MORE. pf carries sub-sample positions throughout
            // and TemplateBin's ppg_*_auto fields are now double, as is the
            // pointer-to-member table that reaches them (PulseAutoGlyph::idx is
            // `double TemplateBin::*`), so the fiducials go in as measured. The
            // rnd() that used to sit here was the last place the pulse
            // detector's precision was thrown away -- and it cost marker
            // display precision only, since the derived indices ba..si and the
            // T80/P50 intervals were always computed on the fractional values.
            b.ppg_peak_auto = pf.peak;
            b.ppg_onset_auto = pf.onset;
            b.ppg_peak2_auto = pf.peak2;
            b.ppg_end_auto = pf.end;
            b.ppg_dicrotic_auto = pf.dicrotic;     b.ppg_dicrotic_found_auto = pf.notch_found;
            b.ppg_t80_auto = pf.t80;
            b.ppg_t50_auto = pf.t50;
            b.ppg_u_auto = pf.u;
            b.ppg_v_auto = pf.v;
            b.ppg_w_auto = pf.w;
            b.ppg_a_auto = pf.a;
            b.ppg_b_auto = pf.b;
            b.ppg_c_auto = pf.c;
            b.ppg_d_auto = pf.d;
            b.ppg_e_auto = pf.e;
            b.ppg_f_auto = pf.f;
            b.ppg_p1_auto = pf.p1;
            b.ppg_p2_auto = pf.p2;
            b.ppg_ba_auto = pf.ba;  b.ppg_ca_auto = pf.ca;  b.ppg_da_auto = pf.da;
            b.ppg_ea_auto = pf.ea;  b.ppg_fa_auto = pf.fa;
            b.ppg_agi_auto = pf.agi;  b.ppg_ri_auto = pf.ri;  b.ppg_si_auto = pf.si;
            b.ppg_found_mask_auto = pf.foundMask;
            b.ppg_dn_tier_auto = pf.dn_tier;  b.ppg_dn_confidence_auto = pf.dn_confidence;
            // Derived doubles (no glyph, not movable): the upslope point at
            // the 80%-downslope level and the width between them. Set here,
            // inside the pf scope, alongside the other pf-derived values.
            b.ppg_t80_rise = pf.t80_rise;
            b.ppg_pw80 = pf.pw80;
        }

        // ---- seed the movable bars once (only when unset) ------------------
        if (b.ppg_onset < 0) b.ppg_onset = b.ppg_onset_auto;
        if (b.ppg_dicrotic < 0) b.ppg_dicrotic = b.ppg_dicrotic_auto;
        if (b.ppg_peak2 < 0) b.ppg_peak2 = b.ppg_peak2_auto;
        if (b.ppg_t80 < 0) b.ppg_t80 = b.ppg_t80_auto;
        if (b.ppg_end < 0) b.ppg_end = b.ppg_end_auto;
        // Auto-only bars: always refreshed.
        b.ppg_peak = b.ppg_peak_auto;
        b.ppg_t50 = b.ppg_t50_auto;
    }

    // ---- ECG (per channel) ---------------------------------------------
    // Ranges defined from landmarks only. Anchor-fit (per spec) finds
    // the precise location. Marker = glyph; user drags from here.
    ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };
    for (int c = 0; c < 3; ++c) {
        const auto& ecg = chs[c]->ecgTemplate_raw;
        if (ecg.empty()) {
            b.bad_r_ch[c] = true;
            // Slot 0's set for this lead. Cleared rather than left alone: an
            // empty channel has no waveform for a landmark to sit on.
            b.slotMarks(c, 0, anchor) = tbank::BankMarkerSet{};
            b.r_peak_ch[c] = -1;
            b.p_peak_auto_ch[c] = b.q_onset_auto_ch[c] = b.r_peak_auto_ch[c] = -1;
            b.s_end_auto_ch[c] = b.t_end_auto_ch[c] = -1;
            b.p_begin_auto_ch[c] = -1;
            b.q_peak_auto_ch[c] = -1;
            b.q_onset_found_auto_ch[c] = false;
            continue;
        }

        // ONE DETECTOR, shared with the bank templates and the archive. This
        // block used to inline the six finder calls and was the only one of the
        // three copies that refined the R anchor first -- which is why the
        // viewer, the bank columns and bin_archive reported three different
        // answers for the same bin. detect_template_landmarks does the
        // refinement and the ordered calls; see feature_marks.hpp.
        const FeatureMarks::TemplateLandmarks lmRaw =
            FeatureMarks::detect_template_landmarks(
                ecg, chs[c]->r_col_raw, sampleRate);

        // Mask to what THIS alignment may report. The finders all run -- one
        // call, cheap -- but a landmark the alignment smeared is reported ABSENT
        // rather than at whatever position the smeared median produced.
        // See landmark_admissibility.hpp.
        // The per-anchor admissibility mask used to live here, to blank the
        // bars this alignment may not report. Bars are no longer seeded in this
        // function (see below), and the glyph *_auto_ch fields are never masked
        // -- a glyph is a measurement reported on every alignment -- so there is
        // nothing left to mask. detect_template_landmarks' result is used raw.

        b.p_peak_auto_ch[c] = lmRaw.p_peak;
        b.q_peak_auto_ch[c] = lmRaw.q_peak;
        b.q_onset_auto_ch[c] = lmRaw.q_onset;
        b.q_onset_found_auto_ch[c] = lmRaw.q_onset_found;
        b.r_peak_auto_ch[c] = lmRaw.r_peak;
        b.s_end_auto_ch[c] = lmRaw.s_end;
        b.t_end_auto_ch[c] = lmRaw.t_end;
        b.p_begin_auto_ch[c] = lmRaw.p_begin;

        // NO PER-ALIGNMENT BAR SEEDING HERE ANY MORE. Each pass used to seed
        // its owned bar from its OWN averaged waveform, so the four bars came
        // from four different traces and could land out of order once
        // frame-shifted together, and disagree with the R-aligned fiducial
        // crosses. All four bars are now seeded once, from the single R-aligned
        // detection, in loadSubject after the seeding loop. This pass still
        // writes the flat *_auto_ch glyph fields (above) and r_peak_ch (below).
        //
        // R falls back to the unrefined column rather than -1: it is the
        // alignment anchor every other landmark is expressed against, so the
        // bin needs SOME R even when refinement could not run.
        b.r_peak_ch[c] = (lmRaw.r_peak >= 0.0)
            ? lmRaw.r_peak
            : std::clamp(static_cast<double>(chs[c]->r_col_raw),
                0.0, static_cast<double>(ecg.size()) - 1.0);
    }


    // ---- Arterial (ABP / ART / ART_PULM) --------------------------------
    const int rFirstArt = (sampleRate > 0.0)
        ? static_cast<int>(std::llround(0.3 * sampleRate)) : 0;
    const int rSecondArt = rFirstArt + ((sampleRate > 0.0)
        ? static_cast<int>(std::llround(0.75 * sampleRate)) : 0);

    auto seedArterial = [&](const std::vector<double>& trace, uint8_t& issue,
        double& onset, double& peak, double& dicrotic, double& peak2, double& end,
        double& onset_auto, double& peak_auto, double& dic_auto, double& p2_auto, double& end_auto)
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
            const FeatureMarks::PpgFiducials pf = FeatureMarks::detect_ppg_fiducials(trace, static_cast<int>(trace.size()), sampleRate, NAN);
            onset_auto = pf.onset; peak_auto = pf.peak; dic_auto = pf.dicrotic;
            p2_auto = pf.peak2; end_auto = pf.end;
            if (onset < 0) onset = pf.onset;
            if (peak < 0) peak = pf.peak;
            if (dicrotic < 0) dicrotic = pf.dicrotic;
            if (peak2 < 0) peak2 = pf.peak2;
            if (end < 0) end = pf.end;
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

FeatureMarks::TemplateLandmarks FeatureMarks::detect_template_landmarks(
    const std::vector<double>& tmplIn, int nominal_r_col, double sampleRate)
{
    TemplateLandmarks out;
    const int n = static_cast<int>(tmplIn.size());
    if (n < 2 || nominal_r_col < 0 || nominal_r_col >= n || sampleRate <= 0.0)
        return out;

    // EDGE MARGIN, ONE PLACE, ALL FINDERS. A steep noise spike at the very
    // first or last drawn sample can look like the largest deviation in the
    // trace, so an argmax-based finder (the P peak especially) locks onto it
    // instead of the real wave. Every finder below already walks in from NaN
    // edges, so blanking a margin at each end of a working copy makes that edge
    // invisible to all of them at once -- no per-finder change. R is refined on
    // the ORIGINAL trace (it is interior by construction and its column must not
    // move); every landmark search then runs on the margined copy. 40 ms is the
    // same order as the finders' own bracket slack.
    const int margin = std::max(1, static_cast<int>(std::lround(0.040 * sampleRate)));
    std::vector<double> tmpl = tmplIn;
    if (2 * margin < n) {
        for (int i = 0; i < margin; ++i)
            tmpl[i] = tmpl[n - 1 - i] = std::numeric_limits<double>::quiet_NaN();
    }

    // R REFINED AGAINST THIS WAVEFORM, not inherited. This is the step the bank
    // and archive paths were missing, and it is the search origin for every
    // finder below -- so omitting it moved every landmark, not just R. Refined
    // on tmplIn (un-margined) so R still anchors even if it sits near an edge.
    const int seed = std::clamp(nominal_r_col, 0, n - 1);
    double r = subsample_refine::symmetricExtremum(tmplIn, seed, 5.0);
    if (std::isnan(r) || r < 0.0 || r > static_cast<double>(n - 1))
        r = static_cast<double>(seed);   // refinement failed; nominal stands
    const int r_anchor = static_cast<int>(r);

    // CALL ORDER IS LOAD-BEARING: the J-point feeds T-offset, the P-peak feeds
    // P-onset, and the Q-peak feeds Q-onset. Each of those parameters exists to
    // stop a second, slightly different search for the same landmark, so
    // reordering these lines changes the answers even though every call looks
    // independent.
    const double j = FeatureMarks::compute_j_point(tmpl, sampleRate, r_anchor);
    const double te = FeatureMarks::compute_t_end(tmpl, sampleRate, r_anchor, j);

    const double qp = FeatureMarks::compute_q_peak(tmpl, r_anchor, sampleRate);
    bool qFound = false;
    const double q = FeatureMarks::compute_q_onset(tmpl, sampleRate, r_anchor, qp, &qFound);
    // P-PEAK UPPER BRACKET, ALWAYS VALID. The P-peak search is an argmax over
    // [0, bracket]; the only way it returns -1 on a template that has a P
    // region is a bracket that is itself invalid. q_onset is the bracket we
    // want, but it is -1 on a monophasic-R (no-Q) beat -- and passing that
    // straight through let a failed Q-onset erase the P wave, which is the bug
    // this fixes. Fall back to a fixed 50 ms before R so the window is never
    // empty and the argmax always runs.
    double pHi = q;
    if (!(pHi >= 0.0)) pHi = r - 0.050 * sampleRate;
    const double pp = FeatureMarks::compute_p_peak(tmpl, 0.0, pHi, sampleRate);
    const double pb = FeatureMarks::compute_p_begin(tmpl, sampleRate, r_anchor, pp);

    // Out-of-range is folded to -1 (absent), NOT clamped to an edge column. A
    // landmark pinned to column 0 is indistinguishable from one genuinely found
    // there, and downstream would integrate over a window that does not exist.
    auto keep = [&](double x) {
        if (std::isnan(x) || x < 0.0 || x > static_cast<double>(n - 1))
            return -1.0;
        return x;
        };

    out.r_peak = r;
    out.q_peak = keep(qp);   // -1 on a monophasic R is the RIGHT answer
    out.q_onset = keep(q);
    out.s_end = keep(j);
    out.t_end = keep(te);
    out.p_peak = keep(pp);   // -1 on a ventricular template is the RIGHT answer
    out.p_begin = keep(pb);
    // A flag that outlives its position would be a lie, so it is anded with the
    // position surviving keep().
    out.q_onset_found = qFound && (out.q_onset >= 0.0);
    out.valid = true;
    return out;
}

void FeatureMarks::seed_bank_template(const std::vector<double>& tmpl, int r_col,
    double sampleRate, AnchorType anchor, tbank::BankMarkerSet& out)
{
    out = tbank::BankMarkerSet{};          // all -1
    const TemplateLandmarks lm =
        FeatureMarks::detect_template_landmarks(tmpl, r_col, sampleRate);
    if (!lm.valid) return;

    // SAME MASK seed_all applies. See landmark_admissibility.hpp.
    const auto msk = landmark_admit::maskFor(anchor);

    // BankMarkerSet is double, so lm's sub-sample positions go in as they are.
    // -1 still means absent.
    //
    // BARS ONLY. p_peak is not seeded here, and no longer exists on
    // BankMarkerSet at all: it is a reactive glyph bracketed by the P-onset and
    // Q-onset bars, so every reader calls FeatureMarks::reactive_ecg on the bar
    // set instead. A detector-sourced copy stored alongside was a second answer
    // that drifted from the X on screen the moment either bracket bar moved.
    if (msk.q_onset) out.q_onset = lm.q_onset;
    if (msk.s_end)   out.s_end = lm.s_end;
    if (msk.t_end)   out.t_end = lm.t_end;
    if (msk.p_begin) out.p_begin = lm.p_begin;
}


void FeatureMarks::seed_pulse_bank_template(const std::vector<double>& tmpl,
    double ppgRate, tbank::BankPulseMarkerSet& out, double heightMeters)
{
    out = tbank::BankPulseMarkerSet{};
    const int W = static_cast<int>(tmpl.size());
    if (W < 3 || ppgRate <= 0.0) return;

    // SAME DETECTOR, THIS TEMPLATE'S OWN WAVEFORM. seed_all runs this on
    // b.ppgTemplate; the viewer draws ppg_bank slots. Those are different
    // pulses, so a foot measured on one is not a minimum on the other -- which
    // is why the foot glyph sat nowhere near a trough.
    const PpgFiducials pf = detect_ppg_fiducials(tmpl, W, ppgRate, heightMeters);

    // BARS AND FROZEN AUTOS ONLY, both sub-sample. peak / peak2 / t50 / t80
    // left BankPulseMarkerSet with the int fields: they are auto-only glyphs
    // that markerAtX never hands out, and the *_auto values below are what
    // BinPlotWidget::overridePulseGlyphs paints. The ix() that used to round
    // all of these is gone -- PpgFiducials is double throughout.
    out.onset_auto = pf.onset;        out.onset = pf.onset;
    out.peak_auto = pf.peak;
    out.dicrotic_auto = pf.dicrotic;  out.dicrotic = pf.dicrotic;
    out.peak2_auto = pf.peak2;
    out.end_auto = pf.end;            out.end = pf.end;
    out.notch_found = pf.notch_found;
}