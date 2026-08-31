/*feature_marks.cpp -- implementations for FeatureMarks.
See feature_marks.hpp for the public interface*/

#include "feature_marks.hpp"
#include "template_marking_gui\template_marking_bin_io.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <vector>
#include <functional>
#include "anchor_fit.hpp"
#include "subsample_refine.hpp"
#include "ppg_derivative.hpp"
#include "ppg_dicrotic.hpp"

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
namespace { constexpr double Q_PEAK_WIN_S = 0.050; constexpr double Q_ONSET_WIN_S = 0.050; }
// S-peak lives within 50 ms after R; the J-point within 50 ms after the S-peak.
namespace { constexpr double S_PEAK_WIN_S = 0.050; constexpr double J_POINT_WIN_S = 0.050; }
// Minimum depth below the PQ isoelectric level for a trough to count as a real
// Q wave, in the template's own amplitude units (raw mV). A strict interior
// minimum alone is not enough: where the P wave's decaying tail meets the R
// upstroke there is always a shallow interior minimum, so without a depth test
// a beat with no Q at all is accepted and the Q peak lands on the PQ segment.
//
// NOTE the measured depth is the trough below the PQ median, which runs well
// under the Q wave's nominal amplitude because the R wave's rising shoulder
// fills the notch in -- a 0.10 mV Q measures ~0.085 mV deep. So this threshold
// is a floor for "is there a notch at all", not a Q-amplitude criterion.
namespace { constexpr double Q_MIN_DEPTH_MV = 0.01; }

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
    if (b_iso - u[qSeed] < Q_MIN_DEPTH_MV) return -1;   // too shallow to be a Q

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
// for onsets/offsets. Bracketed by t_begin/t_end (the T-wave's position varies
// with heart rate, so it can't be derived from R alone). Both the movable bar
// (T_PEAK locator) and the auto glyph call this; the sigma=15 refinement is
// folded in here so neither caller re-refines.
double FeatureMarks::compute_t_peak(const std::vector<double>& v, double tBegin, double tEnd) {
    const int N = static_cast<int>(v.size());
    if (!(tBegin >= 0.0) || !(tEnd >= 0.0) || tBegin >= N || tEnd >= N) return -1.0;
    // The bracket bars are themselves sub-sample, so the search runs over the
    // integer columns strictly inside them and the RESULT is clamped to the
    // fractional bracket -- the seed search needs samples, the answer does not.
    const double loD = std::min(tBegin, tEnd), hiD = std::max(tBegin, tEnd);
    const int lo = static_cast<int>(std::ceil(loD));
    const int hi = static_cast<int>(std::floor(hiD));
    if (hi <= lo) return loD;

    // Bracketed by the T-begin/T-end bars, so it tracks them live. T polarity is
    // INDEPENDENT of QRS polarity, so this must not assume an upright T (a plain
    // argmax picks the highest shoulder when the T is inverted). Both bracket
    // ends sit at the T's feet, i.e. on baseline, so their mean is a local
    // isoelectric reference and the peak is the sample furthest from it.
    const double B = 0.5 * (v[lo] + v[hi]);
    int best = lo; double bd = -1.0;
    for (int i = lo; i <= hi; ++i)
        if (!std::isnan(v[i]) && std::abs(v[i] - B) > bd) { bd = std::abs(v[i] - B); best = i; }

    // Refine on a copy oriented so the peak is a maximum (asymmetric extremum:
    // cubic on Gaussian-weighted samples, analytic derivative, sigma = 15), then
    // clamp to the BRACKET rather than the array -- a bracketed landmark must
    // stay bracketed no matter where the refinement drifts.
    std::vector<double> u = v;
    if (v[best] < B) for (auto& x : u) x = -x;
    return std::clamp(subsample_refine::asymmetricExtremum(u, best, 15.0), loD, hiD);
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

// Q-onset finder, structured exactly like compute_s_end (spec steps 1-5,
// mirrored for an onset):
//   1) Q-peak search range = [R - 50 ms, R], via compute_q_peak.
//   2) Q-peak = Gaussian-weighted quadratic (symmetricExtremum), sigma = 4.
//   4) Q-onset search range = [Q-peak - 50 ms, Q-peak].
//   5) 40-sample window, 4x cubic upsample, fit-and-select via transitionAnchor;
//      onset anchor placed at 0-20% of the onset (fraction 0.10, baseline side).
double FeatureMarks::compute_q_onset(const std::vector<double>& v, double fs, int r_idx) {
    const int N = static_cast<int>(v.size());
    if (r_idx <= 0 || r_idx >= N || N < 4) return -1.0;
    auto cl = [&](int i) { return std::clamp(i, 0, N - 1); };
    auto cld = [&](double d) { return std::clamp(d, 0.0, static_cast<double>(N - 1)); };
    auto ms = [&](double s) { return static_cast<int>(std::lround(s * fs)); };

    const bool is_positive = FeatureMarks::qrs_positive_at(v, r_idx);
    std::vector<double> u = v;
    if (!is_positive) for (auto& x : u) x = -x;

    // Steps 1-2: the Q peak, from the single canonical finder. -1 means no
    // strict interior trough (monophasic R, no Q) -> R-upstroke fallback below.
    const double qPeakD = FeatureMarks::compute_q_peak(v, r_idx, fs);
    const int qPeak = (qPeakD >= 0.0) ? cl(static_cast<int>(std::floor(qPeakD))) : -1;

    if (qPeakD >= 0.0) {
        // Step 4: Q-onset search range = [Q-peak - 50 ms, Q-peak].
        const int lo = cl(qPeak - ms(Q_ONSET_WIN_S));
        const int hi = cl(qPeak);
        if (hi - lo >= 4) {
            // Baseline = left edge (PQ baseline, the level the onset rises FROM).
            const double baseline = u[lo];
            // Step 5: 40-sample 4x cubic-upsample transition fit-and-select.
            // Onset anchor at 0-20% (fraction 0.10, baseline side).
            return cld(subsample_refine::transitionAnchor(u, qPeak, 0.10, 40, baseline, lo, hi));
        }
        return cld(qPeakD);
    }

    // No Q trough (monophasic R): R-upstroke onset. Walk left from R down the
    // steep rise to where the slope flattens to < 10% of the peak upstroke
    // slope. Scan window is 50 ms before R.
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
double FeatureMarks::compute_t_end(const std::vector<double>& v, double fs, int r_col,
    double tBeginIn) {
    const int N = static_cast<int>(v.size());
    if (N < 4 || r_col < 0 || r_col >= N) return -1.0;
    auto cl = [&](int i) { return std::clamp(i, 0, N - 1); };
    auto cld = [&](double d) { return std::clamp(d, 0.0, static_cast<double>(N - 1)); };

    const double t_begin = (tBeginIn >= 0.0) ? tBeginIn
        : compute_t_begin(v, fs, r_col);
    if (t_begin < 0.0) return -1.0;

    // Window: [T-begin + 100 ms, T-begin + 350 ms].
    const int lo0 = cl(static_cast<int>(std::lround(t_begin + 0.100 * fs)));
    const int hi = cl(static_cast<int>(std::lround(t_begin + 0.350 * fs)));
    if (hi <= lo0 + 3) return cld(t_begin);

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
    auto fit = anchor_fit::selectAnchorModel(v, lo, hi);
    const int seed = std::clamp(static_cast<int>(std::round(
        anchor_fit::anchorAtFraction(fit, lo, hi, B, E, 0.10))), 0, N - 1);
    return cld(subsample_refine::transitionAnchor(v, seed, 0.10, 40, B, lo, hi));
}

// P begin: anchor-fit elbow on the ascending onset of the P wave. The
// human-editable P-onset marker. Mirrors compute_q_onset (f = 0.10).
double FeatureMarks::compute_p_begin(const std::vector<double>& v, double fs, int r_idx,
    double pPeakIn) {
    const int N = static_cast<int>(v.size());
    if (N < 4) return -1.0;
    // Detect the P peak first, unless the caller already has it -- the onset
    // window is bracketed on the peak, so it can't be found without one.
    const int pUser = (pPeakIn >= 0.0) ? (int)std::lround(pPeakIn)
        : (int)std::lround(FeatureMarks::detect_p_peak(v, r_idx, fs));
    if (pUser < 0 || pUser >= N) return -1.0;
    const int w = win_005s(fs);
    const bool is_positive = FeatureMarks::qrs_positive_at(v, r_idx);
    std::vector<double> u = v;
    if (!is_positive) for (auto& x : u) x = -x;
    const int lo = std::max(0, pUser - w);
    const int hi = std::min(N - 1, pUser + w / 4);
    if (hi - lo < 4) return std::clamp(pUser, 0, N - 1);
    // B = pre-P baseline: left edge of window (the level the onset rises FROM).
    const double B = u[lo];
    // Spec I-3: P-onset is a transition onset -- 40-sample window, 4x cubic
    // upsample, fit-and-select via transitionAnchor; anchor at 0-20% of the
    // onset (fraction 0.10, baseline side). Mirrors compute_q_onset.
    return std::clamp(
        subsample_refine::transitionAnchor(u, pUser, 0.10, 40, B, lo, hi),
        0.0, static_cast<double>(N - 1));
}

// T-onset: the T wave's foot. Same shape as Q-onset and the J-point -- window
// bounded by an adjacent landmark, baseline from the window edge the wave
// departs from, anchor at 10% of the rise, then the 4x-upsampled
// fit-and-select. Window is [J-point, J-point + 100 ms].
double FeatureMarks::compute_t_begin(const std::vector<double>& v, double fs, int r_idx,
    double jPointIn) {
    const int N = static_cast<int>(v.size());
    if (N < 4 || r_idx < 0 || r_idx >= N) return -1.0;
    auto cl = [&](int i) { return std::clamp(i, 0, N - 1); };
    auto cld = [&](double d) { return std::clamp(d, 0.0, static_cast<double>(N - 1)); };

    // The J-point bounds the window: T-begin is always after it. Reuse the
    // caller's value when given -- recomputing costs a whole transitionAnchor.
    const double j_point = (jPointIn >= 0.0) ? jPointIn
        : compute_j_point(v, fs, r_idx);
    if (j_point < 0.0) return -1.0;

    const bool is_positive = FeatureMarks::qrs_positive_at(v, r_idx);
    std::vector<double> u = v;
    if (!is_positive) for (auto& x : u) x = -x;

    // Window: [J-point, J-point + 100 ms].
    const int lo = cl(static_cast<int>(std::ceil(j_point)));
    const int hi = cl(lo + static_cast<int>(std::lround(0.100 * fs)));
    if (hi - lo < 4) return cld(j_point);

    // B = recovered ST level at the left edge: the level the T departs FROM.
    const double baseline = u[lo];

    // 4x cubic-upsample transition fit-and-select; onset anchor at 10% of the
    // rise (baseline side), matching every other onset and offset here.
    return cld(subsample_refine::transitionAnchor(u, lo, 0.10, 40, baseline, lo, hi));
}

AnchorLocator make_anchor_locator(AnchorType type, int r_col, double fs) {
    switch (type) {
    case AnchorType::R_PEAK:  return [r_col](const std::vector<double>&) {
        return static_cast<double>(r_col); };
    case AnchorType::Q_ONSET:
        // Single canonical finder: detect_q_begin supplies the seed;
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
    case AnchorType::P_PEAK:  return [r_col, fs](const std::vector<double>& b) {
        return FeatureMarks::detect_p_peak(b, r_col, fs); };
    case AnchorType::P_ONSET:
        // compute_p_begin detects the P peak itself, so this is one call.
        return [r_col, fs](const std::vector<double>& b) {
            return FeatureMarks::compute_p_begin(b, fs, r_col);
            };
    case AnchorType::T_PEAK:  // no detect_t_peak: bracket via t_begin/t_end
        // One chain, each landmark bounding the next: J-point -> T-begin ->
        // T-end, then T-peak as the extremum between the two.
        return [r_col, fs](const std::vector<double>& b) {
            const double j = FeatureMarks::compute_j_point(b, fs, r_col);
            const double tbD = FeatureMarks::compute_t_begin(b, fs, r_col, j);
            const double teD = FeatureMarks::compute_t_end(b, fs, r_col, tbD);
            // The brackets stay fractional: compute_t_peak takes doubles.
            return FeatureMarks::compute_t_peak(b, tbD, teD);
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

FeatureMarks::ReactiveEcg FeatureMarks::reactive_ecg(
    const std::vector<double>& ecg, int t_begin, int t_end)
{
    ReactiveEcg r;
    // compute_t_peak folds in the sigma=15 asymmetric-extremum refinement, so
    // this is the fully-refined T-peak (the same one the bar path uses).
    if (static_cast<int>(ecg.size()) >= 3)
        r.t_peak = compute_t_peak(ecg, static_cast<double>(t_begin),
            static_cast<double>(t_end));
    return r;
}

FeatureMarks::ReactivePpg FeatureMarks::reactive_ppg(
    const std::vector<double>& ppg, int onset, int peak, int end)
{
    ReactivePpg r;
    if (static_cast<int>(ppg.size()) < 3) return r;
    if (onset >= 0 && peak > onset) r.t50 = amplitude_crossing(ppg, onset, peak, 0.50);
    if (peak >= 0 && end > peak) {
        r.t80 = amplitude_crossing(ppg, peak, end, 0.80);
        // T80_rise / PW80 at t80's OWN absolute level (see detect_ppg_
        // fiducials for the rationale): upslope crossing of the same value.
        if (onset >= 0 && peak > onset) {
            const double vp = sample_at(ppg, static_cast<double>(peak));
            const double ve = sample_at(ppg, static_cast<double>(end));
            if (std::isfinite(vp) && std::isfinite(ve)) {
                const double target = vp + 0.80 * (ve - vp);
                const double xr = crossing_at_level(ppg, onset, peak, target);
                if (xr >= 0.0) {
                    r.t80_rise = xr;
                    if (r.t80 >= 0.0 && r.t80 > r.t80_rise) r.pw80 = r.t80 - r.t80_rise;
                }
            }
        }
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

FeatureMarks::PpgFiducials FeatureMarks::detect_ppg_fiducials(const std::vector<double>& v, int W, double ppgRate,
    double heightMeters)
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
    // (a dead `w20 = 0.020 * ppgRate` sat here, declared and never read)
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

    // ---- DICROTIC NOTCH (three-tier, E-5) + DIASTOLIC PEAK -----------------
    // The notch comes from ppg_dicrotic::detectDicroticNotch (Tier 1 IEM ->
    // Tier 2 Windkessel -> Tier 3 absent), recording which tier answered. The
    // detector returns only the notch, so the diastolic peak (peak2) still comes
    // from the 4-knot cubic-spline fit over [amp15, amp90]. A not-found notch
    // leaves a midpoint placeholder flagged not-found (draws "o" not "x").
    {
        // amplitude_crossing and cubicSplineNotch work on columns, so the
        // BOUNDS are floored/ceiled from the fractional positions. The values
        // they produce are kept as-is; only the search grid is integral.
        const int peakCol = iFloor(g.peak), endCol = iCeil(g.end);
        // Column bounds for the spline search: the crossing itself is fractional,
        // but cubicSplineNotch needs integer knots, so these two are floored
        // explicitly rather than by implicit conversion.
        int amp15 = iFloor(amplitude_crossing(v, peakCol, endCol, 0.15));
        if (amp15 < 0) amp15 = peakCol;
        int amp90 = iFloor(amplitude_crossing(v, peakCol, endCol, 0.90));
        if (amp90 < 0) amp90 = cl(iFloor(0.5 * (g.peak + g.end)));
        const int lo = cl(amp15);
        const int hi = std::clamp(amp90, lo, Wc - 1);

        // ---- DIASTOLIC PEAK FIRST, then the notch bounded by it ------------
        //
        // THE ORDER HERE IS REVERSED from what it was, and the dependency it
        // used to imply was never real. The comment said peak2 was "first local
        // max AFTER the notch", but cubicSplineNotch(v, lo, hi, ...) fits over
        // [amp15, amp90] and the notch is not an input to it -- the notch
        // appeared only in the ACCEPTANCE TEST afterwards. So peak2 was always
        // computable first, and computing it first is what lets it bound the
        // notch search.
        //
        // Why that is better than correcting afterwards. Bounding the search
        // means a notch past the diastolic peak is never a candidate: the IEM
        // residual's local minima beyond the peak are not examined at all.
        // Correcting afterwards meant the detector could return such a point,
        // report a tier and a confidence for it, and then have peak2 moved to
        // accommodate it -- so a bad notch displaced a good peak, and the tier
        // and confidence recorded in the archive described a landmark that had
        // been overruled.
        int splineDiastolic = -1;
        subsample_refine::cubicSplineNotch(v, lo, hi, &splineDiastolic);
        const bool splineOk = (splineDiastolic > lo && splineDiastolic < hi);
        g.peak2_found = splineOk;
        g.peak2 = splineOk ? static_cast<double>(splineDiastolic)
            : cld(0.5 * (g.peak + hi));

        // Notch: three-tier detector. The template spans ~one cardiac cycle, so
        // RR is its visible duration. Enhancement is left off (gain 0) pending
        // real-record validation -- the 0.15 default fills clear notches.
        //
        // dnWindowHiSample carries the diastolic peak in as the search ceiling.
        // Only when the spline actually found one: the midpoint fallback above
        // is a placeholder, not a measurement, and bounding a real detector by a
        // placeholder would let a fabricated position suppress a genuine notch.
        // With no peak2, the detector keeps its 70%-of-RR ceiling and the
        // post-hoc invariant below is what holds the ordering.
        const double rrSeconds = double(Wc) / ppgRate;
        ppg_dicrotic::PpgConfig dnCfg;
        dnCfg.dnEnhanceGain = 0.0;
        if (splineOk) dnCfg.dnWindowHiSample = splineDiastolic;
        const ppg_dicrotic::DnResult dn =
            ppg_dicrotic::detectDicroticNotch(v, ppgRate, peakCol, rrSeconds, dnCfg);

        g.notch_found = (dn.tier != ppg_dicrotic::DnResult::ABSENT) && dn.index > 0;
        // DnResult carries a sub-sample refinement of its own (subSample, from
        // Phase 1 I-3). It was being discarded in favour of the integer index;
        // use it when finite.
        g.dicrotic = g.notch_found
            ? (std::isfinite(dn.subSample) ? dn.subSample
                : static_cast<double>(dn.index))
            : cld(0.5 * (lo + hi));
        g.dn_tier = static_cast<int>(dn.tier);
        g.dn_confidence = dn.confidence;

        // Belt and braces, and not redundant. The bound above constrains the
        // TIER 1/2 search, but the not-found path still assigns a midpoint
        // placeholder 0.5*(lo + hi) that is computed from the window and not
        // from peak2, and can therefore land past it. This catches that case,
        // and any future path that sets g.dicrotic without going through the
        // bounded search. It should now be a no-op on real data -- if the
        // notch_found counter shows it firing, the bound is not being applied.
        order_notch_before_peak2(v, g.dicrotic, g.peak2, g.peak2_found,
            static_cast<double>(hi));
    }

    // ---- T80 / P50: amplitude crossings (the same helper the GUI's reactive
    // T80/P50 glyphs call, so the two can't disagree). ----------------------
    g.t80 = amplitude_crossing(v, iFloor(g.peak), iCeil(g.end), 0.80);
    if (g.t80 < 0) g.t80 = cld(0.5 * (g.peak + g.end));

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
    g.p50 = amplitude_crossing(v, iFloor(g.onset), iCeil(g.peak), 0.50);
    if (g.p50 < 0) g.p50 = cld(0.5 * (g.onset + g.peak));

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

int FeatureMarks::detect_q_begin(const std::vector<double>& ecg_signal, int r_idx) {
    const int N = static_cast<int>(ecg_signal.size());
    const bool is_positive = FeatureMarks::qrs_positive_at(ecg_signal, r_idx);

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

// P peak, sub-sample refined. Coarse seed = argmax over the PR-region window
// (bounded below); refined with the asymmetric-extremum method (cubic fit on
// Gaussian-weighted samples, analytic derivative, sigma = 12, the P-peak
// sigma). Returns a FLOAT so callers that report/seed off the P peak inherit
// the sub-sample position (int-index callers round).
double FeatureMarks::detect_p_peak(const std::vector<double>& ecg_signal, int r_idx, double fs) {
    const int N = static_cast<int>(ecg_signal.size());
    if (N < 3) return 0.0;

    const bool is_positive = FeatureMarks::qrs_positive_at(ecg_signal, r_idx);
    std::vector<double> upright = ecg_signal;
    if (!is_positive) for (auto& x : upright) x = -x;

    // P-peak search window = [Q-begin - 200 ms, Q-begin]. Q-begin (detect_q_begin
    // is r_idx-relative, ~50 ms before R) is the right bound; the P wave sits in
    // the 200 ms preceding it (a full PR interval). Bounding the left edge this
    // way also stops the argmax from reaching back across an earlier beat's much-
    // larger QRS spike when called with a later r_idx on a multi-beat array.
    const int q_est = detect_q_begin(ecg_signal, r_idx);
    const int win = std::max(1, static_cast<int>(std::lround(0.200 * fs)));
    const int lo = std::max(0, q_est - win);
    const int hi = std::max(lo + 1, std::min(q_est, N));
    int best = -1; double bv = -std::numeric_limits<double>::infinity();
    for (int i = lo; i < hi; ++i)
        if (!std::isnan(upright[i]) && upright[i] > bv) { bv = upright[i]; best = i; }
    const int seed = (best >= 0) ? best : lo;
    // Refine: asymmetric extremum on the UPRIGHT signal (cubic on Gaussian-
    // weighted samples, analytic derivative), sigma = 12.
    return subsample_refine::asymmetricExtremum(upright, seed, 12.0);
}

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

    const int p_idx = (pPeakIn >= 0.0) ? (int)std::lround(pPeakIn)
        : (int)std::lround(detect_p_peak(ecg_signal, r_idx, fs));
    if (p_idx < 0 || p_idx >= N - 1)
        return std::clamp(p_idx + 1, 0, N - 1);

    // Post-P baseline: a short window just after the peak (P is much
    // shorter than T, so this window is smaller than detect_s_end's).
    const int pb_lo = std::min(p_idx + 20, N - 1);
    const int pb_hi = std::min(p_idx + 50, N);
    double baseline = upright[p_idx];
    if (pb_hi - pb_lo >= 5) {
        std::vector<double> w(upright.begin() + pb_lo, upright.begin() + pb_hi);
        std::nth_element(w.begin(), w.begin() + w.size() / 2, w.end());
        baseline = w[w.size() / 2];
    }

    const double p_val = upright[p_idx];
    const double depth = baseline - p_val;
    if (depth <= 0.0)
        return std::clamp(p_idx + 15, 0, N - 1);

    const double target = p_val + 0.90 * depth;
    const int hi = std::min(p_idx + 60, N);
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

// Diastolic peak (peak2): mirror of detect_ppg_dicrotic, searching for a
// local MAXIMUM (the diastolic bump) instead of a minimum, starting after
// the dicrotic notch (physiologically, diastolic peak follows the notch).
int FeatureMarks::detect_ppg_peak2(const std::vector<double>& pulse) {
    const int N = static_cast<int>(pulse.size());
    const int sysPeak = std::clamp((int)std::lround(detect_ppg_peak(pulse)), 0, std::max(0, N - 1));
    const int dic = detect_ppg_dicrotic(pulse, sysPeak);
    const int end = detect_ppg_end(pulse);
    if (dic < 0 || end < 0 || end - dic < 10)
        return std::clamp((dic + end) / 2, 0, N - 1);

    const int margin = std::max(2, (end - dic) / 10);
    const int lo = dic + margin;
    const int hi = end - 1;
    if (hi - lo < 3)
        return std::clamp(dic + (end - dic) / 3, 0, N - 1);

    int best = -1;
    double bestVal = -1e300;
    for (int i = lo + 1; i < hi; ++i) {
        if (pulse[i] >= pulse[i - 1] && pulse[i] >= pulse[i + 1]) {
            if (pulse[i] > bestVal) { bestVal = pulse[i]; best = i; }
        }
    }
    if (best < 0) return std::clamp(dic + (end - dic) / 3, 0, N - 1);
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
            const int p = FeatureMarks::detect_ppg_upstroke_peak(v, lo, hi + 1);
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

        peak = cl(static_cast<int>(std::lround(subsample_refine::symmetricExtremum(v, peak, 8.0))));

        // DIASTOLIC PEAK FIRST HERE TOO, so it can bound the notch search --
        // matching the PPG path exactly. These two were previously found by
        // independent searches over the same cycle in the other order, with no
        // ordering test between them, so an arterial pulse could carry a notch
        // after its diastolic peak and nothing anywhere would notice.
        if (peak2 < 0 && end > onset) {
            const int base = cl(onset);
            std::vector<double> cyc(v.begin() + base, v.begin() + cl(end) + 1);
            const int seed = FeatureMarks::detect_ppg_peak2(cyc);
            const double refined = subsample_refine::asymmetricExtremum(cyc, seed, 10.0);
            peak2 = cl(base + static_cast<int>(std::lround(refined)));
        }
        if (dicrotic < 0 && cl(end) > cl(onset) + 2) {
            const int base = cl(onset);
            // Truncate the cycle AT the diastolic peak. detect_ppg_dicrotic
            // takes no window argument, so the bound is expressed by shortening
            // its input -- the same effect as dnWindowHiSample on the PPG side,
            // reached the only way this detector allows.
            int top = cl(end);
            if (peak2 > base + 2 && peak2 < top) top = peak2;
            std::vector<double> cyc(v.begin() + base, v.begin() + top + 1);
            const int peakInCyc = std::clamp(peak - base, 0, (int)cyc.size() - 1);
            const int seed = FeatureMarks::detect_ppg_dicrotic(cyc, peakInCyc);
            const double refined = subsample_refine::asymmetricExtremum(cyc, seed, 10.0);
            dicrotic = cl(base + static_cast<int>(std::lround(refined)));
        }
        // Same backstop as the PPG path, same function, so the channels cannot
        // drift apart on the invariant. Covers the case where dicrotic or peak2
        // arrived pre-set from a caller rather than from the searches above.
        if (dicrotic >= 0 && peak2 >= 0)
            FeatureMarks::order_notch_before_peak2(v, dicrotic, peak2, cl(end));
        //the ppg foot is transition anchor
        onset = cl(static_cast<int>(std::lround(subsample_refine::transitionAnchor(v, onset, 0.0, 40, std::numeric_limits<double>::quiet_NaN(),
            onset, std::min(onset + 40, n - 1)))));
    }

    inline int clampToVisible(int idx, int visN) {
        return std::clamp(idx, 0, visN - 1);
    }

} // anonymous

void FeatureMarks::seed_all(TemplateBin& b, double sampleRate, double ppgRate, AnchorType anchor,
    double heightMeters) {
    // Per-anchor ECG user markers are seeded into this anchor's set.
    TemplateBin::MarkerSet& mk = b.marks(anchor);

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
            // pf carries sub-sample positions throughout. TemplateBin's
            // ppg_*_auto fields are int, and they are reached by
            // pointer-to-member tables typed `int TemplateBin::*` (the CSV
            // emission table below) as well as by the GUI marker path, so
            // widening them is a change to those tables and to the marking
            // file's layout -- not a local edit. Until that happens the
            // narrowing is done HERE, once, explicitly and visibly, rather
            // than by silent implicit conversion at fifteen assignments.
            //
            // Everything computed from the fiducials before this point (the
            // derived indices ba..si, the notch tier and confidence, and the
            // T80/P50 intervals) uses the fractional positions, so the
            // rounding costs marker DISPLAY precision, not feature precision.
            auto rnd = [](double x) {
                return (x < 0.0) ? -1 : static_cast<int>(std::lround(x));
                };
            b.ppg_peak_auto = rnd(pf.peak);
            b.ppg_onset_auto = rnd(pf.onset);
            b.ppg_peak2_auto = rnd(pf.peak2);           b.ppg_peak2_found_auto = pf.peak2_found;
            b.ppg_end_auto = rnd(pf.end);
            b.ppg_dicrotic_auto = rnd(pf.dicrotic);     b.ppg_dicrotic_found_auto = pf.notch_found;
            b.ppg_t80_auto = rnd(pf.t80);
            b.ppg_t50_auto = rnd(pf.p50);
            b.ppg_u_auto = rnd(pf.u);
            b.ppg_v_auto = rnd(pf.v);
            b.ppg_w_auto = rnd(pf.w);
            b.ppg_a_auto = rnd(pf.a);
            b.ppg_b_auto = rnd(pf.b);
            b.ppg_c_auto = rnd(pf.c);
            b.ppg_d_auto = rnd(pf.d);
            b.ppg_e_auto = rnd(pf.e);
            b.ppg_f_auto = rnd(pf.f);
            b.ppg_p1_auto = rnd(pf.p1);
            b.ppg_p2_auto = rnd(pf.p2);
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
            mk.p_peak_ch[c] = mk.q_begin_ch[c] = -1;
            mk.s_end_ch[c] = mk.t_begin_ch[c] = mk.t_end_ch[c] = mk.p_begin_ch[c] = -1;
            b.r_peak_ch[c] = -1;
            b.p_peak_auto_ch[c] = b.q_begin_auto_ch[c] = b.r_peak_auto_ch[c] = -1;
            b.s_end_auto_ch[c] = b.t_begin_auto_ch[c] = b.t_end_auto_ch[c] = -1;
            continue;
        }

        // Was `const double visN` -- a typo: every use is an int parameter
        // (clampToVisible takes int), so it narrowed on every call.
        const int visN = std::max(static_cast<int>(ecg.size()), 2);
        auto cl = [&](int x) { return clampToVisible(x, visN); };
        auto cld = [&](double d) { return std::clamp(d, 0.0, static_cast<double>(visN - 1)); };

        //The R peak which defines the template is the baseline for every other 
        const double r_peak_location = cld(subsample_refine::symmetricExtremum(ecg, cl(chs[c]->r_col_raw), 5.0));
        // The finders below take an INTEGER search anchor (their r_idx walks
        // the sample grid); the positions they return are doubles and are kept
        // as such. Truncating here is explicit so it is not mistaken for a
        // rounded RESULT.
        const int r_anchor = static_cast<int>(r_peak_location);

        // No shared band table or shared B_iso here any more: every landmark
        // below comes from its own canonical finder, each of which establishes
        // its own window and baseline internally. (compute_q_peak computes the
        // PQ isoelectric median itself, for its 0.1 mV depth test.)

        // Q-onset: single canonical finder (same one the movable bar uses via
        // the Q_ONSET locator). detect_q_begin supplies the seed; compute_q_onset
        // does the Q-peak search + sigma=4 refine + I-3 transitionAnchor
        // refinement internally.
        const double q_auto = cld(FeatureMarks::compute_q_onset(
            ecg, sampleRate, r_anchor));

        // S-end: mirror of Q-onset. Cap the LEFT bound at the S peak (first
        // turning point walking right from R = the S extremum) so E is the S,
        // not the R tail. Then follow the spec crossing over [S peak, ST
        // baseline]. To land at the END of the S wave -- the J-point, where the
        // wave has recovered to baseline (the spec's "80-100% of an offset") --
        // the target sits near baseline. In the literal formula L = B + f*(E-B)
        // that is f = 0.10 (== 90% recovered from the extremum), NOT 0.90 which
        // would sit at the trough. Recovery is monotonic S -> baseline, so the
        // forward crossing lands at the J-point.
        // S-end / J-point: single canonical finder (same one the movable bar
        // uses via the J_POINT locator). detect_s_end supplies the seed;
        // compute_s_end does the S-peak-capped window + fit-and-select + I-3
        // transitionAnchor refinement internally.
        const double s_auto = cld(FeatureMarks::compute_j_point(
            ecg, sampleRate, r_anchor));
        // T-onset: single canonical finder, same onset algorithm as Q-onset and
        // the J-point. s_auto IS this channel's J-point, so it bounds the window
        // and no second compute_j_point call is needed.
        const double tp_auto = cld(FeatureMarks::compute_t_begin(
            ecg, sampleRate, r_anchor, s_auto));

        // T-end: mirror of S-end. Use the shared T peak as the LEFT bound so E
        // is the T, then the spec crossing over [T peak, post-T baseline] with a
        // near-baseline target so the forward crossing lands where the T returns
        // to baseline = the actual T offset (not up at the T peak). Literal
        // L = B + f*(E-B) with f=0.10 == 90% recovered = "80-100% of an offset".
        // T-end: single canonical finder (same one the movable bar uses via
        // the T_PEAK-pass T-end locator). Window bounded by T-begin;
        // compute_t_end does the fit-and-select + I-3 transitionAnchor
        // refinement internally. Identical to the bar path.
        // tp_auto IS this channel's T-begin, so it bounds the window directly.
        const double te_auto = cld(FeatureMarks::compute_t_end(
            ecg, sampleRate, r_anchor, tp_auto));
        // P-peak: single canonical finder (asymmetric-extremum refinement,
        // sigma=12, folded into detect_p_peak). Same source the P-onset seed
        // uses below, so the reported peak and the onset seed are identical.
        // ONE detect_p_peak call: the reported peak and the P-onset seed are the
        // same value, so calling it twice was pure duplication (it also runs
        // detect_q_begin internally, so each call was two searches).
        const double p_auto_refined = cld(FeatureMarks::detect_p_peak(ecg, r_anchor, sampleRate));
        // P-begin: single canonical finder (same one the movable bar uses via
        // the P_ONSET locator). The peak above supplies the seed; compute_p_begin
        // does the window + fit-and-select + I-3 transitionAnchor refinement
        // internally. Identical to the bar path.
        const double pb_auto = cld(FeatureMarks::compute_p_begin(
            ecg, sampleRate, r_anchor, p_auto_refined));

        // Auto fields always updated.
        b.p_peak_auto_ch[c] = p_auto_refined;
        b.q_begin_auto_ch[c] = q_auto;
        b.r_peak_auto_ch[c] = r_peak_location;
        b.s_end_auto_ch[c] = s_auto;
        b.t_begin_auto_ch[c] = tp_auto;
        b.t_end_auto_ch[c] = te_auto;
        b.p_begin_auto_ch[c] = pb_auto;

        // User fields (per-anchor): only seed when unset for THIS anchor.
        // R peak is auto-only (flat) so it's always overwritten with fresh auto.
        // MarkerSet stores integer sample indices, so the double auto values
        // are rounded here (display/storage stays int; the *_auto_ch fields
        // above keep the sub-sample doubles).
        if (mk.p_peak_ch[c] < 0)  mk.p_peak_ch[c] = (int)std::lround(p_auto_refined);
        if (mk.q_begin_ch[c] < 0) mk.q_begin_ch[c] = (int)std::lround(q_auto);
        b.r_peak_ch[c] = (int)std::lround(r_peak_location);
        if (mk.s_end_ch[c] < 0)   mk.s_end_ch[c] = (int)std::lround(s_auto);
        if (mk.t_begin_ch[c] < 0) mk.t_begin_ch[c] = (int)std::lround(tp_auto);
        if (mk.t_end_ch[c] < 0)   mk.t_end_ch[c] = (int)std::lround(te_auto);
        if (mk.p_begin_ch[c] < 0) mk.p_begin_ch[c] = (int)std::lround(pb_auto);
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

// ---------------------------------------------------------------------------
// Ordering invariant: notch before diastolic peak. See feature_marks.hpp.
// ---------------------------------------------------------------------------
bool FeatureMarks::order_notch_before_peak2(const std::vector<double>& pulse,
    double notch, double& peak2, bool& peak2_found, double hi)
{
    // Nothing to enforce if either is absent. An absent notch does not
    // invalidate a diastolic peak: the notch is damped to nothing in stiff
    // vessels (E-5 tier 3 records exactly that) while the reflected wave
    // remains, so the pair {absent, present} is a real observation.
    if (notch < 0.0 || peak2 < 0.0) return true;
    if (peak2 > notch) return true;

    // Violated. Re-search for the first local maximum strictly after the notch,
    // in whole samples, then hand the seed back to the same sub-sample
    // refinement the original finder used so the returned position is of the
    // same kind as the one it replaces.
    const int n = static_cast<int>(pulse.size());
    const int lo = static_cast<int>(std::floor(notch)) + 1;
    const int top = std::min(n - 2, static_cast<int>(std::floor(hi)));
    for (int i = std::max(1, lo); i <= top; ++i) {
        if (std::isnan(pulse[i - 1]) || std::isnan(pulse[i]) || std::isnan(pulse[i + 1]))
            continue;
        if (pulse[i] >= pulse[i - 1] && pulse[i] >= pulse[i + 1]) {
            const double refined = subsample_refine::asymmetricExtremum(pulse, i, 10.0);
            // The refinement can pull the position back across the notch. If it
            // does, keep the integer seed rather than reintroduce the violation
            // this function exists to remove.
            peak2 = (refined > notch) ? refined : static_cast<double>(i);
            peak2_found = true;
            return true;
        }
    }

    // No maximum after the notch. Absent, and said so.
    peak2 = -1.0;
    peak2_found = false;
    return false;
}

bool FeatureMarks::order_notch_before_peak2(const std::vector<double>& pulse,
    int notch, int& peak2, int hi)
{
    double p2 = (peak2 >= 0) ? static_cast<double>(peak2) : -1.0;
    bool found = (peak2 >= 0);
    const bool ok = order_notch_before_peak2(pulse, static_cast<double>(notch),
        p2, found, static_cast<double>(hi));
    peak2 = (p2 >= 0.0) ? static_cast<int>(std::lround(p2)) : -1;
    return ok;
}
