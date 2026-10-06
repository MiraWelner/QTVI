/*feature_marks.cpp -- implementations for FeatureMarks.
See feature_marks.hpp for the public interface*/

#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstdio>
#include <cstdint>
#include <limits>
#include <numeric>
#include <vector>
#include <functional>

#include "feature_marks.hpp"
#include "sample_extent.hpp"
#include "fiducial_marker_finding\template_marking_bin_io.hpp"
#include "fiducial_marker_finding\anchor_view.hpp"
#include "subsample_refine.hpp"
#include "ppg_derivative.hpp"
#include "ppg_dicrotic.hpp"

namespace {
    // Search spans, in seconds.
    const double Q_PEAK_WIN_S = 0.04;
    const double Q_ONSET_WIN_S = 0.02;
    const double S_PEAK_WIN_S = 0.10;
    const double J_POINT_WIN_S = 0.10;
    const double T_END_WIN_LO_S = 0.100;
    const double T_END_WIN_HI_S = 0.500;
    //how deep does Q have to be for it to be found - else it is a circle and marked not found
    const double Q_MIN_DEPTH = 0.0075;

    //this carries over if the user clicked 'lead reversed' and the trace is inverted
    std::vector<double> upright_copy(const std::vector<double>& v, double sgn) {
        std::vector<double> u = v;
        if (sgn < 0.0) for (auto& x : u) x = -x;
        return u;
    }
}

double FeatureMarks::sample_at(const std::vector<double>& v, double p) {
    //get the amplitude at given subsampled position. Doesn't do the subsampling itself, just returns it
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

double FeatureMarks::find_q_peak(const std::vector<double>& ecg, int r_idx, double fs, double sgn, curve_fit::PeakFitMode peakMode) {
    /* The Q peak is found by: search the Q_PEAK_WIN_S window with the rightmost
       bound being the R peak, require the trough to be strictly interior, then
       gate it on depth relative to R before sub-sample refinement. */
    const int N = static_cast<int>(ecg.size());
    if (r_idx <= 0 || r_idx >= N) return -1;
    auto cl = [&](int i) { return std::clamp(i, 0, N - 1); };

    const std::vector<double> u = upright_copy(ecg, sgn);

    const int qp_lo = cl(r_idx - static_cast<int>(std::lround(Q_PEAK_WIN_S * fs)));
    const int qp_hi = cl(r_idx);
    int qSeed = qp_hi;
    double qv = std::numeric_limits<double>::infinity();
    for (int i = qp_lo; i <= qp_hi; ++i)
        if (!std::isnan(u[i]) && u[i] < qv) { qv = u[i]; qSeed = i; }
    // Strictly interior. Landing on qp_hi means the trace descends monotonically
    // into R with no notch -- no Q wave. Landing on qp_lo means the real minimum
    // is probably outside the window. An all-NaN window leaves qSeed at qp_hi,
    // which this same test rejects.
    if (qSeed <= qp_lo || qSeed >= qp_hi) return -1;

    double b_iso;
    {
        // PQ isoelectric level: median of [R - 100 ms, R - 40 ms]. Median, not
        // mean, so a P wave or a spike inside the span does not drag the level.
        const int a = cl(r_idx - static_cast<int>(std::lround(0.100 * fs)));
        const int b = cl(r_idx - static_cast<int>(std::lround(Q_PEAK_WIN_S * fs)));
        std::vector<double> s;
        for (int i = std::min(a, b); i <= std::max(a, b); ++i)
            if (!std::isnan(u[i])) s.push_back(u[i]);
        if (s.empty()) return -1;
        std::nth_element(s.begin(), s.begin() + s.size() / 2, s.end());
        b_iso = s[s.size() / 2];
    }
    const double rAmp = u[r_idx] - b_iso;
    if (rAmp > 0.0 && (b_iso - u[qSeed]) < Q_MIN_DEPTH * rAmp) return -1;   // too shallow to be a Q
    return std::clamp(upsample_for_fit::find_peak(u, qSeed,
        upsample_for_fit::peak_sigma::Q,
        upsample_for_fit::peak_halfwidth::Q, peakMode),
        0.0, static_cast<double>(N - 1));
}

// S = the trough after R on the upright copy. Search window is S_PEAK_WIN_S, so
// it needs no s_end bound -- the single S-trough source, used for the s_end
// detection and for |R|+|S| normalization.
double FeatureMarks::find_s_peak(const std::vector<double>& ecg, int r_idx, double fs, double sgn, curve_fit::PeakFitMode peakMode) {
    const int N = static_cast<int>(ecg.size());
    if (r_idx < 0 || r_idx >= N - 1)
        return static_cast<double>(std::clamp(r_idx + 1, 0, std::max(0, N - 1)));
    auto cl = [&](int i) { return std::clamp(i, 0, N - 1); };

    const std::vector<double> u = upright_copy(ecg, sgn);

    // Search range: [R, R + S_PEAK_WIN_S]. S is the trough just after R, so on
    // the upright copy it is the minimum.
    const int sp_lo = cl(r_idx);
    const int sp_hi = cl(r_idx + static_cast<int>(std::lround(S_PEAK_WIN_S * fs)));

    int sSeed = sp_lo;
    double sv = std::numeric_limits<double>::infinity();
    for (int i = sp_lo; i <= sp_hi; ++i)
        if (!std::isnan(u[i]) && u[i] < sv) { sv = u[i]; sSeed = i; }

    //return best fit (cubic or quadratic)
    return std::clamp(upsample_for_fit::find_peak(u, sSeed,
        upsample_for_fit::peak_sigma::S,
        upsample_for_fit::peak_halfwidth::S, peakMode),
        0.0, static_cast<double>(N - 1));
}


// NO sgn PARAMETER, AND THAT IS NOT AN OMISSION. A T wave inverted in an
// upright lead is ordinary pathology, not lead reversal, so this decides
// deflection direction locally -- from the sample at the bracket extremum
// versus the mean of the bracket ends. Lead polarity is the wrong input here.
// ---- THE REFINEMENT, WITH ITS CANDIDATES KEPT ---------------------------
//
// find_p_peak and find_t_peak both search a flipped copy when the deflection
// is negative, so the contest runs on `u` and its curves are upside down
// relative to the trace the caller draws. Negating the coefficients puts them
// back: eval() is a polynomial in the coefficients, so negating every one
// negates the curve, and `position` is a column and needs no conversion.
//
// Same shape as template_viewer_focus.cpp's `cf /= eref` for the normalized
// panel -- a unit conversion on the coefficients rather than on the curve.
static double refine_peak_keeping_candidates(const std::vector<double>& u,
    int seed, double sigma, int halfWidth, curve_fit::PeakFitMode peakMode,
    bool flipped, upsample_for_fit::PeakCandidates* cand)
{
    if (!cand)
        return upsample_for_fit::find_peak(u, seed, sigma, halfWidth, peakMode);

    *cand = upsample_for_fit::peakCandidates(u, seed, sigma, halfWidth, peakMode);
    if (flipped)
        for (auto& f : cand->draw)
            for (double& cf : f.coeff) cf = -cf;
    return cand->placement;
}

double FeatureMarks::find_t_peak(const std::vector<double>& v, double bracketSEnd, double bracketTEnd, curve_fit::PeakFitMode peakMode,
    upsample_for_fit::PeakCandidates* cand) {
    const int N = static_cast<int>(v.size());
    if (N < 1) return -1.0;

    int fFin = -1, lFin = -1;
    if (!sample_extent::finiteExtent(v, fFin, lFin))
        return -1.0;             // entirely NaN: no data
    const double loD = std::min(bracketSEnd, bracketTEnd);
    const double hiD = std::max(bracketSEnd, bracketTEnd);
    int lo = std::clamp(static_cast<int>(std::ceil(loD)), fFin, lFin);
    int hi = std::clamp(static_cast<int>(std::floor(hiD)), fFin, lFin);
    if (hi < lo) hi = lo;

    int a = lo, b = hi;
    const bool anyFinite = sample_extent::trimToFinite(v, a, b);
    (void)anyFinite;   // each caller's own sentinel follows
    if (b < a) return static_cast<double>(lo);   // window was a NaN gap

    // Baseline = mean of the (finite) window ends; both sit at the T's feet, so
    // their mean is a local isoelectric reference. The peak is the sample
    // furthest from it -- polarity-independent, so an inverted T works too.
    const double B = 0.5 * (v[a] + v[b]);
    int best = a; double bd = -1.0;
    for (int i = a; i <= b; ++i)
        if (!std::isnan(v[i]) && std::abs(v[i] - B) > bd) { bd = std::abs(v[i] - B); best = i; }

    std::vector<double> u = v;
    const bool flipped = (v[best] < B);
    if (flipped) for (auto& x : u) x = -x;
    return std::clamp(refine_peak_keeping_candidates(u, best,
        upsample_for_fit::peak_sigma::T, upsample_for_fit::peak_halfwidth::T,
        peakMode, flipped, cand), static_cast<double>(lo), static_cast<double>(hi));
}

// Local deflection direction, same reasoning as find_t_peak: an inverted P is a
// real waveform on an upright lead. No sgn parameter.
double FeatureMarks::find_p_peak(const std::vector<double>& v, double loIn, double hiIn, double fs, curve_fit::PeakFitMode peakMode,
    upsample_for_fit::PeakCandidates* cand)
{
    //found coarsely by the max deviation from the qonset,
    const int N = static_cast<int>(v.size());
    if (N < 1 || fs <= 0.0) return -1.0;

    int fFin = -1, lFin = -1;
    if (!sample_extent::finiteExtent(v, fFin, lFin))
        return -1.0;             // entirely NaN: no data

    const double loD = std::min(loIn, hiIn);
    const double hiD = std::max(loIn, hiIn) - 0.020 * fs;
    int lo = std::clamp(static_cast<int>(std::ceil(loD)), fFin, lFin);
    int hi = std::clamp(static_cast<int>(std::floor(hiD)), fFin, lFin);
    if (hi < lo) hi = lo;

    int a = lo, b = hi;
    const bool anyFinite = sample_extent::trimToFinite(v, a, b);
    (void)anyFinite;   // each caller's own sentinel follows
    if (b < a) return static_cast<double>(fFin);   // window was a NaN gap
    const double B = v[b];
    int best = a; double bd = -1.0;
    for (int i = a; i <= b; ++i) {
        if (std::isnan(v[i])) continue;
        const double d = std::abs(v[i] - B);
        if (d > bd) { bd = d; best = i; }
    }
    std::vector<double> u = v;
    const bool flipped = (v[best] < B);
    if (flipped) for (double& x : u) x = -x;
    return std::clamp(refine_peak_keeping_candidates(u, best,
        upsample_for_fit::peak_sigma::P, upsample_for_fit::peak_halfwidth::P,
        peakMode, flipped, cand), static_cast<double>(fFin), static_cast<double>(lFin));
}

double FeatureMarks::find_j_point(const std::vector<double>& v, double fs, int r_col, double sgn, upsample_for_fit::TransitionCandidates* candOut, curve_fit::FitMode mode) {
    const int N = static_cast<int>(v.size());
    if (r_col < 0 || r_col >= N - 1 || N < 4) return -1.0;
    auto cl = [&](int i) { return std::clamp(i, 0, N - 1); };
    auto cld = [&](double d) { return std::clamp(d, 0.0, static_cast<double>(N - 1)); };
    auto ms = [&](double s) { return static_cast<int>(std::lround(s * fs)); };
    const std::vector<double> u = upright_copy(v, sgn);
    const double sPeakD = FeatureMarks::find_s_peak(v, r_col, fs, sgn);
    const int sPeak = cl(static_cast<int>(std::floor(sPeakD)));
    const int lo = cl(sPeak);

    int hi = cl(sPeak + ms(J_POINT_WIN_S));
    // NOT FOUND, not the S peak: a window too short to fit used to return
    // the S peak as J, a substitute that looked like a measurement.
    if (hi - lo < 4) return -1.0;
    {
        const int sm = std::max(1, ms(0.006));
        auto slopeAt = [&](int i) {
            const int a = std::max(lo, i - sm), b = std::min(hi, i + sm);
            if (b <= a || std::isnan(u[a]) || std::isnan(u[b])) return 0.0;
            return (u[b] - u[a]) / static_cast<double>(b - a);
            };
        double s0 = 0.0; int k = lo + 1;
        for (; k <= hi; ++k) { s0 = slopeAt(k); if (std::abs(s0) > 1e-12) break; }
        for (int j = k + 1; j <= hi && k <= hi; ++j) {
            if (s0 * slopeAt(j) < 0.0) { if (j > lo + 3) hi = j; break; }
        }
    }
    if (hi - lo < 4) return -1.0;   // not found -- see above
    const double baseline = u[hi];
    return cld(upsample_for_fit::upsample_transition_and_curve_fit(u, sPeak, 0.10, 40, baseline, lo, hi, candOut, mode));
}

double FeatureMarks::find_q_onset(const std::vector<double>& v, double fs, int r_idx, double sgn, double qPeakIn, bool* measured,
    upsample_for_fit::TransitionCandidates* candOut, curve_fit::FitMode mode) {
    if (measured) *measured = false;   // set true only on the fit path below

    const int N = static_cast<int>(v.size());
    if (r_idx <= 0 || r_idx >= N || N < 4) return -1.0;
    auto cl = [&](int i) { return std::clamp(i, 0, N - 1); };
    auto cld = [&](double d) { return std::clamp(d, 0.0, static_cast<double>(N - 1)); };
    auto ms = [&](double s) { return static_cast<int>(std::lround(s * fs)); };

    const std::vector<double> u = upright_copy(v, sgn);

    // Steps 1-2: the Q peak, from the single canonical finder unless the caller
    // supplied one. -1 means no strict interior trough deep enough to be a Q
    // wave (monophasic R) -> R-upstroke fallback below.
    const double qPeakD = (qPeakIn >= 0.0) ? qPeakIn
        : FeatureMarks::find_q_peak(v, r_idx, fs, sgn);

    if (qPeakD >= 0.0) {
        // transitionAnchor needs an integer seed and integer bounds, so the
        // SEARCH runs on columns. floor, not lround: hi is this same column, so
        // rounding up would put the seed outside its own window when the
        // sub-sample peak sits just below an integer.
        const int qPeak = cl(static_cast<int>(std::floor(qPeakD)));

        // Step 4: Q-onset search range = [Q-peak - Q_ONSET_WIN_S, Q-peak].
        const int lo = cl(qPeak - ms(Q_ONSET_WIN_S));
        const int hi = cl(qPeak);
        if (hi - lo >= 4) {
            // Baseline = left edge (PQ baseline, the level the onset rises FROM).
            const double baseline = u[lo];
            // Step 5: 40-sample 4x cubic-upsample transition fit-and-select.
            // Onset anchor at 0-20% (fraction 0.10, baseline side).
            if (measured) *measured = true;
            return cld(upsample_for_fit::upsample_transition_and_curve_fit(u, qPeak, 0.10, 40, baseline, lo, hi, candOut, mode));
        }
        // Window too short to fit: NOT FOUND. It used to return the Q peak as
        // the onset -- a position that was not an onset measurement.
        return -1.0;
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

// ---- T ONSET ---------------------------------------------------------------
//
// THE ST SEGMENT'S RIGHT EDGE, AND THE T BAND'S LEFT. There was no finder for
// it, so every consumer that needed it left it at -1: the ST band came out
// empty and the T band either vanished or silently started at J, swallowing
// the ST segment.
//
// Same construction as the other onsets (Q onset, P begin): a transition fit
// on the rising side of the wave, anchored at 10% of the way from the level it
// rises FROM to its apex.
//   * B, the level it rises from, is the ST level: the median of the first
//     20 ms after J, so a sloped or shifted ST is the reference rather than
//     the PQ baseline the T never returns to on an elevated ST.
//   * E is the T apex: the sample furthest from B between J and T end, by
//     distance, so an inverted T is found the same way.
//   * The fit runs on [J, apex], so the only crossing it can find is the
//     upslope's.
// Bracketed by [J, apex] on the way out, so T onset can never precede the
// QRS end or follow the T peak. A T apex within 4 samples of J (no ST
// segment at all) returns J.
double FeatureMarks::find_t_begin(const std::vector<double>& v, double fs, double j_point, double t_end,
    upsample_for_fit::TransitionCandidates* candOut, curve_fit::FitMode mode) {
    const int N = static_cast<int>(v.size());
    if (N < 8 || !(fs > 0.0) || !(j_point >= 0.0) || !(t_end > j_point)) return -1.0;
    const int lo = std::clamp(static_cast<int>(std::ceil(j_point)), 0, N - 1);
    const int hiT = std::clamp(static_cast<int>(std::floor(t_end)), 0, N - 1);
    if (hiT - lo < 6) return -1.0;

    std::vector<double> st;
    const int stHi = std::min(hiT, lo + std::max(2, static_cast<int>(std::lround(0.020 * fs))));
    for (int i = lo; i <= stHi; ++i) if (std::isfinite(v[i])) st.push_back(v[i]);
    if (st.empty()) return -1.0;
    std::nth_element(st.begin(), st.begin() + st.size() / 2, st.end());
    const double B = st[st.size() / 2];

    int ePos = lo; double bestDist = 0.0;
    for (int i = lo; i <= hiT; ++i) {
        if (!std::isfinite(v[i])) continue;
        const double d = std::abs(v[i] - B);
        if (d > bestDist) { bestDist = d; ePos = i; }
    }
    if (ePos - lo < 4 || !(bestDist > 0.0)) return static_cast<double>(lo);

    const double tb = upsample_for_fit::upsample_transition_and_curve_fit(
        v, (lo + ePos) / 2, 0.10, 40, B, lo, ePos, candOut, mode);
    if (!std::isfinite(tb)) return -1.0;
    return std::clamp(tb, static_cast<double>(lo), static_cast<double>(ePos));
}

double FeatureMarks::find_t_end(const std::vector<double>& v, double fs, int r_col, double j_point,
    upsample_for_fit::TransitionCandidates* candOut, curve_fit::FitMode mode) {
    const int N = static_cast<int>(v.size());
    int lo0 = 0, hi = 0;
    FeatureMarks::t_end_window(v, fs, r_col, j_point, lo0, hi);
    // B = post-T baseline at the ceiling.
    const double B = v[hi];
    int ePos = lo0;
    double bestDist = 0.0;
    for (int i = lo0; i <= hi; ++i) {
        if (std::isnan(v[i])) continue;
        const double d = std::abs(v[i] - B);
        if (d > bestDist) { bestDist = d; ePos = i; }
    }
    const double E = v[ePos];
    const int lo = std::min(ePos, hi - 3);
    int fitHi = hi;
    {
        const double band = 0.05 * std::abs(E - B);
        const int hold = std::max(2, static_cast<int>(std::lround(0.040 * fs)));
        const int margin = std::max(2, static_cast<int>(std::lround(0.080 * fs)));
        for (int i = ePos + 1; i + hold <= hi; ++i) {
            bool settled = true;
            for (int k = i; k <= i + hold; ++k) {
                if (std::isfinite(v[k]) && std::abs(v[k] - B) > band) {
                    settled = false;
                    i = k;   // resume past the sample that broke the hold
                    break;
                }
            }
            if (settled) { fitHi = std::min(hi, i + margin); break; }
        }
        if (fitHi <= lo + 3) fitHi = hi;
    }
    const auto fit = curve_fit::selectBestFit(v, lo, fitHi);
    const double af = curve_fit::anchorAtFraction(fit, lo, fitHi, B, E, 0.05);
    const int seed = std::clamp(static_cast<int>(std::round(af)), 0, N - 1);
    const double te = upsample_for_fit::upsample_transition_and_curve_fit(v, seed, 0.02, 40, B, lo, fitHi, candOut, mode);
    return std::clamp(te, 0.0, static_cast<double>(N - 1));
}

double FeatureMarks::find_p_begin(const std::vector<double>& v, double fs, int r_idx, double sgn, double pPeakIn, upsample_for_fit::TransitionCandidates* candOut, curve_fit::FitMode mode) {
    const int N = static_cast<int>(v.size());
    const int fFin = sample_extent::firstFinite(v);
    if (fFin < 0) return -1.0;   // entirely NaN: no data
    double pPeak = pPeakIn;
    if (!(pPeak >= 0.0)) {
        const double q = FeatureMarks::find_q_onset(v, fs, r_idx, sgn);
        const double hi = (q >= 0.0) ? q : static_cast<double>(r_idx);
        pPeak = FeatureMarks::find_p_peak(v, std::max(static_cast<double>(fFin), static_cast<double>(r_idx) - 0.300 * fs), hi, fs);
    }
    const int pUser = std::clamp(static_cast<int>(std::lround(pPeak)), fFin, N - 1);
    const int w = std::max(1, static_cast<int>(std::lround(0.150 * fs))); //150ms window before the p peak - the p wave can be wide

    const std::vector<double> u = upright_copy(v, sgn);

    const int lo = std::max(fFin, pUser - w);
    const int hi = std::min(pUser, N - 1);
    const double B = u[std::clamp(lo, 0, N - 1)];
    const double pb = upsample_for_fit::upsample_transition_and_curve_fit(u, pUser, 0.10, 40, B, lo, hi, candOut, mode);

    if (!std::isfinite(pb)) return -1.0;
    return pb;
}

// THE LOCATORS CANNOT KNOW LEAD POLARITY. make_anchor_locator is applied per
// BEAT, from the alignment pass, which has no channel index to ask -- so the
// locators pass +1. An inverted lead's per-beat anchor search is therefore
// still done on the recorded polarity; the TEMPLATE detectors, which do have a
// channel, get the real sign. Threading it here means a sign on AnchorLocator's
// construction, which the alignment caller would have to supply per channel.
AnchorLocator make_anchor_locator(AnchorType type, int r_col, double fs) {
    switch (type) {
    case AnchorType::R_PEAK:  return [r_col](const std::vector<double>&) {
        return static_cast<double>(r_col); };
    case AnchorType::Q_ONSET:
        // Single canonical finder: find_q_onset does the Q-peak search, the
        // refine and the transitionAnchor refinement, and folds in the
        // R-upstroke fallback for monophasic-R (no-Q) beats.
        return [r_col, fs](const std::vector<double>& b) {
            return FeatureMarks::find_q_onset(b, fs, r_col, 1.0);
            };
    case AnchorType::J_POINT:
        // One call: find_j_point re-derives the S peak itself and places the
        // anchor via transitionAnchor. Sub-sample, not rounded.
        return [r_col, fs](const std::vector<double>& b) {
            return FeatureMarks::find_j_point(b, fs, r_col, 1.0);
            };
    case AnchorType::P_ONSET:
        // find_p_begin detects the P peak itself, so this is one call.
        return [r_col, fs](const std::vector<double>& b) {
            return FeatureMarks::find_p_begin(b, fs, r_col, 1.0);
            };
    case AnchorType::T_END:
        // The landmark chain the template detector runs: J bounds the T
        // search, so T end is found from this beat's own J.
        return [r_col, fs](const std::vector<double>& b) {
            const double j = FeatureMarks::find_j_point(b, fs, r_col, 1.0);
            if (!(j >= 0.0)) return -1.0;
            return FeatureMarks::find_t_end(b, fs, r_col, j);
            };
    }
    return [](const std::vector<double>&) { return -1.0; };
}


FeatureMarks::ReactiveEcg FeatureMarks::update_t_and_p_location(const std::vector<double>& ecg, double p_begin, double q_onset, double s_end, double t_end, double sampleRate, curve_fit::PeakFitMode peakMode)
{
    //called when the tend or pbegin bars move such that the ppeak and tpeak need to be recalculated
    ReactiveEcg r;
    r.p_peak = find_p_peak(ecg, p_begin, q_onset, sampleRate, peakMode, &r.p_peak_cand);
    r.t_peak = find_t_peak(ecg, s_end, t_end, peakMode, &r.t_peak_cand);
    return r;
}
FeatureMarks::ReactivePpg FeatureMarks::update_ppg_markings(const std::vector<double>& ppg, double onset, double peak, double dicrotic, double end)
{
    // Called when the onset, peak, dicrotic or end bars move, so that t50, t80,
    // t80_rise, pw80 and peak2 are recalculated.
    ReactivePpg r;

    // The crossing search steps through whole samples, so the bars are rounded
    // to give it a start and end column. The results are still sub-sample: the
    // crossing is interpolated between the two samples either side of it.
    int integer_onset = std::lround(onset);
    int integer_peak = std::lround(peak);
    int integer_end = std::lround(end);

    // Signal heights at the bars themselves (sample_at interpolates).
    const double onset_height = sample_at(ppg, onset);
    const double peak_height = sample_at(ppg, peak);
    const double end_height = sample_at(ppg, end);

    // t50: halfway up the rising edge.
    const double t50_height = onset_height + 0.50 * (peak_height - onset_height);
    r.t50 = signal_location_at_height(ppg, integer_onset, integer_peak, t50_height);

    // t80 and t80_rise: the same height on both edges, so pw80 is a true width.
    const double t80_height = peak_height + 0.80 * (end_height - peak_height);
    r.t80 = signal_location_at_height(ppg, integer_peak, integer_end, t80_height);
    r.t80_rise = signal_location_at_height(ppg, integer_onset, integer_peak, t80_height);
    if (r.t80 >= 0.0 && r.t80_rise >= 0.0 && r.t80 > r.t80_rise)
        r.pw80 = r.t80 - r.t80_rise;

    r.peak2 = detect_ppg_peak2(ppg, integer_peak, r.t80, integer_end);
    return r;
}

double FeatureMarks::signal_location_at_height(const std::vector<double>& v, int a, int b, double target) {
	//for features like t50, t80, t80_rise, pw80, and peak2, find the sub-sample location of a target height between two bars. 
    //Returns -1 if the target is not crossed between the bars.
    const int N = static_cast<int>(v.size());
    if (a < 0 || b <= a || b >= N) return -1.0;
    if (std::isnan(v[a]) || std::isnan(v[b]) || std::isnan(target)) return -1.0;

    const bool rising = v[b] >= v[a];
    for (int i = a + 1; i <= b; ++i) {
        const double prev = v[i - 1], cur = v[i];
        if (std::isnan(prev) || std::isnan(cur)) continue;

        const bool crossed = rising ? (prev < target && cur >= target)
            : (prev > target && cur <= target);
        if (crossed)
            return (i - 1) + (target - prev) / (cur - prev);
    }
    return -1.0;
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

// ---- EVERY PULSE LANDMARK EXISTS ------------------------------------------
//
// Detection first, then every landmark it did not find (-1) placed at the
// centre of its search range and flagged in .placeholder (kPhPpg*). The peak
// first, since the other three ranges hang off it. Only the four that are
// drawn as bars or the peak glyph are placed; the derived points (t50, t80,
// peak2, ...) are recomputed from these by reactive_ppg wherever they are
// shown.
FeatureMarks::PpgFiducials FeatureMarks::detect_ppg_fiducials(const std::vector<double>& v, int W, double ppgRate, double heightMeters, curve_fit::PeakFitMode peakMode, double measuredNotchCol)
{
    const int N = static_cast<int>(v.size());
    if (N < 1) return PpgFiducials{};          // no samples at all: nothing to place on

    // ---- DETECTION ---------------------------------------------------------
    // Immediately-invoked so its early returns end detection only; the
    // placeholder pass below always runs.
    PpgFiducials g = [&]() -> PpgFiducials {
        PpgFiducials g;
        if (N < 2) return g;                   // std::clamp(W, 2, N) needs N >= 2
        const int Wc = std::clamp(W, 2, N);    // visible window; nothing is ever placed past Wc-1
        // Fractional clamp, and column floor/ceil for the helpers that still need
        // an integer search grid. Positions themselves stay fractional.
        auto cld = [&](double x) { return std::clamp(x, 0.0, static_cast<double>(Wc - 1)); };
        auto iFloor = [&](double x) {
            return std::clamp(static_cast<int>(std::floor(x)), 0, Wc - 1);
            };
        auto iCeil = [&](double x) {
            return std::clamp(static_cast<int>(std::ceil(x)), 0, Wc - 1);
            };

        // Skip any leading NaN run: the template's first samples sit before the
        // first R (construction pads by `pad` seconds), so a partial pulse there
        // can win the slope gate and drag every landmark onto the left edge.
        int lo0 = std::max(0, sample_extent::firstFinite(v));
        if (lo0 > Wc) lo0 = Wc;

        int pkSeed = FeatureMarks::detect_ppg_upstroke_peak(v, lo0, Wc);
        if (pkSeed < 0) {
            // No detectable upstroke - use argmax for ppg peak
            int best = -1; double bv = -std::numeric_limits<double>::infinity();
            for (int i = lo0; i < Wc; ++i)
                if (!std::isnan(v[i]) && v[i] > bv) { bv = v[i]; best = i; }
            pkSeed = best;
        }
        if (pkSeed < 0) return g;              // all-NaN window: nothing to mark

        const int    peakHw = upsample_for_fit::pulse_window::peakHalfwidth(ppgRate);
        const double peakSg = upsample_for_fit::pulse_window::peakSigma(ppgRate);
        g.peak_cand = upsample_for_fit::peakCandidates(v, pkSeed,
            peakSg, peakHw, peakMode);
        // SEED STANDS WHEN THE CONTEST HAS NO ANSWER, exactly as the ECG path
        // leaves the detector's column standing. cld() on a -1 placement would
        // clamp it to column 0, which is inside the template's leading NaN pad --
        // see the note on refine_trough below for what that does downstream.
        g.peak = (g.peak_cand.valid && g.peak_cand.placement >= 0.0)
            ? cld(g.peak_cand.placement)
            : static_cast<double>(pkSeed);
        // The contest's own winner and placement stand; the panel colours the
        // curve that actually placed the mark, and reports "(rejected)" only on a
        // model the operator forced.
        g.peak_cand.placement = g.peak;
        if (g.peak < 3) { return g; }
        const int    footHw = upsample_for_fit::pulse_window::footHalfwidth(ppgRate);
        const double footSg = upsample_for_fit::pulse_window::footSigma(ppgRate);
        auto refine_trough = [&](int seed, int searchLo,
            upsample_for_fit::PeakCandidates& out) -> double {
                out = upsample_for_fit::peakCandidates(v, seed,
                    footSg, footHw, peakMode);
                const bool fitRan = out.valid && (out.placement >= 0.0);
                double pos = cld(std::max(
                    fitRan ? cld(out.placement)
                    : static_cast<double>(seed),
                    static_cast<double>(searchLo)));
                if (std::isnan(sample_at(v, pos))) return -1.0;
                // The placement follows the searchLo clamp so the panel's green
                // curve reports the column the mark is actually at; the WINNER is
                // the contest's, untouched.
                if (fitRan) out.placement = pos;
                return pos;
            };

        const int coarse_seed_for_onset = trough_in(v, lo0, iFloor(g.peak) - 1);
        g.onset = refine_trough(coarse_seed_for_onset, lo0, g.onset_cand);
        {
            const int seed = trough_in(v, std::min(iCeil(g.peak) + 1, Wc - 1), Wc - 1);
            // NOT FOUND IS -1, NOT Wc - 1. The old fallback pinned the end to the
            // LAST COLUMN OF THE ARRAY whenever no trough resolved after the peak
            // -- and Wc is the full template length, padded far tail included, so
            // that column sits outside the drawn extent by construction. The bar
            // then landed past the right edge of the pulse, where the marker loop
            // drops it: invisible and unclickable, and no amount of re-seeding
            // could rescue it because the freshly detected value was itself out of
            // range.
            //
            // -1 is the sentinel every consumer of PpgFiducials already handles,
            // for exactly the reason the refine_trough note above gives about 0:
            // Wc - 1 is a position, and a wrong one. An unresolvable end now reads
            // as no mark rather than as a mark at the wall.
            g.end = (seed >= 0) ? refine_trough(seed, lo0, g.end_cand) : -1.0;
        }
        {
            const double vp = sample_at(v, g.peak);
            const double ve = sample_at(v, g.end);
            g.t80_rise_y = vp + 0.80 * (ve - vp);
            g.t80 = signal_location_at_height(v, iFloor(g.peak), iCeil(g.end), g.t80_rise_y);
        }
        if (g.t80 < 0) g.t80 = cld(0.5 * (g.peak + g.end));
        {
            const bool have = (measuredNotchCol >= 0.0)
                && (measuredNotchCol <= static_cast<double>(Wc - 1));
            g.notch_found = have;
            if (have) {
                g.dicrotic = measuredNotchCol;
                g.dn_tier = 1;
            }
            else {
                // HALFWAY BETWEEN PEAK AND t80, so the circle and the bar have
                // somewhere to sit; notch_found = false says it is not a
                // measurement. Half the interval rather than a fixed offset:
                // peak + 0.12 * rate is 120 ms at any heart rate and on a short
                // cycle lands past t80, outside the interval it belongs in.
                // -1 when t80 collapses onto the peak, leaving nothing to bisect.
                g.dicrotic = (g.t80 > g.peak) ? cld(0.5 * (g.peak + g.t80)) : -1.0;
                g.dn_tier = 3;
            }
            g.dn_confidence = std::numeric_limits<double>::quiet_NaN();
        }

        g.peak2 = detect_ppg_peak2(v, iFloor(g.dicrotic), g.t80, iFloor(g.end));
        {
            const double vo = sample_at(v, g.onset);
            g.t50 = signal_location_at_height(v, iFloor(g.onset), iCeil(g.peak), vo + 0.50 * (sample_at(v, g.peak) - vo));
        }
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
        }();

    // ---- PLACEHOLDERS ------------------------------------------------------
    const int Wc = std::clamp(W, 1, N);
    const double first = std::clamp(static_cast<double>(std::max(0, sample_extent::firstFinite(v))),
        0.0, static_cast<double>(Wc - 1));
    const double last = static_cast<double>(Wc - 1);
    auto centre = [&](double lo, double hi) {
        lo = std::clamp(lo, first, last);
        hi = std::clamp(hi, first, last);
        if (hi < lo) std::swap(lo, hi);
        return std::round(0.5 * (lo + hi));
        };
    if (!(g.peak >= 0.0)) { g.peak = centre(first, last);         g.placeholder |= kPhPpgPeak; }
    if (!(g.onset >= 0.0)) { g.onset = centre(first, g.peak);      g.placeholder |= kPhPpgOnset; }
    if (!(g.end >= 0.0)) { g.end = centre(g.peak, last);         g.placeholder |= kPhPpgEnd; }
    if (!(g.dicrotic >= 0.0)) { g.dicrotic = centre(g.peak, g.end);   g.placeholder |= kPhPpgDicrotic; }
    return g;
}



// =========================================================================
// Movable (auto-detected seeds)
// =========================================================================

// P-end: walk forward from the P peak until the signal recovers to within
// 10% of a post-P baseline estimate. Needed for the PQ segment (spec: "end
// of P to immediately before Q onset"), which is a genuinely different
// landmark from P-onset -- no prior detector existed for it.
int FeatureMarks::find_p_end(const std::vector<double>& ecg_signal, int r_idx, double fs, double sgn,
    double pPeakIn) {
    const int N = static_cast<int>(ecg_signal.size());
    // GUARD FIRST. The polarity lookup used to sit above this, so it ran with
    // an out-of-range r_idx on the early-return path.
    if (r_idx < 0 || r_idx >= N)
        return std::clamp(r_idx, 0, std::max(0, N - 1));

    const std::vector<double> upright = upright_copy(ecg_signal, sgn);

    const double p_idx_d = (pPeakIn >= 0.0) ? pPeakIn
        : find_p_peak(ecg_signal, 0.0, find_q_onset(ecg_signal, fs, r_idx, sgn), fs);
    if (p_idx_d < 0.0 || p_idx_d >= static_cast<double>(N - 1))
        return static_cast<int>(p_idx_d) + 1;

    // One integer column, so the windows below are not silently narrowing a
    // double at every use. Truncation matches what the double indexing did.
    const int p_idx = static_cast<int>(p_idx_d);

    // Post-P baseline: a short window just after the peak (P is much shorter
    // than T, so this window is smaller than the T-offset equivalent).
    //
    // SAMPLES, NOT MILLISECONDS -- these four offsets (20, 50, 15, 60) do not
    // scale with fs, so at 1000 Hz the window is 20-50 ms and at 250 Hz it is
    // 80-200 ms. fs is already a parameter; converting them is a detector
    // change and belongs on its own.
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
int FeatureMarks::detect_ppg_upstroke_peak(const std::vector<double>& v, int lo, int hi)
{
    const int n = static_cast<int>(v.size());
    if (hi <= 0 || hi > n) hi = n;
    lo = std::max(0, lo);
    if (hi - lo < 5) return -1;

    const int h = std::max(3, (hi - lo) / 50);

    // Smoothed central-difference derivative; NaN-safe.
    //
    // SIZED TO THE WINDOW, NOT THE SIGNAL. This allocated and NaN-filled a
    // vector as long as the whole trace and then used [lo, hi) of it. Called
    // once per landmark that is nothing; called once per BEAT by
    // ppg_realign::upstrokePctCol, which is how a peak-anchored re-stack pays
    // for it -- a few hundred allocations of a few thousand doubles per
    // column, all but the window wasted. The foot-anchored path returns before
    // reaching it (pct == 0 needs no crossing), which is why only peak
    // alignment felt slow.
    //
    // `at(i)` DOES THE INDEX SHIFT IN ONE PLACE. Every loop below still runs
    // in signal coordinates, so the arithmetic cannot drift between them, and
    // the results are identical to the full-length version.
    const int dlo = lo, dn = hi - lo;
    std::vector<double> d(static_cast<size_t>(dn),
        std::numeric_limits<double>::quiet_NaN());
    const auto at = [&d, dlo](int i) -> double& {
        return d[static_cast<size_t>(i - dlo)];
        };

    for (int i = std::max(lo, h); i + h < hi; ++i) {
        if (std::isnan(v[i - h]) || std::isnan(v[i + h])) continue;
        at(i) = (v[i + h] - v[i - h]) / (2.0 * h);
    }

    double maxSlope = 0.0;
    for (int i = lo; i < hi; ++i)
        if (!std::isnan(at(i)) && at(i) > maxSlope) maxSlope = at(i);
    if (maxSlope <= 0.0) return -1;                 // flat / no rise
    const double gate = 0.25 * maxSlope;
    const int minRun = std::max(2, h);

    int anchor = -1;
    for (int i = lo; i < hi; ++i) {
        if (std::isnan(at(i)) || at(i) < gate) continue;
        int j = i, bestJ = i, held = 0;
        double bestD = at(i);
        while (j < hi && (std::isnan(at(j)) || at(j) >= gate)) {
            if (!std::isnan(at(j))) {
                ++held;
                if (at(j) > bestD) { bestD = at(j); bestJ = j; }
            }
            ++j;
        }
        if (held >= minRun) { anchor = bestJ; break; }
        i = j;                                      // too brief: keep looking
    }
    if (anchor < 0) return -1;

    int pk = -1;
    for (int i = anchor; i + 1 < hi; ++i) {
        if (std::isnan(at(i))) continue;
        if (at(i) <= 0.0) { pk = i; break; }
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
// upstroke (was argmax, i.e. the tallest peak); refined by find_peak -- THE
// SAME FUNCTION EVERY ECG PEAK USES, guarded quadratic vs cubic on BIC with
// the 5-point parabola as the fallback. Returns a FLOAT position so downstream
// fiducials that key off the peak (onset/t80/dicrotic/end brackets) inherit
// the sub-sample peak.
//
// This called quadratic_fit directly, with 8.0 written out where the pulse
// sigma belongs and no contest at all -- a third placement rule for one
// landmark, next to detect_ppg_fiducials' and the ECG's.
double FeatureMarks::detect_ppg_peak(const std::vector<double>& pulse,
    double ppgRate, curve_fit::PeakFitMode peakMode) {
    if (pulse.empty()) return 0.0;
    const int seed = detect_ppg_upstroke_peak(pulse);
    if (seed < 0) return 0.0;
    const int hw = upsample_for_fit::pulse_window::peakHalfwidth(ppgRate);
    const double pos = upsample_for_fit::find_peak(pulse, seed,
        upsample_for_fit::pulse_window::sigma(hw), hw, peakMode);
    // Seed stands on a failed refinement, as everywhere else.
    return (std::isfinite(pos) && pos >= 0.0)
        ? pos : static_cast<double>(seed);
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

void FeatureMarks::seed_all(time_bin& b, double sampleRate, double ppgRate,
    double abpRate, double artRate, double artPulmRate, AnchorType anchor,
    const LeadPolarity& pol, double heightMeters,
    curve_fit::FitMode fitMode, curve_fit::PeakFitMode peakMode) {
    // Per-anchor ECG user markers are seeded into this anchor's set.
    // (No bin-wide marker handle. Landmarks are per (lead, slot, anchor) now --
    // see TemplateBin::slotMarks -- so the set is fetched inside the per-channel
    // loop below, where the lead is known.)

    // ---- PPG ------------------------------------------------------------
    if (b.ppgTemplate.empty()) {
        b.bad_ppg = 2;
        // NO BIN-LEVEL BARS TO CLEAR. They are per slot now; a slot with no
        // pulse gets its -1s from applyTemplateToWidget's else branch and from
        // seed_pulse_bank_template declining to run.
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
        // OPERATOR-BAD PULSE: nothing to clear at bin level any more, and the
        // *_auto positions are deliberately left alone -- they are the
        // original detections and the verdict does not unmake them. The bars
        // this used to blank are per slot; badPulseMarked() on the slot is
        // what suppresses them.
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

        // ---- PPG autodetect positions: ONE call, single source of truth
        // for every fiducial. Nothing here recomputes anything independently
        // -- the GUI's frozen glyphs read straight from these same *_auto
        // fields (see BinPlotWidget::captureGlyphSnapshot), so "peak" (etc.)
        // can never mean two different things in two different places.
        {
            const auto pf = FeatureMarks::detect_ppg_fiducials(v, W, ppgRate, heightMeters, peakMode);
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
            // the 80%-downslope level and the width between them.
            // (t80_rise / pw80 were cached on the bin here. Both are pure
            //  functions of the three bars plus the peak, so the panel and the
            //  CSV derive them per slot through reactive_ppg instead.)
        }

        // ---- NO BAR SEEDING HERE ANY MORE ---------------------------------
        //
        // This seeded the bin-level bars from the bin-level autos. The bars are
        // per slot, and seedOneBin -> seed_pulse_bank_template seeds each slot
        // from that slot's OWN detection on that slot's OWN waveform -- which
        // is the reason the bars moved. Seeding a bin-wide position into a
        // slot's bar was the defect, not the mechanism.
        //
        // The *_auto fields above are still filled: they are the bin-wide
        // detection, still in the markings CSV's auto-feature columns and
        // still written to the .bin.
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

        // ONE DETECTOR, shared with the bank templates and the archive, and
        // THIS CHANNEL'S OWN POLARITY -- pol is indexed by c so the sign and
        // the lead cannot get out of step.
        const FeatureMarks::TemplateLandmarks lmRaw = FeatureMarks::detect_template_landmarks(
            ecg, chs[c]->r_col_raw, sampleRate, pol.sign(c), fitMode, peakMode);

        b.p_peak_auto_ch[c] = lmRaw.p_peak;
        b.q_peak_auto_ch[c] = lmRaw.q_peak;
        b.q_onset_auto_ch[c] = lmRaw.q_onset;
        b.q_onset_found_auto_ch[c] = lmRaw.q_onset_found;
        b.r_peak_auto_ch[c] = lmRaw.r_peak;
        b.s_end_auto_ch[c] = lmRaw.s_end;
        b.t_end_auto_ch[c] = lmRaw.t_end;
        b.p_begin_auto_ch[c] = lmRaw.p_begin;

        // R falls back to the unrefined column rather than -1: it is the
        // alignment anchor every other landmark is expressed against, so the
        // bin needs SOME R even when refinement could not run.
        b.r_peak_ch[c] = (lmRaw.r_peak >= 0.0)
            ? lmRaw.r_peak
            : std::clamp(static_cast<double>(chs[c]->r_col_raw),
                0.0, static_cast<double>(ecg.size()) - 1.0);
    }


    // ---- Arterial (ABP / ART / ART_PULM) --------------------------------
    // No polarity dimension: a pressure waveform has a physical sign.
    // THE RATE IS A PARAMETER, not the captured sampleRate. It used to be the
    // latter, which is the bug: sampleRate is the ECG rate and these are
    // arterial channels with their own.
    auto seedArterial = [&](const std::vector<double>& trace, double rate,
        uint8_t& issue,
        double& onset, double& peak, double& dicrotic, double& peak2, double& end,
        double& onset_auto, double& peak_auto, double& dic_auto, double& p2_auto, double& end_auto)
        {
            // AN ABSENT RATE IS AN ABSENT CHANNEL, reported the way an empty
            // trace is -- which gets issue = 2 right, as the bail alone would
            // not.
            if (trace.empty() || !(rate > 0.0)) {
                issue = 2;
                onset = peak = dicrotic = peak2 = end = -1;
                onset_auto = peak_auto = dic_auto = p2_auto = end_auto = -1;
                return;
            }
            if (issue == 1) {
                onset = peak = dicrotic = peak2 = end = -1;
                return;
            }
            const FeatureMarks::PpgFiducials pf = FeatureMarks::detect_ppg_fiducials(
                trace, static_cast<int>(trace.size()), rate, NAN, peakMode);
            onset_auto = pf.onset; peak_auto = pf.peak; dic_auto = pf.dicrotic;
            p2_auto = pf.peak2; end_auto = pf.end;
            if (onset < 0) onset = pf.onset;
            if (peak < 0) peak = pf.peak;
            if (dicrotic < 0) dicrotic = pf.dicrotic;
            if (peak2 < 0) peak2 = pf.peak2;
            if (end < 0) end = pf.end;
        };

    seedArterial(b.abpTemplate, abpRate, b.abp_issue,
        b.abp_onset, b.abp_peak, b.abp_dicrotic, b.abp_peak2, b.abp_end,
        b.abp_onset_auto, b.abp_peak_auto, b.abp_dicrotic_auto,
        b.abp_peak2_auto, b.abp_end_auto);
    seedArterial(b.artTemplate, artRate, b.art_issue,
        b.art_onset, b.art_peak, b.art_dicrotic, b.art_peak2, b.art_end,
        b.art_onset_auto, b.art_peak_auto, b.art_dicrotic_auto,
        b.art_peak2_auto, b.art_end_auto);
    seedArterial(b.artPulmTemplate, artPulmRate, b.art_pulm_issue,
        b.art_pulm_onset, b.art_pulm_peak, b.art_pulm_dicrotic,
        b.art_pulm_peak2, b.art_pulm_end,
        b.art_pulm_onset_auto, b.art_pulm_peak_auto, b.art_pulm_dicrotic_auto,
        b.art_pulm_peak2_auto, b.art_pulm_end_auto);
}

FeatureMarks::TemplateLandmarks FeatureMarks::detect_template_landmarks(
    const std::vector<double>& tmplIn, int nominal_r_col, double sampleRate, double sgn,
    curve_fit::FitMode fitMode, curve_fit::PeakFitMode peakMode)
{
    TemplateLandmarks out;
    const int n = static_cast<int>(tmplIn.size());
    if (n < 2 || nominal_r_col < 0 || nominal_r_col >= n || sampleRate <= 0.0)
        return out;
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
    double r = upsample_for_fit::find_peak(tmplIn, seed,
        upsample_for_fit::peak_sigma::R,
        upsample_for_fit::peak_halfwidth::R, peakMode);
    if (std::isnan(r) || r < 0.0 || r > static_cast<double>(n - 1))
        r = static_cast<double>(seed);   // refinement failed; nominal stands
    const int r_anchor = static_cast<int>(r);

    const double j = FeatureMarks::find_j_point(tmpl, sampleRate, r_anchor, sgn, &out.s_end_cand, fitMode);
    const double te = FeatureMarks::find_t_end(tmpl, sampleRate, r_anchor, j, &out.t_end_cand, fitMode);
    const double qp = FeatureMarks::find_q_peak(tmpl, r_anchor, sampleRate, sgn, peakMode);
    bool qFound = false;
    const double q = FeatureMarks::find_q_onset(tmpl, sampleRate, r_anchor, sgn, qp, &qFound, &out.q_onset_cand, fitMode);
    double pHi = q;
    if (!(pHi >= 0.0)) pHi = r - 0.050 * sampleRate;
    const double p_peak = FeatureMarks::find_p_peak(tmpl, std::max(0.0, static_cast<double>(r_anchor) - 0.300 * sampleRate), pHi, sampleRate, peakMode);
    const double pb = FeatureMarks::find_p_begin(tmpl, sampleRate, r_anchor, sgn, -1.0, &out.p_begin_cand, fitMode);

    // Out-of-range is folded to -1 (absent), NOT clamped to an edge column. A
    // landmark pinned to column 0 is indistinguishable from one genuinely found
    // there, and downstream would integrate over a window that does not exist.
    auto keep = [&](double x) {
        if (std::isnan(x) || x < 0.0 || x > static_cast<double>(n - 1))
            return -1.0;
        return x;
        };

    // R AGAIN, NOW THAT ITS BRACKETS EXIST. The pass above refined the
    // alignment's nominal column, and on a negative complex the refiner walks
    // uphill away from the true R. With q_onset and s_end in hand, R is the
    // sample furthest from the mean of those two ends -- the same rule
    // find_t_peak uses on s_end/t_end. LOCAL, so it is left alone by the sgn
    // threading: it is deciding which way THIS complex deflects between its own
    // brackets, not which way the lead was recorded.
    if (q >= 0.0 && j > q) {
        const int qa = std::clamp((int)std::lround(q), 0, n - 1);
        const int jb = std::clamp((int)std::lround(j), 0, n - 1);
        if (jb > qa && !std::isnan(tmpl[qa]) && !std::isnan(tmpl[jb])) {
            const double B = 0.5 * (tmpl[qa] + tmpl[jb]);
            int best = qa; double bd = -1.0;
            for (int i = qa; i <= jb; ++i) {
                if (std::isnan(tmpl[i])) continue;
                const double d = std::abs(tmpl[i] - B);
                if (d > bd) { bd = d; best = i; }
            }
            std::vector<double> u = tmpl;
            if (tmpl[best] < B) for (double& x : u) x = -x;
            const double rr = upsample_for_fit::find_peak(
                u, best, upsample_for_fit::peak_sigma::R,
                upsample_for_fit::peak_halfwidth::R, peakMode);
            if (std::isfinite(rr) && rr >= 0.0 && rr <= (double)(n - 1)) r = rr;
        }
    }

    out.r_peak = r;
    out.q_peak = keep(qp);   // -1 on a monophasic R is the RIGHT answer
    out.q_onset = keep(q);
    out.s_end = keep(j);
    out.t_end = keep(te);
    out.p_peak = keep(p_peak);   // -1 on a ventricular template is the RIGHT answer
    out.p_begin = keep(pb);
    // A flag that outlives its position would be a lie, so it is anded with the
    // position surviving keep().
    out.q_onset_found = qFound && (out.q_onset >= 0.0);
    out.valid = true;
    // NO -1 LEAVES HERE: a bar not found is placed at the centre of the range
    // it was searched in, and flagged (TemplateLandmarks::placeholder).
    return with_bar_placeholders(out, tmplIn, sampleRate);
}

// find_t_end's window: [J + 100 ms, J + 700 ms], trimmed to the last real
// sample. Fixed from J -- nothing to do with the next R. ONE DEFINITION, used
// by find_t_end and by the T-end placeholder, so a placeholder sits at the
// centre of exactly the range T end was searched in.
void FeatureMarks::t_end_window(const std::vector<double>& v, double fs, int r_col,
    double j_point, int& lo, int& hi)
{
    const int N = static_cast<int>(v.size());
    lo = hi = 0;
    if (N < 4 || r_col < 0 || r_col >= N) return;
    auto cl = [&](int i) { return std::clamp(i, 0, N - 1); };
    const int lastFin = sample_extent::lastFinite(v);
    if (lastFin < 0) return;
    lo = cl(static_cast<int>(std::lround(j_point + 0.100 * fs)));
    hi = std::min(cl(static_cast<int>(std::lround(j_point + 0.700 * fs))), lastFin);
}

FeatureMarks::TemplateLandmarks FeatureMarks::with_bar_placeholders(
    const TemplateLandmarks& in, const std::vector<double>& tmpl, double fs)
{
    TemplateLandmarks lm = in;
    lm.placeholder = 0;
    const int n = static_cast<int>(tmpl.size());
    if (!lm.valid || n < 2 || fs <= 0.0 || !(lm.r_peak >= 0.0)) return lm;
    const double R = lm.r_peak;
    // The usable columns: detect_template_landmarks blanks 40 ms at each end.
    const double edge = std::max(1.0, std::round(0.040 * fs));
    const double first = std::min(edge, static_cast<double>(n - 1));
    const double last = std::max(first, static_cast<double>(n - 1) - edge);
    auto centre = [&](double lo, double hi) {
        lo = std::clamp(lo, first, last);
        hi = std::clamp(hi, first, last);
        if (hi < lo) std::swap(lo, hi);
        return std::round(0.5 * (lo + hi));
        };
    if (!(lm.q_onset >= 0.0)) {
        lm.q_onset = (lm.q_peak >= 0.0)
            ? centre(lm.q_peak - Q_ONSET_WIN_S * fs, lm.q_peak)
            : centre(R - (Q_PEAK_WIN_S + Q_ONSET_WIN_S) * fs, R);
        lm.placeholder |= kPhQOnset;
    }
    if (!(lm.s_end >= 0.0)) {
        lm.s_end = centre(R, R + (S_PEAK_WIN_S + J_POINT_WIN_S) * fs);
        lm.placeholder |= kPhSEnd;
    }
    if (!(lm.t_end >= 0.0)) {
        int lo = 0, hi = 0;
        t_end_window(tmpl, fs, static_cast<int>(std::lround(R)), lm.s_end, lo, hi);
        lm.t_end = centre(lo, std::max(lo, hi));
        lm.placeholder |= kPhTEnd;
    }
    if (!(lm.p_begin >= 0.0)) {
        lm.p_begin = (lm.p_peak >= 0.0)
            ? centre(lm.p_peak - 0.150 * fs, lm.p_peak)
            : centre(R - 0.450 * fs, std::min(lm.q_onset, R - 0.050 * fs));
        lm.placeholder |= kPhPBegin;
    }
    return lm;
}

void FeatureMarks::seed_bank_template(const std::vector<double>& tmpl, int r_col,
    double sampleRate, double sgn, AnchorType anchor, tbank::BankMarkerSet& out,
    curve_fit::FitMode fitMode, curve_fit::PeakFitMode peakMode)
{
    out = tbank::BankMarkerSet{};          // all -1
    const TemplateLandmarks lm =
        FeatureMarks::detect_template_landmarks(tmpl, r_col, sampleRate, sgn, fitMode, peakMode);
    if (!lm.valid) return;
    if (anchor_view::showsBar(anchor, anchor_view::q_begin)) out.q_onset = lm.q_onset;
    if (anchor_view::showsBar(anchor, anchor_view::j_point)) out.s_end = lm.s_end;
    if (anchor_view::showsBar(anchor, anchor_view::t_end))   out.t_end = lm.t_end;
    // lm.p_begin is the -1 call now (see detect_template_landmarks), so the
    // override that used to live here is gone: one source again.
    if (anchor_view::showsBar(anchor, anchor_view::p_begin)) out.p_begin = lm.p_begin;
}


void FeatureMarks::seed_pulse_bank_template(const std::vector<double>& tmpl, double ppgRate, tbank::BankPulseMarkerSet& out, double heightMeters,
    curve_fit::PeakFitMode peakMode, double measuredNotchCol) {
    out = tbank::BankPulseMarkerSet{};
    const int W = static_cast<int>(tmpl.size());
    if (W < 3 || ppgRate <= 0.0) return;

    const PpgFiducials pf = detect_ppg_fiducials(tmpl, W, ppgRate, heightMeters, peakMode,
        measuredNotchCol);
    out.onset_auto = pf.onset;
    out.onset = pf.onset;
    out.peak_auto = pf.peak;
    out.dicrotic_auto = pf.dicrotic;
    out.dicrotic = pf.dicrotic;
    out.peak2_auto = pf.peak2;
    out.end_auto = pf.end;
    out.end = pf.end;
    out.notch_found = pf.notch_found;
    out.auto_placeholder = pf.placeholder;
}