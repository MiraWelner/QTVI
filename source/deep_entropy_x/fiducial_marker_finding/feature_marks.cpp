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
double FeatureMarks::find_t_peak(const std::vector<double>& v, double bracketSEnd, double bracketTEnd, curve_fit::PeakFitMode peakMode) {
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
    if (v[best] < B) for (auto& x : u) x = -x;
    return std::clamp(upsample_for_fit::find_peak(u, best, upsample_for_fit::peak_sigma::T, upsample_for_fit::peak_halfwidth::T, peakMode), static_cast<double>(lo), static_cast<double>(hi));
}

// Local deflection direction, same reasoning as find_t_peak: an inverted P is a
// real waveform on an upright lead. No sgn parameter.
double FeatureMarks::find_p_peak(const std::vector<double>& v, double loIn, double hiIn, double fs, curve_fit::PeakFitMode peakMode)
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
    if (v[best] < B) for (double& x : u) x = -x;
    return std::clamp(upsample_for_fit::find_peak(u, best, upsample_for_fit::peak_sigma::P, upsample_for_fit::peak_halfwidth::P, peakMode), static_cast<double>(fFin), static_cast<double>(lFin));
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
    if (hi - lo < 4) return cld(sPeakD);
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
    if (hi - lo < 4) return cld(sPeakD);
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

double FeatureMarks::find_t_end(const std::vector<double>& v, double fs, int r_col, double j_point,
    upsample_for_fit::TransitionCandidates* candOut, curve_fit::FitMode mode) {
    const int N = static_cast<int>(v.size());
    if (N < 4 || r_col < 0 || r_col >= N) return -1.0;
    auto cl = [&](int i) { return std::clamp(i, 0, N - 1); };
    auto cld = [&](double d) { return std::clamp(d, 0.0, static_cast<double>(N - 1)); };

    // Window: [T-begin + 100 ms, T-begin + 350 ms], so the landmark chain runs
    // J-point -> T-begin -> T-end, each bounding the next.
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

    if (std::isnan(B) || std::isnan(E)) return -1.0;
    auto fit = curve_fit::selectBestFit(v, lo, hi);
    const double af = curve_fit::anchorAtFraction(fit, lo, hi, B, E, 0.02);
    if (!std::isfinite(af)) return -1.0;
    const int seed = std::clamp(static_cast<int>(std::round(af)), 0, N - 1);
    const double te = upsample_for_fit::upsample_transition_and_curve_fit(v, seed, 0.02, 40, B, lo, hi, candOut, mode);
    if (!std::isfinite(te)) return -1.0;
    return cld(te);
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

FeatureMarks::ReactiveEcg FeatureMarks::reactive_ecg(const std::vector<double>& ecg, double p_begin, double q_onset, double s_end, double t_end, double sampleRate,
    curve_fit::PeakFitMode peakMode)
{
    ReactiveEcg r;
    if (static_cast<int>(ecg.size()) < 3) return r;
    // Both reactive peaks are locally-polarised finders, so neither needs sgn.
    r.p_peak = find_p_peak(ecg, p_begin, q_onset, sampleRate, peakMode);
    r.t_peak = find_t_peak(ecg, s_end, t_end, peakMode);
    return r;
}

FeatureMarks::ReactivePpg FeatureMarks::reactive_ppg(const std::vector<double>& ppg, double onset, double peak, double dicrotic, double end)
{
    ReactivePpg r;
    if (static_cast<int>(ppg.size()) < 3) return r;

    // THE BRACKETS ROUND, THE RESULTS DO NOT. The bars arrive sub-sample, but
    // amplitude_crossing and crossing_at_level take an integer search grid and
    // return an interpolated position inside it, so a bracket half a sample
    // either way does not move the answer.
    auto win = [](double x) {
        return (x < 0.0) ? -1 : static_cast<int>(std::lround(x));
        };
    const int iOnset = win(onset), iPeak = win(peak), iEnd = win(end);

    if (iOnset >= 0 && iPeak > iOnset)
        r.t50 = amplitude_crossing(ppg, iOnset, iPeak, 0.50);
    if (iPeak >= 0 && iEnd > iPeak) {
        r.t80 = amplitude_crossing(ppg, iPeak, iEnd, 0.80);
        // T80_rise / PW80 at t80's OWN absolute level (see
        // detect_ppg_fiducials for the rationale): upslope crossing of the
        // same value. The two amplitudes come from the UNROUNDED bar positions
        // -- a level is a measurement, not a window, and sample_at interpolates.
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

FeatureMarks::PpgFiducials FeatureMarks::detect_ppg_fiducials(const std::vector<double>& v, int W, double ppgRate, double heightMeters,
    curve_fit::PeakFitMode peakMode)
{
    PpgFiducials g;
    const int N = static_cast<int>(v.size());
    if (N < 3) { fprintf(stderr, "[ppg] bail: N=%d\n", N); return g; }
    const int Wc = std::clamp(W, 2, N);   // visible window; nothing is ever placed past Wc-1
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

    // ---- THE SAME CONTEST THE ECG PEAKS RUN ----------------------------
    //
    // peakCandidates runs bestPeakExtremumFit -- guarded quadratic vs guarded
    // cubic, lower BIC wins, 5-point parabola when both trip the residual
    // guard -- and `placement` is its answer. That is exactly what an ECG peak
    // takes (BinPlotWidget::detectedLandmarks: `if (pc.valid &&
    // pc.placement >= 0.0) fid = pc.placement;`), so the two channels now
    // differ only in their seed, sigma and half-width.
    //
    // WHAT THIS REPLACED, AND WHY IT HAD TO GO. The placement used to be
    // draw[0].position with `winner = 0` written over the contest's choice:
    // draw[] are UNGUARDED fits, so the pulse peak was the raw weighted
    // quadratic's vertex whether or not it fitted, and the guard could not
    // demote it because nothing read the guarded result. Two consequences,
    // both visible in the focus panel: a quadratic that failed the 10%
    // residual test still placed the mark (the "Quadratic (rejected)" in
    // green), and the Fit-Peaks radio did nothing on a pulse.
    //
    // THIS MOVES DETECTED PULSE PEAKS. That is the intended effect -- a
    // skewed pulse apex is what the cubic is for, and a misfitted one is what
    // the 5-point fallback is for -- but stored *_auto columns and the pulse
    // CSVs will differ from a build before this change.
    //
    // THE WINDOW IS A DURATION, resolved against this channel's own rate --
    // see upsample_for_fit::pulse_window for why a fixed sample count was
    // fitting a tenth of a pulse on one recording and half of systole on
    // another. The focus panel resolves the same two numbers from the same
    // helpers, so the drawn curve is still the fit that placed the mark.
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

    // A systolic peak with no room for a foot before it is a head fragment,
    // not a pulse. Bail rather than pile every landmark at sample 0.
    if (g.peak < 3) { return g; }

    // Foot and end of cycle: the cubic's vertex on the rising / falling
    // shoulder. Both landmarks use one lambda: they had two, identical down to
    // the constants, which is two places for the pulse foot's window to drift.
    // winner = 1 names the cubic as the placing model.
    //
    // ---- WHY THIS RETURNS -1 AND DOES NOT CLAMP ------------------------
    //
    // cld() IS NOT SAFE ON A FIT RESULT. peakCandidates bails early --
    // n < 5, seed out of range, halfWidth < 3 -- and returns a
    // DEFAULT-CONSTRUCTED PeakCandidates: valid = false, draw[].position = -1.
    // Running that -1 through cld() is std::clamp(-1.0, 0.0, Wc-1), which is
    // 0.0. So a fit that never ran reported the foot AT SAMPLE 0, and sample 0
    // is inside the leading NaN pad every template carries (construction pads
    // `pad` seconds before the first R). That is not a cosmetic error:
    //
    //   onset = 0  ->  footY = sample_y(tmpl, 0) = NaN
    //              ->  calculate_perfusion_index returns NaN for EVERY sample
    //              ->  the whole normalized pulse trace is NaN
    //              ->  nothing is drawn, the pulse axis falls back to its
    //                  0..1 default, and drawFeatureGlyphs paints every
    //                  fiducial at the axis floor
    //
    // -- i.e. an empty panel with a row of X's along the bottom, on a slot
    // holding hundreds of clean beats. And it was STICKY: BankPulseMarkerSet::
    // isUnset() tests onset < 0, so the laundered 0 counted as a real
    // detection and the lazy re-seed in showPage never fired again.
    //
    // The peak never had this failure because it is validated immediately
    // after placement (`if (g.peak < 3) return g;`). The trough path lost that
    // guard when the two lambdas were consolidated into this one. -1 is the
    // sentinel every consumer of PpgFiducials already handles as "not found";
    // 0 is a position, and a wrong one.
    //
    // THE CONTEST HERE TOO, for the same reason as the peak above. This used
    // to read draw[1].position -- the unguarded CUBIC -- and then write
    // `winner = 1` over the contest's choice, so a trough was placed by the
    // cubic whether or not the cubic fitted it, and "Cubic (rejected)" in the
    // focus panel was a guard the detector had already ignored.
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
    g.dicrotic = cld(g.peak + 0.12 * ppgRate);
    if (g.end > g.peak && g.dicrotic >= g.end)
        g.dicrotic = cld(g.peak + 0.5 * (g.end - g.peak));   // midway peak->end
    g.notch_found = false;
    g.dn_tier = 0;
    g.dn_confidence = 0.0;

    g.t80 = amplitude_crossing(v, iFloor(g.peak), iCeil(g.end), 0.80);
    if (g.t80 < 0) g.t80 = cld(0.5 * (g.peak + g.end));
    g.peak2 = detect_ppg_peak2(v, iFloor(g.dicrotic), g.t80, iFloor(g.end));
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
    const LeadPolarity& pol, double heightMeters,
    curve_fit::FitMode fitMode, curve_fit::PeakFitMode peakMode) {
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
            const FeatureMarks::PpgFiducials pf = FeatureMarks::detect_ppg_fiducials(trace, static_cast<int>(trace.size()), sampleRate, NAN, peakMode);
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
    return out;
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
    curve_fit::PeakFitMode peakMode) {
    out = tbank::BankPulseMarkerSet{};
    const int W = static_cast<int>(tmpl.size());
    if (W < 3 || ppgRate <= 0.0) return;

    const PpgFiducials pf = detect_ppg_fiducials(tmpl, W, ppgRate, heightMeters, peakMode);
    out.onset_auto = pf.onset;
    out.onset = pf.onset;
    out.peak_auto = pf.peak;
    out.dicrotic_auto = pf.dicrotic;
    out.dicrotic = pf.dicrotic;
    out.peak2_auto = pf.peak2;
    out.end_auto = pf.end;
    out.end = pf.end;
    out.notch_found = pf.notch_found;
}