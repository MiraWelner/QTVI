// ppg_dicrotic.cpp -- DeepEntropyX Phase 2, E-5. See ppg_dicrotic.hpp for the
// method, the implemented/not-implemented table, and every deviation.

#include "ppg_dicrotic.hpp"
#include "stats_utils.hpp"                 // Biquad, butterworth_lowpass, filtfilt
#include "template_generation/ppg_pipeline.hpp"   // the one IEM and the one spline
#include "ppg_derivative.hpp"              // the one Savitzky-Golay bank

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace ppg_dicrotic {

    namespace {
        constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

        /// Solve M c = v in place for a small dense symmetric system, by
        /// Gauss-Jordan with partial pivoting. False when M is singular.
        bool solveInPlace(std::vector<std::vector<double>>& M, std::vector<double>& v) {
            const int n = static_cast<int>(v.size());
            for (int col = 0; col < n; ++col) {
                int piv = col;
                for (int r = col + 1; r < n; ++r)
                    if (std::fabs(M[r][col]) > std::fabs(M[piv][col])) piv = r;
                if (std::fabs(M[piv][col]) < 1e-300) return false;
                std::swap(M[piv], M[col]);
                std::swap(v[piv], v[col]);
                const double d = M[col][col];
                for (int c = col; c < n; ++c) M[col][c] /= d;
                v[col] /= d;
                for (int r = 0; r < n; ++r) {
                    if (r == col) continue;
                    const double f = M[r][col];
                    if (f == 0.0) continue;
                    for (int c = col; c < n; ++c) M[r][c] -= f * M[col][c];
                    v[r] -= f * v[col];
                }
            }
            return true;
        }
    }

    int DnConfig::sgWindow(double fs) const {
        if (sg_nc > 0) return (sg_nc % 2 == 1) ? sg_nc : sg_nc + 1;
        // ~0.098 s: 25 at 256 Hz, 49 at 500 Hz, per E-5.10.
        int n = static_cast<int>(std::lround(0.098 * fs));
        if (n % 2 == 0) ++n;
        // At least polyOrder+3, forced odd: a window at or below the
        // polynomial order is an exact (useless) fit.
        int floorN = sg_poly_order + 3;
        if (floorN % 2 == 0) ++floorN;
        return std::max(n, floorN);
    }

    // =========================================================================
    // Shared primitives
    // =========================================================================

    std::vector<double> bandLimit(const std::vector<double>& pulse, double fs,
        double cutoffHz)
    {
        constexpr int kOrder = 4;
        if (pulse.size() < 4 || !(fs > 0.0) || !(cutoffHz > 0.0)
            || cutoffHz >= fs / 2.0) return pulse;

        // NaN-safe: fill for the filter run, restore after. A single unfilled
        // NaN would contaminate every output sample, since filtfilt is IIR.
        std::vector<bool> wasNan(pulse.size(), false);
        std::vector<double> xi = pulse;
        std::size_t nValid = 0;
        for (std::size_t i = 0; i < xi.size(); ++i) {
            if (std::isnan(xi[i])) wasNan[i] = true;
            else ++nValid;
        }
        if (nValid == 0) return pulse;
        if (nValid < xi.size()) {
            double last = kNaN;
            for (std::size_t i = 0; i < xi.size(); ++i)
                if (!wasNan[i]) { last = xi[i]; break; }
            for (std::size_t i = 0; i < xi.size(); ++i) {
                if (wasNan[i]) xi[i] = last;
                else last = xi[i];
            }
            for (std::size_t i = 0; i < xi.size(); ) {
                if (!wasNan[i]) { ++i; continue; }
                std::size_t j = i;
                while (j < xi.size() && wasNan[j]) ++j;
                if (i > 0 && j < xi.size()) {
                    const double a = xi[i - 1], b = xi[j];
                    const double denom = static_cast<double>(j - (i - 1));
                    for (std::size_t k = i; k < j; ++k)
                        xi[k] = a + (b - a) * static_cast<double>(k - (i - 1)) / denom;
                }
                i = j;
            }
        }

        const std::vector<Biquad> sos = butterworth_lowpass(kOrder, cutoffHz, fs);
        std::vector<double> y = filtfilt(sos, xi);
        for (std::size_t i = 0; i < y.size(); ++i)
            if (wasNan[i]) y[i] = kNaN;
        return y;
    }

    std::vector<double> minMaxScale(const std::vector<double>& window) {
        double lo = std::numeric_limits<double>::infinity();
        double hi = -std::numeric_limits<double>::infinity();
        std::size_t nFinite = 0;
        for (const double v : window) {
            if (std::isnan(v)) continue;
            ++nFinite;
            if (v < lo) lo = v;
            if (v > hi) hi = v;
        }
        if (nFinite < 2 || !(hi > lo)) return window;

        std::vector<double> out(window.size());
        const double range = hi - lo;
        for (std::size_t i = 0; i < window.size(); ++i)
            out[i] = std::isnan(window[i]) ? window[i] : (window[i] - lo) / range;
        return out;
    }

    // =========================================================================
    // E-5.3 Envelope knots
    // =========================================================================
    //
    // CORRECT: sign change of d2 locates extrema of d1. See the header on why
    // the sign of d2 must not be used -- that version returns 0 knots on real
    // data, silently.
    std::vector<int> extremaOfFirstDerivative(const std::vector<double>& d1,
        const std::vector<double>& d2,
        std::vector<int>& maxima,
        std::vector<int>& minima)
    {
        (void)d1;   // not read in this pass; see the header note
        maxima.clear();
        minima.clear();
        std::vector<int> all;
        for (std::size_t i = 1; i < d2.size(); ++i) {
            if (std::isnan(d2[i - 1]) || std::isnan(d2[i])) continue;
            if ((d2[i - 1] > 0.0) == (d2[i] > 0.0)) continue;   // no sign change
            const int idx = (std::fabs(d2[i - 1]) < std::fabs(d2[i])) ? (int)i - 1 : (int)i;
            (d2[i - 1] > 0.0 ? maxima : minima).push_back(idx);   // + to - is a maximum of d1
            all.push_back(idx);
        }
        return all;
    }

    // E-5.1 steps 1 and 2: band-limit, then normalize. See the header for why
    // the order matters to the gate.
    std::vector<double> prepareWindow(const std::vector<double>& raw, double fs,
        const DnConfig& cfg)
    {
        return minMaxScale(bandLimit(raw, fs, cfg.ppg_prefilter_hz));
    }

    // =========================================================================
    // E-5.1 step 3: the artifact gate
    // =========================================================================
    ArtifactGateResult artifactGate(const std::vector<double>& window, double fs,
        const DnConfig& cfg, const RailBounds& rail)
    {
        ArtifactGateResult out;
        const std::size_t n = window.size();

        // ---- SUBSTITUTED criterion 1: dropout run / rail / constant --------
        int curRun = 0;
        double curVal = kNaN;
        for (std::size_t i = 0; i < n; ++i) {
            const double v = window[i];
            if (std::isnan(v)) { curRun = 0; curVal = kNaN; continue; }
            curRun = (curRun > 0 && v == curVal) ? curRun + 1 : 1;
            curVal = v;
            if (curRun > out.longestRunSamples) out.longestRunSamples = curRun;
            if (rail.known()
                && ((std::isfinite(rail.lo) && v <= rail.lo)
                    || (std::isfinite(rail.hi) && v >= rail.hi)))
                out.atRail = true;
        }
        if (fs > 0.0) {
            const double runMs = out.longestRunSamples * 1000.0 / fs;
            out.constantRun = runMs > cfg.dropout_run_ms;   // "longer than", strict
        }
        else {
            // No rate to convert ms -> samples; the only rate-free reading of
            // the rule is "the whole window is one run".
            out.constantRun = (n > 0 && out.longestRunSamples >= (int)n);
        }
        {
            // Independent of run length: a window that never varies except for
            // NaN gaps resets curRun at every gap and would dodge the check
            // above without this one.
            double lo = std::numeric_limits<double>::infinity();
            double hi = -std::numeric_limits<double>::infinity();
            std::size_t nFinite = 0;
            for (const double v : window) {
                if (std::isnan(v)) continue;
                ++nFinite;
                if (v < lo) lo = v;
                if (v > hi) hi = v;
            }
            out.constantWindow = (nFinite > 0 && !(hi > lo));
        }

        // ---- UNCHANGED criteria 2 and 3: peaks above the 75th percentile of
        // the window amplitude, via minMaxScale so the amplitude threshold and
        // the [0,1] scaling used elsewhere share one definition.
        {
            // Plain local-maximum test, NO prominence/width/distance
            // constraint: E-5.1 requires matching scipy.signal.find_peaks with
            // default settings, or the counts driving the gate differ from the
            // published ones.
            const std::vector<double> scaled = minMaxScale(window);
            int count = 0, all = 0;
            for (std::size_t i = 1; i + 1 < scaled.size(); ++i) {
                const double a = scaled[i - 1], b = scaled[i], c = scaled[i + 1];
                if (std::isnan(a) || std::isnan(b) || std::isnan(c)) continue;
                if (!(b > a && b > c)) continue;
                ++all;
                if (b > 0.75) ++count;
            }
            out.nPeaksAboveP75 = count;
            out.nPeaksAll = all;
            out.tooFewPeaks = count <= cfg.dn_gate_peak_min;
            out.tooManyPeaks = count > cfg.dn_gate_peak_max;
        }

        out.reject = out.constantRun || out.atRail || out.constantWindow
            || out.tooFewPeaks || out.tooManyPeaks;
        return out;
    }

    // =========================================================================
    // E-5.8 Prominence
    // =========================================================================
    double valleyProminence(const std::vector<double>& nonStationary, int valley) {
        const int n = static_cast<int>(nonStationary.size());
        if (valley < 1 || valley >= n - 1) return kNaN;
        const double v = nonStationary[valley];
        if (!std::isfinite(v)) return kNaN;

        // Walk out each way to the highest point reached before the series
        // falls back below this valley (the key saddle on that side).
        auto sideMax = [&](int step) {
            double best = v;
            for (int i = valley + step; i >= 0 && i < n; i += step) {
                const double x = nonStationary[i];
                if (!std::isfinite(x)) break;
                if (x < v) break;            // a deeper basin: this side closes
                if (x > best) best = x;
            }
            return best;
            };
        const double l = sideMax(-1), r = sideMax(+1);
        return std::min(l, r) - v;
    }

    int valleyWidthSamples(const std::vector<double>& nonStationary, int valley) {
        const int n = static_cast<int>(nonStationary.size());
        if (valley < 1 || valley >= n - 1) return 0;
        const double v = nonStationary[valley];
        if (!std::isfinite(v)) return 0;

        // Out each way to the local maximum bounding this valley, then back to
        // where the series first reached half that height: the half-prominence
        // crossing. Half-prominence rather than, say, the zero crossing,
        // because a notch riding on the broad diastolic bowl never returns to
        // zero on its right-hand side -- a zero-referenced width would measure
        // the bowl on every pulse and separate nothing.
        auto side = [&](int step) {
            double best = v;
            int at = valley;
            for (int i = valley + step; i >= 0 && i < n; i += step) {
                const double x = nonStationary[i];
                if (!std::isfinite(x) || x < v) break;
                if (x > best) { best = x; at = i; }
            }
            if (at == valley) return 0;
            const double half = v + 0.5 * (best - v);
            for (int i = valley; i != at; i += step)
                if (nonStationary[i] >= half) return std::abs(i - valley);
            return std::abs(at - valley);
            };
        return side(-1) + side(+1);
    }

    // =========================================================================
    // E-5.5 Tier 1 localization
    // =========================================================================
    DnResult tier1Dn(const std::vector<double>& nonStationary, double fs,
        int sysPeak, int diastolicEnd, const DnConfig& cfg)
    {
        DnResult out;   // UNDETECTED, index -1
        const int n = static_cast<int>(nonStationary.size());
        if (n < 5 || !(fs > 0.0) || sysPeak < 0) return out;

        // At least dn_window_lo_ms after the systolic peak; upper bound is the
        // DIASTOLIC ENDPOINT, not a fraction of RR.
        const int lo = sysPeak + static_cast<int>(std::lround(cfg.dn_window_lo_ms * 1e-3 * fs));
        const int hi = std::min(n - 1, (diastolicEnd > sysPeak) ? diastolicEnd : n - 1);
        if (hi - lo < 3 || lo < 1) return out;

        // THE FIRST valley meeting the conditions, not the deepest.
        for (int i = std::max(lo, 1); i < hi; ++i) {
            const double a = nonStationary[i - 1], b = nonStationary[i], c = nonStationary[i + 1];
            if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c)) continue;
            if (!(b <= a && b <= c)) continue;   // not a valley
            if (!(b < 0.0)) continue;            // source condition 2

            // NO THIRD CONDITION. E-5.5: the first valley meeting both is
            // taken, and there is no prominence gate.
            out.index = i;
            out.tier = DnResult::IEM;
            out.prominence = valleyProminence(nonStationary, i);
            // Parabolic vertex through the three samples, for the sub-sample
            // position. Skipped when the three are collinear (flat valley
            // floor), where a vertex is not defined.
            const double den = (a - 2.0 * b + c);
            if (std::fabs(den) > 1e-300) {
                const double delta = 0.5 * (a - c) / den;
                if (std::fabs(delta) <= 1.0) out.subSample = static_cast<double>(i) + delta;
            }
            return out;
        }
        return out;   // Tier 3: no qualifying valley
    }

    // =========================================================================
    // E-5.6 Tier 2: flow reconstruction
    // =========================================================================
    double fitTauWindow(const std::vector<double>& raw, double fs,
        const std::vector<int>& sysPeaks,
        const std::vector<int>& cycleEnds)
    {
        const int n = static_cast<int>(raw.size());
        if (n < 8 || !(fs > 0.0) || sysPeaks.empty()) return kNaN;

        const double pmin = *std::min_element(raw.begin(), raw.end());

        // POOLED DIASTOLIC SEGMENTS. The source says "pooled diastolic
        // segments" without fixing where diastole starts; the notch is not
        // known yet, so the LATTER HALF of each systolic-peak-to-cycle-end
        // interval is used, which is safely past the notch on a normal pulse.
        // Named here as a choice, not a published value.
        std::vector<double> t, v;
        const std::size_t nCyc = std::min(sysPeaks.size(), cycleEnds.size());
        for (std::size_t c = 0; c < nCyc; ++c) {
            const int a = sysPeaks[c], b = cycleEnds[c];
            if (a < 0 || b <= a || b >= n) continue;
            const int start = a + (b - a) / 2;
            for (int i = start; i <= b; ++i) {
                const double y = raw[i] - pmin;
                if (!std::isfinite(y) || y <= 0.0) continue;
                t.push_back(static_cast<double>(i - start) / fs);
                v.push_back(y);
            }
        }
        if (t.size() < 8) return kNaN;

        // Log-linear seed: log(v) = log(A) - t/tau.
        double st = 0, sl = 0, stt = 0, stl = 0;
        const double m = static_cast<double>(t.size());
        for (std::size_t i = 0; i < t.size(); ++i) {
            const double lv = std::log(v[i]);
            st += t[i]; sl += lv; stt += t[i] * t[i]; stl += t[i] * lv;
        }
        const double den = m * stt - st * st;
        if (std::fabs(den) < 1e-300) return kNaN;
        const double slope = (m * stl - st * sl) / den;
        if (!(slope < 0.0)) return kNaN;             // not a decay
        double tau = -1.0 / slope;
        double A = std::exp((sl - slope * st) / m);

        // Two or three Gauss-Newton iterations on the UNTRANSFORMED
        // exponential: the log transform makes the noise multiplicative and
        // over-weights the low-amplitude tail.
        for (int it = 0; it < 3; ++it) {
            double jaa = 0, jat = 0, jtt = 0, ra = 0, rt = 0;
            for (std::size_t i = 0; i < t.size(); ++i) {
                const double e = std::exp(-t[i] / tau);
                const double res = v[i] - A * e;
                const double dA = -e;
                const double dT = -A * t[i] * e / (tau * tau);
                jaa += dA * dA; jat += dA * dT; jtt += dT * dT;
                ra += dA * res; rt += dT * res;
            }
            std::vector<std::vector<double>> M{ { jaa, jat }, { jat, jtt } };
            std::vector<double> g{ -ra, -rt };
            if (!solveInPlace(M, g)) break;
            A += g[0];
            tau += g[1];
            if (!(tau > 0.0) || !std::isfinite(tau) || !std::isfinite(A)) return kNaN;
        }
        // Physiological guard: a time constant outside this range is a failed
        // fit, not a measurement (Hoeksel's are near 550 ms).
        if (!(tau > 0.05 && tau < 5.0)) return kNaN;
        return tau;
    }

    DnResult flowReconstructionDn(const std::vector<double>& raw, double fs,
        int sysPeak, double tauWindow, const DnConfig& cfg)
    {
        DnResult out;
        const int n = static_cast<int>(raw.size());
        if (n < 8 || !(fs > 0.0) || sysPeak < 0 || !(tauWindow > 0.0)) return out;

        const ppg_deriv::DerivBank D = ppg_deriv::buildDerivatives(
            raw, fs, cfg.sgHalfWidth(fs), cfg.sgHalfWidth(fs), cfg.sg_poly_order);
        if (D.d1.empty()) return out;
        const double pmin = *std::min_element(raw.begin(), raw.end());

        // g(t) is proportional to flow: g = dP/dt + P/tau. Compliance scales g
        // but cannot move its extremum, so only tau is needed.
        std::vector<double> g(n, kNaN);
        for (int i = 0; i < n; ++i) {
            if (!std::isfinite(D.d1[i]) || !std::isfinite(raw[i])) continue;
            g[i] = D.d1[i] + (raw[i] - pmin) / tauWindow;
        }

        const int lo = sysPeak + static_cast<int>(std::lround(cfg.dn_window_lo_ms * 1e-3 * fs));
        const int hi = std::min(n, lo + static_cast<int>(std::lround(0.60 * fs)));
        if (hi - lo < 5 || lo < 0) return out;

        // MINIMUM OF THE FIRST NEGATIVE DIP. Plain minimum, not minimum of
        // |g|: min|g| selects where flow is near zero, not where backflow is
        // greatest. And NOT a zero-crossing -- see the header.
        int best = -1;
        for (int i = lo; i < hi; ++i) {
            if (!std::isfinite(g[i])) continue;
            if (best < 0 || g[i] < g[best]) best = i;
        }
        if (best < 0) return out;

        out.index = best;
        out.tier = DnResult::FLOW;
        return out;
    }

    // =========================================================================
    // E-5.7 The windowed pass and the join contract
    // =========================================================================
    WindowedPass detectDicroticNotchWindowed(const std::vector<double>& signal,
        double fs,
        const std::vector<int>& sysPeaks,
        const std::vector<int>& cycleEnds,
        const DnConfig& cfg,
        const RailBounds& rail)
    {
        WindowedPass P;
        const int N = static_cast<int>(signal.size());
        const int nCyc = static_cast<int>(std::min(sysPeaks.size(), cycleEnds.size()));
        P.perCycle.assign(nCyc, DnResult{});
        P.coveredCentrally.assign(nCyc, 0);
        if (N < 16 || !(fs > 0.0) || nCyc == 0) return P;

        const int winLen = std::max(16, static_cast<int>(std::lround(cfg.iem_window_s * fs)));
        const double ov = (cfg.dn_window_overlap > 0.0 && cfg.dn_window_overlap < 1.0)
            ? cfg.dn_window_overlap : 0.5;
        const int advance = std::max(1, static_cast<int>(std::lround(winLen * ov)));

        // One candidate per (cycle, covering window). Kept rather than reduced
        // on the fly because the join needs BOTH of a cycle's estimates: the
        // nearest-center one to report, and the other one to measure
        // disagreement against.
        struct Cand {
            DnResult r;
            double distToCenter = 0.0;
            bool central = false;
            bool rejected = false;
        };
        std::vector<std::vector<Cand>> cand(nCyc);

        for (int w0 = 0; w0 < N; w0 += advance) {
            const int w1 = std::min(N, w0 + winLen);
            const int L = w1 - w0;
            if (L < 16) break;
            ++P.nWindows;

            const double center = 0.5 * (w0 + w1);
            const int cenLo = w0 + L / 4;          // middle half: see the header
            const int cenHi = w0 + (3 * L) / 4;

            // Which cycles this window can speak for at all: fully inside it.
            std::vector<int> mine;
            for (int c = 0; c < nCyc; ++c)
                if (sysPeaks[c] >= w0 && cycleEnds[c] <= w1
                    && cycleEnds[c] > sysPeaks[c] && sysPeaks[c] >= 0)
                    mine.push_back(c);
            if (mine.empty()) { if (w1 >= N) break; else continue; }

            const std::vector<double> raw(signal.begin() + w0, signal.begin() + w1);

            // ---- E-5.1 STEPS 1-3, IN ORDER --------------------------------
            //
            // ORDER MATTERS, and getting it wrong rejects most of a real
            // record. The gate's peak-count criteria (3 to 10 peaks above 75%
            // of the window's range) presuppose a clean signal: a 4-second
            // window holds about 4.5 cycles at a normal rate, so the band is
            // really "1 to 2 peaks per cycle". Counted on the RAW window,
            // high-frequency ripple on each systolic peak contributes extra
            // local maxima and the count inflates past the ceiling.
            //
            // Measured on MESA 3014843 (500 Hz, median RR 883 ms): counting on
            // the raw window gives 6 to 17 peaks and rejects 9 of 14 windows;
            // counting on the 16 Hz band-limited window gives 5 to 11 and
            // rejects 1. The rejected windows were not artefactual -- the
            // substituted criterion (dropout run / rail / constant) never
            // fired on any of them.
            const std::vector<double> seg = prepareWindow(raw, fs, cfg);
            const ArtifactGateResult gate = artifactGate(seg, fs, cfg, rail);
            if (gate.reject) {
                ++P.nWindowsRejected;
                for (const int c : mine) {
                    Cand k;
                    k.rejected = true;
                    k.central = (sysPeaks[c] >= cenLo && sysPeaks[c] <= cenHi);
                    k.distToCenter = std::fabs(sysPeaks[c] - center);
                    cand[c].push_back(k);
                }
                if (w1 >= N) break;
                continue;
            }

            // ---- one IEM per window, shared by every cycle in it ----------
            // Already band-limited AND min-max normalized by prepareWindow:
            // iem_beta is an absolute energy threshold (E-5.1 step 2).
            const std::vector<double>& filt = seg;
            const ppg_pipeline::IemEnvelope iem =
                ppg_pipeline::iemEnvelope(filt, fs, cfg.iem_max_iter, cfg.iem_beta);

            // tau FITTED ACROSS THE WHOLE WINDOW (E-5.6), from every cycle in
            // it -- not per pulse. Computed once, lazily, and only if some
            // cycle actually falls through to Tier 2.
            double tau = 0.0;
            bool tauTried = false;

            for (const int c : mine) {
                Cand k;
                k.central = (sysPeaks[c] >= cenLo && sysPeaks[c] <= cenHi);
                k.distToCenter = std::fabs(sysPeaks[c] - center);

                DnResult r = tier1Dn(iem.nonStationary, fs,
                    sysPeaks[c] - w0, cycleEnds[c] - w0, cfg);
                if (!r.found()) {
                    if (!tauTried) {
                        tauTried = true;
                        std::vector<int> wp, we;
                        for (const int q : mine) { wp.push_back(sysPeaks[q] - w0); we.push_back(cycleEnds[q] - w0); }
                        tau = fitTauWindow(filt, fs, wp, we);
                    }
                    if (std::isfinite(tau))
                        r = flowReconstructionDn(filt, fs, sysPeaks[c] - w0, tau, cfg);
                }
                if (r.index >= 0) {
                    r.index += w0;                                   // -> signal coords
                    if (std::isfinite(r.subSample)) r.subSample += w0;
                }
                k.r = r;
                cand[c].push_back(k);
            }
            if (w1 >= N) break;
        }

        // ---- the join ----------------------------------------------------
        for (int c = 0; c < nCyc; ++c) {
            // NEAREST-CENTER WINS, among the windows that covered this cycle
            // CENTRALLY. A peripheral candidate is never promoted: that is
            // what "discard boundary-only cycles" means.
            const Cand* best = nullptr;
            bool centralRejected = false;
            for (const Cand& k : cand[c]) {
                if (!k.central) continue;
                P.coveredCentrally[c] = 1;
                if (k.rejected) { centralRejected = true; continue; }
                if (!best || k.distToCenter < best->distToCenter) best = &k;
            }
            if (!P.coveredCentrally[c]) { ++P.nBoundaryOnly; continue; }

            if (!best) {
                // Covered centrally only by window(s) the gate rejected.
                // windowRejected, DISTINCT from UNDETECTED (E-5.9): these
                // cycles were never examined, so they are not evidence about
                // notch detectability and must not be counted as if they were.
                if (centralRejected) {
                    P.perCycle[c].windowRejected = true;
                    ++P.nCyclesWindowRejected;
                }
                else ++P.nUndetected;
                continue;
            }
            P.perCycle[c] = best->r;

            // ARCHIVE THE DISAGREEMENT: |t_A - t_B| against the other
            // window's estimate of the SAME cycle, in samples. Only between
            // two real detections -- a difference against a non-detection is
            // not a position difference.
            if (best->r.found()) {
                const Cand* other = nullptr;
                for (const Cand& k : cand[c]) {
                    if (&k == best || k.rejected || !k.r.found()) continue;
                    if (!other || k.distToCenter < other->distToCenter) other = &k;
                }
                if (other)
                    P.perCycle[c].confidence = std::fabs(best->r.position() - other->r.position());
            }

            if (P.perCycle[c].tier == DnResult::IEM) ++P.nTier1;
            else if (P.perCycle[c].tier == DnResult::FLOW) ++P.nTier2;
            else ++P.nUndetected;
        }
        return P;
    }

}  // namespace ppg_dicrotic