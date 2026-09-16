#pragma once
/**
 * @file   subsample_refine.hpp
 * @brief  Sub-sample (floating-point) landmark refinement, per spec I-3.
 *         All functions take an integer SEED (from existing detectors) and
 *         return a double sub-sample position. Every method here operates
 *         on a small local window around the seed and is safe to call with
 *         a seed anywhere in a larger array (no absolute-column-0 searches
 *         -- the class of bug found and fixed in Phase B).
 */
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <functional>
#include "template_anchoring\curve_fit.hpp"

namespace subsample_refine {

    inline constexpr int kWindowHalfWidth = 7;     // 15-point window = seed +- 7
    inline constexpr double kResidualGuardFrac = 0.10;   // 10% of peak amplitude

    // ONE source of truth for each peak's Gaussian-weighting sigma, shared by
    // the DETECTOR (bestPeakExtremum call sites in feature_marks) and the focus
    // panel's redraw. Duplicating these as literals in both places let the drawn
    // fit and the placement fit diverge; centralize so they cannot.
    namespace peak_sigma {
        inline constexpr double R = 5.0;
        inline constexpr double P = 12.0;
        inline constexpr double T = 15.0;
        inline constexpr double Q = 5.0;   // Q-peak: no dedicated detector site; sensible default
    }

    // ---------------------------------------------------------------------
    // Exposed fit result for the peak finders. The finders previously returned
    // only the sub-sample position and discarded the polynomial; this carries
    // the fitted curve out so callers (the focus overlay, the CSV export) can
    // draw it and record its parameters.
    //
    // coeff is in ASCENDING powers of (t - seed), t in samples: value(t) =
    // coeff[0] + coeff[1]*(t-seed) + coeff[2]*(t-seed)^2 [+ coeff[3]*(t-seed)^3].
    // order is 2 (QUADRATIC) or 3 (CUBIC); unused coeffs are 0. A fallback or a
    // degenerate solve reports FIVE_POINT / SEED with order 0 and no usable
    // coeffs -- there is no polynomial to draw or serialize in that case.
    // ---------------------------------------------------------------------
    enum class CurveType { SEED, QUADRATIC, CUBIC, FIVE_POINT };

    struct ExtremumFit {
        double     position = -1.0;          // sub-sample extremum, absolute samples
        CurveType  type = CurveType::SEED;
        int        order = 0;                // 2, 3, or 0 (no polynomial)
        int        seed = 0;                 // absolute index where (t-seed)=0
        double     rss = std::numeric_limits<double>::quiet_NaN();
        int        npts = 0;                 // points the fit used (n, for BIC)
        std::array<double, 4> coeff{ 0,0,0,0 };  // ascending powers of (t-seed)

        // Evaluate the fitted polynomial at absolute sample x. NaN if no curve.
        double eval(double x) const {
            if (order < 2) return std::numeric_limits<double>::quiet_NaN();
            const double t = x - seed;
            double s = 0.0, tp = 1.0;
            for (int i = 0; i <= order; ++i) { s += coeff[i] * tp; tp *= t; }
            return s;
        }
    };

    // ---------------------------------------------------------------------
    // Shared small-matrix helpers (no external dependency beyond <cmath>).
    // ---------------------------------------------------------------------

    // Solve a symmetric normal-equations system A x = b via Gaussian
    // elimination with partial pivoting. n <= 4 in this file's usage.
    inline bool solveLinear(std::vector<std::vector<double>> A, std::vector<double> b,
        std::vector<double>& x) {
        const int n = static_cast<int>(b.size());
        x.assign(n, 0.0);
        for (int col = 0; col < n; ++col) {
            int piv = col;
            double best = std::fabs(A[col][col]);
            for (int r = col + 1; r < n; ++r)
                if (std::fabs(A[r][col]) > best) { best = std::fabs(A[r][col]); piv = r; }
            if (best < 1e-14) return false;
            std::swap(A[col], A[piv]); std::swap(b[col], b[piv]);
            for (int r = col + 1; r < n; ++r) {
                const double f = A[r][col] / A[col][col];
                for (int c = col; c < n; ++c) A[r][c] -= f * A[col][c];
                b[r] -= f * b[col];
            }
        }
        for (int row = n - 1; row >= 0; --row) {
            double s = b[row];
            for (int c = row + 1; c < n; ++c) s -= A[row][c] * x[c];
            x[row] = s / A[row][row];
        }
        return true;
    }

    struct WindowSamples {
        std::vector<double> t;    // local coordinate, seed at t=0
        std::vector<double> y;
        std::vector<double> w;    // Gaussian weight
        int seedAbs = 0;          // absolute index of t=0
        bool ok = false;
    };

    // Gather a (2*halfWidth+1)-point window around `seed`, Gaussian-weighted
    // by `sigma`, skipping NaN samples and out-of-bounds positions. Returns
    // ok=false if fewer than 5 usable samples remain (can't fit anything
    // meaningful).
    inline WindowSamples gatherWindow(const std::vector<double>& signal, int seed,
        double sigma, int halfWidth = kWindowHalfWidth) {
        WindowSamples ws; ws.seedAbs = seed;
        const int N = static_cast<int>(signal.size());
        if (seed < 0 || seed >= N) return ws;
        for (int d = -halfWidth; d <= halfWidth; ++d) {
            const int idx = seed + d;
            if (idx < 0 || idx >= N) continue;
            const double v = signal[idx];
            if (std::isnan(v)) continue;
            ws.t.push_back(static_cast<double>(d));
            ws.y.push_back(v);
            ws.w.push_back(std::exp(-(double)(d * d) / (2.0 * sigma * sigma)));
        }
        ws.ok = ws.t.size() >= 5;
        return ws;
    }

    // ---------------------------------------------------------------------
    // Residual guard: weighted RMS residual of a fit against its window,
    // as a fraction of the window's peak amplitude (max |y - mean(y)|).
    // ---------------------------------------------------------------------
    inline double weightedRmsResidualFrac(const WindowSamples& ws,
        const std::vector<double>& fitted) {
        double wsum = 0.0, wsq = 0.0;
        double mean = 0.0, wsumForMean = 0.0;
        for (size_t i = 0; i < ws.y.size(); ++i) { mean += ws.w[i] * ws.y[i]; wsumForMean += ws.w[i]; }
        mean = (wsumForMean > 0) ? mean / wsumForMean : 0.0;
        double peakAmp = 1e-12;
        for (double v : ws.y) peakAmp = std::max(peakAmp, std::fabs(v - mean));
        for (size_t i = 0; i < ws.y.size(); ++i) {
            const double r = ws.y[i] - fitted[i];
            wsum += ws.w[i] * r * r;
            wsq += ws.w[i];
        }
        const double rms = (wsq > 0) ? std::sqrt(wsum / wsq) : 0.0;
        return rms / peakAmp;
    }

    // 5-point unweighted parabola fallback: fit y = a t^2 + b t + c on the
    // Five-point parabola exposed as a drawable curve: the guaranteed fallback
    // for a peak, so even a broad/flat wave (P) always has a visible parabola.
    // Coeffs ascending in (t - seed): value = c0 + c1*(t-seed) + c2*(t-seed)^2.
    inline ExtremumFit fivePointParabolaFit(const std::vector<double>& signal, int seed) {
        ExtremumFit out;
        out.seed = seed;
        out.position = static_cast<double>(seed);
        const int N = static_cast<int>(signal.size());
        std::vector<double> t, y;
        for (int d = -2; d <= 2; ++d) {
            const int idx = seed + d;
            if (idx < 0 || idx >= N || std::isnan(signal[idx])) continue;
            t.push_back(static_cast<double>(d)); y.push_back(signal[idx]);
        }
        if (t.size() < 3) return out;   // SEED, no curve
        double S0 = 0, S1 = 0, S2 = 0, S3 = 0, S4 = 0, Y0 = 0, Y1 = 0, Y2 = 0;
        for (size_t i = 0; i < t.size(); ++i) {
            const double ti = t[i], ti2 = ti * ti;
            S0 += 1; S1 += ti; S2 += ti2; S3 += ti2 * ti; S4 += ti2 * ti2;
            Y0 += y[i]; Y1 += ti * y[i]; Y2 += ti2 * y[i];
        }
        std::vector<double> sol;
        if (!solveLinear({ {S4,S3,S2},{S3,S2,S1},{S2,S1,S0} }, { Y2,Y1,Y0 }, sol))
            return out;
        out.type = CurveType::FIVE_POINT;
        out.order = 2;
        out.npts = static_cast<int>(t.size());
        out.coeff = { sol[2], sol[1], sol[0], 0.0 };   // c0=c, c1=b, c2=a
        if (std::fabs(sol[0]) >= 1e-12) {
            const double tv = -sol[1] / (2.0 * sol[0]);
            out.position = static_cast<double>(seed) + std::clamp(tv, -2.0, 2.0);
        }
        return out;
    }

    // 5 points nearest the seed (unweighted), return the vertex (extremum)
    // or, for a linear/degenerate fit, the seed itself.
    inline double fivePointParabolaExtremum(const std::vector<double>& signal, int seed) {
        const int N = static_cast<int>(signal.size());
        std::vector<double> t, y;
        for (int d = -2; d <= 2; ++d) {
            const int idx = seed + d;
            if (idx < 0 || idx >= N || std::isnan(signal[idx])) continue;
            t.push_back(static_cast<double>(d)); y.push_back(signal[idx]);
        }
        if (t.size() < 3) return static_cast<double>(seed);
        // Normal equations for y = a t^2 + b t + c.
        double S0 = 0, S1 = 0, S2 = 0, S3 = 0, S4 = 0, Y0 = 0, Y1 = 0, Y2 = 0;
        for (size_t i = 0; i < t.size(); ++i) {
            const double ti = t[i], ti2 = ti * ti;
            S0 += 1; S1 += ti; S2 += ti2; S3 += ti2 * ti; S4 += ti2 * ti2;
            Y0 += y[i]; Y1 += ti * y[i]; Y2 += ti2 * y[i];
        }
        std::vector<double> sol;
        const bool ok = solveLinear({ {S4,S3,S2},{S3,S2,S1},{S2,S1,S0} }, { Y2,Y1,Y0 }, sol);
        if (!ok || std::fabs(sol[0]) < 1e-12) return static_cast<double>(seed);
        const double tVertex = -sol[1] / (2.0 * sol[0]);
        return static_cast<double>(seed) + std::clamp(tVertex, -2.0, 2.0);
    }

    inline ExtremumFit symmetricExtremumFit(const std::vector<double>& signal, int seed,
        double sigma, int halfWidth = kWindowHalfWidth) {
        ExtremumFit out;
        out.seed = seed;
        out.position = static_cast<double>(seed);

        const WindowSamples ws = gatherWindow(signal, seed, sigma, halfWidth);
        if (!ws.ok) return out;   // SEED, no polynomial

        // Weighted least squares for y = a t^2 + b t + c, t relative to seed.
        double S0 = 0, S1 = 0, S2 = 0, S3 = 0, S4 = 0, Y0 = 0, Y1 = 0, Y2 = 0;
        for (size_t i = 0; i < ws.t.size(); ++i) {
            const double wt = ws.w[i], ti = ws.t[i], ti2 = ti * ti;
            S0 += wt; S1 += wt * ti; S2 += wt * ti2; S3 += wt * ti2 * ti; S4 += wt * ti2 * ti2;
            Y0 += wt * ws.y[i]; Y1 += wt * ti * ws.y[i]; Y2 += wt * ti2 * ws.y[i];
        }
        std::vector<double> sol;
        const bool ok = solveLinear({ {S4,S3,S2},{S3,S2,S1},{S2,S1,S0} }, { Y2,Y1,Y0 }, sol);
        if (!ok || std::fabs(sol[0]) < 1e-12) {
            out.type = CurveType::FIVE_POINT;
            out.position = fivePointParabolaExtremum(signal, seed);
            return out;
        }

        std::vector<double> fitted(ws.t.size());
        double rss = 0.0;
        for (size_t i = 0; i < ws.t.size(); ++i) {
            fitted[i] = sol[0] * ws.t[i] * ws.t[i] + sol[1] * ws.t[i] + sol[2];
            const double e = ws.y[i] - fitted[i];
            rss += e * e;
        }
        if (weightedRmsResidualFrac(ws, fitted) > kResidualGuardFrac) {
            out.type = CurveType::FIVE_POINT;
            out.position = fivePointParabolaExtremum(signal, seed);
            return out;
        }

        const double tVertex = -sol[1] / (2.0 * sol[0]);
        out.type = CurveType::QUADRATIC;
        out.order = 2;
        out.rss = rss;
        out.npts = static_cast<int>(ws.t.size());
        out.coeff = { sol[2], sol[1], sol[0], 0.0 };   // c0=c, c1=b, c2=a
        out.position = static_cast<double>(seed)
            + std::clamp(tVertex, (double)-halfWidth, (double)halfWidth);
        return out;
    }

    // ---------------------------------------------------------------------
    // Symmetric extrema: Gaussian-weighted quadratic, vertex = extremum.
    // ---------------------------------------------------------------------
    inline double symmetricExtremum(const std::vector<double>& signal, int seed, double sigma) {
        return symmetricExtremumFit(signal, seed, sigma).position;
    }

    inline ExtremumFit asymmetricExtremumFit(const std::vector<double>& signal, int seed,
        double sigma, int halfWidth = kWindowHalfWidth) {
        ExtremumFit out;
        out.seed = seed;
        out.position = static_cast<double>(seed);

        const WindowSamples ws = gatherWindow(signal, seed, sigma, halfWidth);
        if (!ws.ok) return out;   // SEED, no polynomial

        auto fivePoint = [&]() {
            out.type = CurveType::FIVE_POINT;
            out.order = 0;
            out.position = fivePointParabolaExtremum(signal, seed);
            return out;
            };

        // Weighted least squares for y = a t^3 + b t^2 + c t + d, t rel. seed.
        double S[7] = { 0,0,0,0,0,0,0 };    // sum(w*t^k), k=0..6
        double Y[4] = { 0,0,0,0 };          // sum(w*t^k*y), k=0..3
        for (size_t i = 0; i < ws.t.size(); ++i) {
            const double wt = ws.w[i], ti = ws.t[i];
            double p = wt;
            for (int k = 0; k <= 6; ++k) { S[k] += p; p *= ti; }
            double q = wt;
            for (int k = 0; k <= 3; ++k) { Y[k] += q * ws.y[i]; q *= ti; }
        }
        std::vector<std::vector<double>> A = {
            {S[6],S[5],S[4],S[3]}, {S[5],S[4],S[3],S[2]},
            {S[4],S[3],S[2],S[1]}, {S[3],S[2],S[1],S[0]}
        };
        std::vector<double> sol;
        if (!solveLinear(A, { Y[3],Y[2],Y[1],Y[0] }, sol)) return fivePoint();

        std::vector<double> fitted(ws.t.size());
        double rss = 0.0;
        for (size_t i = 0; i < ws.t.size(); ++i) {
            const double ti = ws.t[i];
            fitted[i] = sol[0] * ti * ti * ti + sol[1] * ti * ti + sol[2] * ti + sol[3];
            const double e = ws.y[i] - fitted[i];
            rss += e * e;
        }
        if (weightedRmsResidualFrac(ws, fitted) > kResidualGuardFrac) return fivePoint();

        // dy/dt = 3a t^2 + 2b t + c = 0
        const double a = sol[0], b = sol[1], c = sol[2];
        double tBest;
        if (std::fabs(a) < 1e-12) {
            if (std::fabs(b) < 1e-12) return fivePoint();
            tBest = -c / (2.0 * b);
        }
        else {
            const double disc = 4.0 * b * b - 12.0 * a * c;
            if (disc < 0.0) return fivePoint();
            const double sq = std::sqrt(disc);
            const double r1 = (-2.0 * b + sq) / (6.0 * a);
            const double r2 = (-2.0 * b - sq) / (6.0 * a);
            const bool r1In = r1 >= -halfWidth && r1 <= halfWidth;
            const bool r2In = r2 >= -halfWidth && r2 <= halfWidth;
            if (r1In && r2In) tBest = (std::fabs(r1) < std::fabs(r2)) ? r1 : r2;
            else if (r1In) tBest = r1;
            else if (r2In) tBest = r2;
            else return fivePoint();
        }
        out.type = CurveType::CUBIC;
        out.order = 3;
        out.rss = rss;
        out.npts = static_cast<int>(ws.t.size());
        out.coeff = { sol[3], sol[2], sol[1], sol[0] };   // c0=d, c1=c, c2=b, c3=a
        out.position = static_cast<double>(seed)
            + std::clamp(tBest, (double)-halfWidth, (double)halfWidth);
        return out;
    }

    // ---------------------------------------------------------------------
    // Asymmetric extrema: cubic fit on Gaussian-weighted samples, solve
    // dy/dt = 3a t^2 + 2b t + c = 0 analytically; pick the root inside the
    // window closest to t=0 (the seed).
    // ---------------------------------------------------------------------
    inline double asymmetricExtremum(const std::vector<double>& signal, int seed, double sigma) {
        return asymmetricExtremumFit(signal, seed, sigma).position;
    }

    // ---------------------------------------------------------------------
    // Best-of quadratic vs cubic for a peak: fit BOTH over the same window and
    // keep the lower BIC. BIC = n*ln(RSS/n) + k*ln(n), with k the parameter
    // count (quadratic 3, cubic 4). Unlike raw RSS -- which a cubic, nesting the
    // quadratic, almost always wins -- BIC charges the cubic ln(n) for its extra
    // parameter, so the cubic is kept only when its lower RSS actually earns it
    // (a genuinely skewed peak). Ties and degenerate cubics fall to the
    // quadratic.
    // ---------------------------------------------------------------------
    inline ExtremumFit bestPeakExtremumFit(const std::vector<double>& signal, int seed,
        double sigma, int halfWidth = kWindowHalfWidth,
        curve_fit::PeakFitMode peakMode = curve_fit::PeakFitMode::Auto) {
        const ExtremumFit q = symmetricExtremumFit(signal, seed, sigma, halfWidth);
        const ExtremumFit c = asymmetricExtremumFit(signal, seed, sigma, halfWidth);
        // FORCED by the Fit-Peaks radio: return that model (fall back only if it
        // degenerated). Parabola = quadratic (symmetric), Cubic = asymmetric,
        // FivePoint = the 5-point fallback parabola.
        if (peakMode == curve_fit::PeakFitMode::FivePoint) return fivePointParabolaFit(signal, seed);
        if (peakMode == curve_fit::PeakFitMode::Parabola) return (q.order >= 2) ? q : c;
        if (peakMode == curve_fit::PeakFitMode::Cubic)    return (c.order >= 2) ? c : q;
        const bool qOk = q.order >= 2, cOk = c.order >= 2;
        if (!cOk) return q;               // cubic degenerate: quadratic (or its fallback)
        if (!qOk) return c;

        auto bic = [](double rss, int n, int k) -> double {
            if (n <= 0) return std::numeric_limits<double>::infinity();
            // Guard a (near-)perfect fit: ln(0) -> -inf would auto-win.
            const double mse = std::max(rss / n, 1e-300);
            return n * std::log(mse) + k * std::log(static_cast<double>(n));
            };
        const double bicQ = bic(q.rss, q.npts, 3);
        const double bicC = bic(c.rss, c.npts, 4);
        return (bicC < bicQ) ? c : q;     // strict: quadratic wins ties
    }
    inline double bestPeakExtremum(const std::vector<double>& signal, int seed, double sigma,
        curve_fit::PeakFitMode peakMode = curve_fit::PeakFitMode::Auto) {
        return bestPeakExtremumFit(signal, seed, sigma, kWindowHalfWidth, peakMode).position;
    }

    // ---------------------------------------------------------------------
    // Maximum-slope points: first derivative, then the SAME Gaussian-
    // weighted quadratic technique applied to the derivative (its vertex
    // is the point of maximum slope, same math as an extremum, just one
    // derivative order up).
    // ---------------------------------------------------------------------
    inline double maxSlopePoint(const std::vector<double>& signal, int seed, double sigma) {
        const int N = static_cast<int>(signal.size());
        std::vector<double> d1(N, std::numeric_limits<double>::quiet_NaN());
        for (int i = 1; i < N - 1; ++i)
            if (!std::isnan(signal[i - 1]) && !std::isnan(signal[i + 1]))
                d1[i] = (signal[i + 1] - signal[i - 1]) / 2.0;
        return symmetricExtremum(d1, seed, sigma);
    }

    // ---------------------------------------------------------------------
    // The three models transitionAnchor tests, exposed for display so the
    // focus panel can draw the EXACT fits that placed the mark (not a re-fit
    // over a different window). Each curve is SAMPLE-INDEXED already: the
    // caller evaluates curve[k](sampleIndex) directly -- the lo offset and the
    // 4x upsample mapping are baked into the closure. `winner` is the index the
    // BIC selector chose (0=piecewise, 1=sigmoid, 2=fractional). Populated only
    // when a caller passes candOut, so the detection hot path pays nothing.
    struct TransitionCandidates {
        std::function<double(double)> curve[4];   // 0=piecewise 1=sigmoid 2=fractional 3=cubic-spline
        double cross[4] = { -1.0, -1.0, -1.0, -1.0 };  // fiducial crossing per model (sample-indexed)
        int  winner = -1;
        bool valid = false;
    };

    // ---------------------------------------------------------------------
    // Transition onsets/offsets: locally upsample a 40-sample window from
    // its native rate to 4x via cubic interpolation, then fit-and-select
    // (Section 4.2 machinery, Phase A) on the upsampled window, returning a
    // sub-sample position in the ORIGINAL sample-rate coordinate.
    // ---------------------------------------------------------------------
    inline double transitionAnchor(const std::vector<double>& signal, int seed,
        double fraction, int windowSamples = 40,
        double externalBaseline = std::numeric_limits<double>::quiet_NaN(),
        double boundLo = -1.0, double boundHi = -1.0,
        TransitionCandidates* candOut = nullptr,
        curve_fit::FitMode mode = curve_fit::FitMode::Auto) {
        const int N = static_cast<int>(signal.size());
        // If the caller supplies explicit bounds (e.g. already correctly
        // one-sided, capped at a known extremum so the window can't cross
        // it), use those directly. Otherwise fall back to the seed +-
        // windowSamples/2 symmetric window, same as before. The symmetric
        // default can accidentally reach across a nearby extremum if the
        // seed sits close to one (the same class of bug found and fixed for
        // compute_s_end/compute_t_end/compute_q_onset in an earlier pass) --
        // callers that already know a safe one-sided range should supply it.
        // Bounds arrive as sub-sample doubles so callers can pass a landmark
        // straight through. The slice below is necessarily integer, so widen
        // outward -- floor the lower, ceil the upper -- and the requested span
        // is always covered rather than clipped.
        int lo, hi;
        if (boundLo >= 0.0 && boundHi >= 0.0 && boundHi > boundLo) {
            lo = std::max(0, static_cast<int>(std::floor(boundLo)));
            hi = std::min(N - 1, static_cast<int>(std::ceil(boundHi)));
        }
        else {
            const int half = windowSamples / 2;
            lo = std::max(0, seed - half);
            hi = std::min(N - 1, seed + half);
        }
        if (hi - lo < 4) return static_cast<double>(seed);

        std::vector<double> local(signal.begin() + lo, signal.begin() + hi + 1);
        for (double& v : local) if (std::isnan(v)) v = local.front();   // guard: no NaN into cubic interp

        // Cubic (Catmull-Rom) interpolation, 4x upsample.
        const int nIn = static_cast<int>(local.size());
        const int upsampleFactor = 4;
        const int nOut = (nIn - 1) * upsampleFactor + 1;
        std::vector<double> up(nOut);
        auto catmullRom = [&](double p0, double p1, double p2, double p3, double t) {
            return 0.5 * ((2 * p1) + (-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t * t + (-p0 + 3 * p1 - 3 * p2 + p3) * t * t * t);
            };
        for (int i = 0; i < nIn - 1; ++i) {
            const double p0 = local[std::max(0, i - 1)];
            const double p1 = local[i];
            const double p2 = local[i + 1];
            const double p3 = local[std::min(nIn - 1, i + 2)];
            for (int k = 0; k < upsampleFactor; ++k) {
                const double t = static_cast<double>(k) / upsampleFactor;
                up[i * upsampleFactor + k] = catmullRom(p0, p1, p2, p3, t);
            }
        }
        up[nOut - 1] = local.back();

        // seed's position on the upsampled grid.
        const int seedUp = (seed - lo) * upsampleFactor;
        // B: caller-supplied baseline (e.g. seed_all's own B_iso, so every
        // landmark shares ONE consistent isoelectric reference) if given,
        // else fall back to the local window's own first sample.
        const double B = std::isnan(externalBaseline) ? up.front() : externalBaseline;
        double E = B, bestDist = 0.0;
        for (double v : up) { const double dd = std::fabs(v - B); if (dd > bestDist) { bestDist = dd; E = v; } }

        auto fit = curve_fit::selectAnchorModel(up, 0, nOut - 1, mode);
        const double anchorUp = curve_fit::anchorAtFraction(fit, 0, nOut - 1, B, E, fraction);
        (void)seedUp;

        // Expose the three tested candidates as sample-indexed closures (only
        // when asked). Each model's evaluator works on the UPSAMPLED grid, so
        // the closure maps a sample x to that grid: up = (x - lo) * factor.
        if (candOut) {
            const int loc = lo, hic = hi, uf = upsampleFactor;
            auto mk = [loc, hic, uf](const curve_fit::FitResult& fr)
                -> std::function<double(double)> {
                if (!fr.f) return {};
                auto f = fr.f;
                return [f, loc, hic, uf](double sample) -> double {
                    // Only within the window it was fit on -- outside, a
                    // fractional model decays to its constant term and would
                    // paint a long flat tail that reads as "frac is flat".
                    if (sample < static_cast<double>(loc)
                        || sample > static_cast<double>(hic))
                        return std::numeric_limits<double>::quiet_NaN();
                    return f((sample - static_cast<double>(loc)) * static_cast<double>(uf));
                    };
                };
            const auto pw = curve_fit::fitPiecewiseLinear(up, 0, nOut - 1);
            const auto sg = curve_fit::fitSigmoid(up, 0, nOut - 1, pw);
            const auto fr = curve_fit::fitFractionalPolynomial(up, 0, nOut - 1);
            const auto sp = curve_fit::fitCubicSpline(up, 0, nOut - 1);
            candOut->curve[0] = mk(pw);
            candOut->curve[1] = mk(sg);
            candOut->curve[2] = mk(fr);
            candOut->curve[3] = mk(sp);
            // Each model's OWN fiducial crossing, mapped back to sample space --
            // the position that model would place the mark at. The focus dotted
            // line and (on selection) the mark itself use cross[winner].
            auto crossOf = [&](const curve_fit::FitResult& fr_) -> double {
                if (!fr_.f) return -1.0;
                const double au = curve_fit::anchorAtFraction(fr_, 0, nOut - 1, B, E, fraction);
                return static_cast<double>(lo) + au / static_cast<double>(upsampleFactor);
                };
            candOut->cross[0] = crossOf(pw);
            candOut->cross[1] = crossOf(sg);
            candOut->cross[2] = crossOf(fr);
            candOut->cross[3] = crossOf(sp);
            switch (fit.type) {
            case curve_fit::FitType::SIGMOID:      candOut->winner = 1; break;
            case curve_fit::FitType::FRACTIONAL:   candOut->winner = 2; break;
            case curve_fit::FitType::CUBIC_SPLINE: candOut->winner = 3; break;
            default:                                candOut->winner = 0; break;
            }
            candOut->valid = true;
        }

        // Back to original-rate coordinate: anchorUp is a position on the
        // 4x grid starting at local[0] == signal[lo].
        return static_cast<double>(lo) + anchorUp / upsampleFactor;
    }


    inline int cubicSplineNotch(const std::vector<double>&, int lo, int, int*) { return lo; }
}  // namespace subsample_refine