/**
 * @file   curve_fit.hpp
 * @brief  Model-selection anchor placement for ECG transition landmarks.
 *
 *         Fits three candidate models (piecewise-linear, sigmoid, fractional
 *         polynomial) to the signal region around a landmark, selects the
 *         best by BIC, and reads the anchor position as a fractional crossing
 *         of the fitted curve between baseline and extremum.
 *
 *         Used by feature_marks.cpp's compute_q_onset, compute_s_end,
 *         compute_t_end, and the P-onset helper. Amplitude landmarks
 *         (R, S, T-peak, P-peak) stay on their existing detectors.
 *
 *         Design notes (from spec Sections 3.7 and 4.2):
 *           - Piecewise-linear and sigmoid are always tried.
 *           - Fractional polynomial is the escalation model, tried only
 *             when neither simple model fits well (gated by POOR_FIT_THRESH).
 *           - The POOR_FIT_THRESH below is a starting estimate; calibrate
 *             empirically from clean expert-marked templates.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-07-26
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <vector>

namespace curve_fit {

    // =========================================================================
    // Tuning constants
    // =========================================================================

    // Raw RSS threshold for escalating to the fractional polynomial. If the
    // better of piecewise-linear and sigmoid has RSS above POOR_FIT_FACTOR * n,
    // the fractional polynomial is tried. Tune from clean-template residuals.
    //
    // Per spec: "const double POOR = 0.05 * n; // tune from clean-template residuals"
    inline constexpr double POOR_FIT_FACTOR = 0.05;

    // Sigmoid fitter iteration budget. 50 iterations is generous for a
    // 4-parameter fit over a 30-100 sample window; early-exit on convergence
    // means most calls finish in 10-20.
    inline constexpr int    SIGMOID_MAX_ITER = 50;
    inline constexpr double SIGMOID_CONV_TOL = 1e-8;

    // =========================================================================
    // FitResult
    // =========================================================================

    // Which candidate model a FitResult came from. Recorded so downstream
    // consumers (e.g. the boundary training log) can label each fit.
    enum class FitType { LINEAR, SIGMOID, FRACTIONAL, CUBIC_SPLINE, FLAT };

    // Forced model selection, driven by the on/offset radio group. Auto = the
    // BIC contest across all candidates (the previous behaviour); any other
    // value returns exactly that model so the anchor is placed from it.
    enum class FitMode { Auto, Linear, CubicSpline, Sigmoid, FracPoly };
    // Peak radio group (Fit Peaks): Auto = BIC quad-vs-cubic; else forced.
    enum class PeakFitMode { Auto, Cubic, Parabola, FivePoint };

    struct FitResult {
        double rss = std::numeric_limits<double>::infinity();
        int    nparams = 0;
        FitType type = FitType::FLAT;
        std::function<double(double)> f;   // evaluator: f(sample_index) -> fitted value
        // Raw fitted parameters, exposed for serialization. Meaning by type:
        //   LINEAR     : {m1, c1, m2, c2, breakpoint}
        //   SIGMOID    : {a, k, t0, c}          (a/(1+exp(-k(t-t0)))+c)
        //   FRACTIONAL : {c0, c1, c2, p1, p2, lo}
        //   FLAT       : {mean}
        std::vector<double> params;
    };

    // =========================================================================
    // BIC
    // =========================================================================

    // Gaussian-error BIC: n * ln(RSS/n) + k * ln(n). Lower is better.
    inline double bic(double rss, int n, int k) {
        if (n <= 0 || rss <= 0.0) return std::numeric_limits<double>::infinity();
        return n * std::log(rss / n) + k * std::log(static_cast<double>(n));
    }

    // =========================================================================
    // Helpers: small linear least-squares solver (normal equations, pivoted
    // Gaussian elimination). Same pattern as feature_marks.cpp's polyfit.
    // =========================================================================

    namespace detail {

        // Solve A*x = b in place (A is m x m row-major, b is m). Returns false
        // if singular.
        inline bool solve_normal(std::vector<double>& A, std::vector<double>& b, int m) {
            for (int i = 0; i < m; ++i) {
                int piv = i;
                for (int r = i + 1; r < m; ++r)
                    if (std::abs(A[r * m + i]) > std::abs(A[piv * m + i])) piv = r;
                if (piv != i) {
                    for (int j = 0; j < m; ++j) std::swap(A[i * m + j], A[piv * m + j]);
                    std::swap(b[i], b[piv]);
                }
                const double diag = A[i * m + i];
                if (std::abs(diag) < 1e-15) return false;
                for (int r = 0; r < m; ++r) {
                    if (r == i) continue;
                    const double f = A[r * m + i] / diag;
                    for (int j = 0; j < m; ++j) A[r * m + j] -= f * A[i * m + j];
                    b[r] -= f * b[i];
                }
            }
            for (int i = 0; i < m; ++i) {
                const double diag = A[i * m + i];
                if (std::abs(diag) > 1e-15) b[i] /= diag; else b[i] = 0.0;
            }
            return true;
        }

        // Linear least-squares fit of basis functions to y[lo..hi]. basis(t,j)
        // returns the j-th basis value at sample index t. Returns coefficients
        // and RSS via out params; returns false if degenerate.
        inline bool linear_ls(const std::vector<double>& y, int lo, int hi, int nbasis,
            std::function<double(double, int)> basis,
            std::vector<double>& coeffs, double& rss)
        {
            const int n = hi - lo + 1;
            const int m = nbasis;
            std::vector<double> A(m * m, 0.0), b(m, 0.0);
            for (int k = 0; k < n; ++k) {
                const double t = static_cast<double>(lo + k);
                const double yv = y[lo + k];
                if (std::isnan(yv)) continue;
                std::vector<double> phi(m);
                for (int j = 0; j < m; ++j) phi[j] = basis(t, j);
                for (int i = 0; i < m; ++i) {
                    b[i] += phi[i] * yv;
                    for (int j = 0; j < m; ++j) A[i * m + j] += phi[i] * phi[j];
                }
            }
            coeffs.resize(m, 0.0);
            if (!solve_normal(A, b, m)) return false;
            coeffs = b;

            rss = 0.0;
            for (int k = 0; k < n; ++k) {
                const double t = static_cast<double>(lo + k);
                const double yv = y[lo + k];
                if (std::isnan(yv)) continue;
                double yhat = 0.0;
                for (int j = 0; j < m; ++j) yhat += coeffs[j] * basis(t, j);
                const double e = yv - yhat;
                rss += e * e;
            }
            return true;
        }

    } // namespace detail

    // =========================================================================
    // Model 1: Piecewise linear
    // =========================================================================
    // Two line segments joined at a breakpoint b. Grid-search b over the
    // interior, least-squares fit each side, keep the b with lowest RSS.
    // Parameters: m1, c1, m2, c2, breakpoint = 5.

    inline FitResult fitPiecewiseLinear(const std::vector<double>& y, int lo, int hi) {
        FitResult best;
        best.nparams = 5;
        best.type = FitType::LINEAR;

        auto lineFit = [&](int a, int b, double& m, double& c) -> double {
            int n = 0;
            double sx = 0, sy = 0, sxx = 0, sxy = 0;
            for (int i = a; i <= b; ++i) {
                if (std::isnan(y[i])) continue;
                const double x = static_cast<double>(i);
                sx += x; sy += y[i]; sxx += x * x; sxy += x * y[i]; ++n;
            }
            if (n < 2) { m = 0.0; c = (n == 1) ? sy : 0.0; return 0.0; }
            const double den = n * sxx - sx * sx;
            m = (std::abs(den) > 1e-15) ? (n * sxy - sx * sy) / den : 0.0;
            c = (sy - m * sx) / n;
            double r = 0.0;
            for (int i = a; i <= b; ++i) {
                if (std::isnan(y[i])) continue;
                const double e = y[i] - (m * i + c);
                r += e * e;
            }
            return r;
            };

        double bestRss = std::numeric_limits<double>::infinity();
        int bestB = (lo + hi) / 2;
        double bm1 = 0, bc1 = 0, bm2 = 0, bc2 = 0;

        // Need at least 2 samples on each side of the breakpoint.
        for (int b = lo + 2; b <= hi - 2; ++b) {
            double mm1, cc1, mm2, cc2;
            const double r = lineFit(lo, b, mm1, cc1) + lineFit(b, hi, mm2, cc2);
            if (r < bestRss) {
                bestRss = r; bestB = b;
                bm1 = mm1; bc1 = cc1; bm2 = mm2; bc2 = cc2;
            }
        }

        best.rss = bestRss;
        best.params = { bm1, bc1, bm2, bc2, static_cast<double>(bestB) };
        const double fm1 = bm1, fc1 = bc1, fm2 = bm2, fc2 = bc2;
        const int fb = bestB;
        best.f = [=](double t) -> double {
            return (t <= fb) ? fm1 * t + fc1 : fm2 * t + fc2;
            };
        return best;
    }

    // =========================================================================
    // Model 2: Sigmoid
    // =========================================================================
    // f(t) = a / (1 + exp(-k * (t - t0))) + c
    // 4 parameters: a, k, t0, c.
    //
    // Fitted with damped Gauss-Newton (Levenberg-Marquardt style).
    // Initialized from the piecewise-linear result: t0 = breakpoint,
    // a = amplitude span across the window, c = baseline level,
    // k = 4 / window_width (transition spans ~the window).

    inline FitResult fitSigmoid(const std::vector<double>& y, int lo, int hi,
        const FitResult& pwInit)
    {
        FitResult result;
        result.nparams = 4;
        result.type = FitType::SIGMOID;
        (void)pwInit;   // init is now data-driven; kept in the signature for callers

        const int n = hi - lo + 1;
        if (n < 5) { result.rss = std::numeric_limits<double>::infinity(); return result; }

        // Collect non-NaN samples.
        std::vector<int> idx;
        idx.reserve(n);
        for (int i = lo; i <= hi; ++i)
            if (!std::isnan(y[i])) idx.push_back(i);
        const int npts = static_cast<int>(idx.size());
        if (npts < 5) { result.rss = std::numeric_limits<double>::infinity(); return result; }

        // Data-driven initialization. Seeding the amplitude from the window
        // ENDPOINTS (pwInit.f(lo/hi)) collapsed to a flat line whenever the
        // transition sat mid-window, because Gauss-Newton then had ~no gradient.
        // Instead: baseline/plateau from robust edge means, t0 at the steepest
        // DATA slope, and k from that actual slope.
        const int edge = std::max(1, npts / 5);
        double yFirst = 0.0, yLast = 0.0;
        for (int i = 0; i < edge; ++i) {
            yFirst += y[idx[i]];
            yLast += y[idx[npts - 1 - i]];
        }
        yFirst /= edge; yLast /= edge;
        double a = yLast - yFirst;     // signed transition amplitude
        double c = yFirst;             // baseline
        double t0 = (lo + hi) / 2.0;
        double maxSlope = 0.0;
        for (size_t j = 1; j < idx.size(); ++j) {
            const int dx = std::max(1, idx[j] - idx[j - 1]);
            const double s = (y[idx[j]] - y[idx[j - 1]]) / dx;
            if (std::abs(s) > std::abs(maxSlope)) {
                maxSlope = s;
                t0 = 0.5 * (idx[j] + idx[j - 1]);
            }
        }
        // sigmoid'(t0) = a*k/4  =>  k = 4*slope/a. Fall back to a gentle slope
        // (correctly signed) if the amplitude is ~0.
        double k = (std::abs(a) > 1e-9)
            ? 4.0 * maxSlope / a
            : ((maxSlope >= 0 ? 1.0 : -1.0) * 4.0 / std::max(1, hi - lo));
        if (std::abs(k) < 1e-6) k = 4.0 / std::max(1, hi - lo);

        // Gauss-Newton with Levenberg damping.
        double lambda = 1e-3;
        double prevRss = std::numeric_limits<double>::infinity();

        auto sigmoid = [](double a, double k, double t0, double c, double t) -> double {
            const double z = -k * (t - t0);
            const double ez = (z > 300.0) ? 1e130 : (z < -300.0) ? 0.0 : std::exp(z);
            return a / (1.0 + ez) + c;
            };

        // Compute RSS for current params.
        auto computeRss = [&](double a_, double k_, double t0_, double c_) -> double {
            double r = 0.0;
            for (int i : idx) {
                const double e = y[i] - sigmoid(a_, k_, t0_, c_, static_cast<double>(i));
                r += e * e;
            }
            return r;
            };

        prevRss = computeRss(a, k, t0, c);

        for (int iter = 0; iter < SIGMOID_MAX_ITER; ++iter) {
            // Build J^T J and J^T r (4x4 normal equations).
            double JtJ[16] = {}, Jtr[4] = {};
            for (int i : idx) {
                const double t = static_cast<double>(i);
                const double z = -k * (t - t0);
                const double ez = (z > 300.0) ? 1e130 : (z < -300.0) ? 0.0 : std::exp(z);
                const double denom = 1.0 + ez;
                const double sig_val = a / denom + c;
                const double residual = y[i] - sig_val;

                // Partial derivatives: df/da, df/dk, df/dt0, df/dc
                const double dfda = 1.0 / denom;
                const double common = a * ez / (denom * denom);
                const double dfdk = common * (t - t0);
                const double dfdt0 = -common * k;
                const double dfdc = 1.0;

                const double J[4] = { dfda, dfdk, dfdt0, dfdc };
                for (int r = 0; r < 4; ++r) {
                    Jtr[r] += J[r] * residual;
                    for (int cc = 0; cc < 4; ++cc)
                        JtJ[r * 4 + cc] += J[r] * J[cc];
                }
            }

            // Levenberg damping: add lambda to diagonal.
            for (int i = 0; i < 4; ++i)
                JtJ[i * 4 + i] *= (1.0 + lambda);

            // Solve 4x4 system for parameter update delta.
            std::vector<double> A(JtJ, JtJ + 16);
            std::vector<double> b(Jtr, Jtr + 4);
            if (!detail::solve_normal(A, b, 4)) break;

            const double na = a + b[0], nk = k + b[1], nt0 = t0 + b[2], nc = c + b[3];
            const double newRss = computeRss(na, nk, nt0, nc);

            if (newRss < prevRss) {
                a = na; k = nk; t0 = nt0; c = nc;
                lambda *= 0.5;
                if (std::abs(prevRss - newRss) < SIGMOID_CONV_TOL * prevRss) {
                    prevRss = newRss;
                    break;
                }
                prevRss = newRss;
            }
            else {
                lambda *= 4.0;
                if (lambda > 1e10) break;   // stuck
            }
        }

        result.rss = prevRss;
        result.params = { a, k, t0, c };
        const double fa = a, fk = k, ft0 = t0, fc = c;
        result.f = [=](double t) -> double {
            return sigmoid(fa, fk, ft0, fc, t);
            };
        return result;
    }

    // =========================================================================
    // Model 3: Fractional polynomial (escalation only)
    // =========================================================================
    // Basis {1, t^p1, t^p2} with powers from the Royston-Altman set
    // {-2, -1, -0.5, 0(=ln), 0.5, 1, 2, 3}. Fit by linear least-squares
    // for each power pair, keep the best by RSS. 3 parameters per pair.

    inline FitResult fitFractionalPolynomial(const std::vector<double>& y, int lo, int hi, double poorRssThreshold = -1.0) {
        static const double powers[] = { -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0 };
        static const int npow = 8;

        FitResult best;
        best.nparams = 3;
        best.type = FitType::FRACTIONAL;

        const int n = hi - lo + 1;
        if (n < 4) return best;

        // Shift t so that t >= 1 (avoids log(0), negative powers of 0).
        // t_shifted = (sample_index - lo) + 1, so t_shifted in [1, n].
        auto fp_basis = [&](double p, double t_shifted) -> double {
            if (p == 0.0) return std::log(t_shifted);
            return std::pow(t_shifted, p);
            };

        // Precompute each of the 8 powers' basis vector over the window ONCE
        // (change #3): fp_basis is a pow/log, and every power reappears in
        // many pairs, so evaluating it per-pair repeats the same 8 vectors
        // ~36x. basisByPow[pi][k] is fp_basis(powers[pi], k+1) for sample
        // lo+k. The linear_ls basis lambda then just indexes into these.
        std::vector<std::vector<double>> basisByPow(npow, std::vector<double>(n));
        for (int pi = 0; pi < npow; ++pi)
            for (int k = 0; k < n; ++k)
                basisByPow[pi][k] = fp_basis(powers[pi], static_cast<double>(k) + 1.0);

        for (int pi = 0; pi < npow; ++pi) {
            for (int pj = pi; pj < npow; ++pj) {
                const double p1 = powers[pi], p2 = powers[pj];

                // t is a sample index in [lo, hi]; k = t - lo indexes the
                // precomputed vectors. (Same ts = (t-lo)+1 convention as
                // fp_basis was called with above.)
                auto basis = [&](double t, int j) -> double {
                    if (j == 0) return 1.0;
                    const int k = static_cast<int>(t) - lo;
                    return (j == 1) ? basisByPow[pi][k] : basisByPow[pj][k];
                    };

                std::vector<double> coeffs;
                double rss;
                if (!detail::linear_ls(y, lo, hi, 3, basis, coeffs, rss)) continue;

                if (rss < best.rss) {
                    best.rss = rss;
                    const double c0 = coeffs[0], c1 = coeffs[1], c2 = coeffs[2];
                    const double pp1 = p1, pp2 = p2;
                    const int flo = lo;
                    best.params = { c0, c1, c2, pp1, pp2, static_cast<double>(flo) };
                    best.f = [=](double t) -> double {
                        const double ts = (t - flo) + 1.0;
                        const double b1 = (pp1 == 0.0) ? std::log(ts) : std::pow(ts, pp1);
                        const double b2 = (pp2 == 0.0) ? std::log(ts) : std::pow(ts, pp2);
                        return c0 + c1 * b1 + c2 * b2;
                        };
                }

                // Early-out (change #1): once a pair fits to the caller's
                // "good enough" threshold, no later pair can change the
                // selectAnchorModel outcome (it only accepts fp if it beats
                // the current winner, and this already clears POOR), so stop
                // grinding through the remaining pairs -- this is exactly the
                // hard-landmark case that was slow.
                if (poorRssThreshold >= 0.0 && best.rss <= poorRssThreshold)
                    return best;
            }
        }
        return best;
        return best;
    }

    // =========================================================================
    // Model selection (tiered, per spec Section 4.2 Step 2)
    // =========================================================================

    // =========================================================================
    // Model 4: natural cubic spline through 3 EVENLY-SPACED knots
    // =========================================================================
    // Knots at the window start / middle / end (x0, x0+h, x0+2h). A natural
    // cubic spline (second derivative 0 at both ends) through those three knots
    // is fully determined by the three knot VALUES v0,v1,v2, which are fit by
    // least squares over the window (cardinal-basis regression: basis j is the
    // spline that is 1 at knot j and 0 at the others). 3 parameters.
    inline FitResult fitCubicSpline(const std::vector<double>& y, int lo, int hi) {
        FitResult r;
        r.type = FitType::CUBIC_SPLINE;
        r.nparams = 3;
        const double x0 = static_cast<double>(lo);
        const double h = (hi - lo) / 2.0;

        auto flat = [&]() {
            double mean = 0.0; int cnt = 0;
            for (int i = lo; i <= hi; ++i) if (!std::isnan(y[i])) { mean += y[i]; ++cnt; }
            if (cnt > 0) mean /= cnt;
            r.type = FitType::FLAT; r.nparams = 1; r.rss = 0.0; r.params = { mean };
            r.f = [=](double) { return mean; };
            return r;
            };
        if (h <= 0.0) return flat();

        // Natural cubic spline value at x for knot values (v0,v1,v2). With the
        // natural end conditions M0 = M2 = 0, the only interior second
        // derivative is M1 = (3/2h^2)(v0 - 2v1 + v2).
        auto spline3 = [x0, h](double x, double v0, double v1, double v2) -> double {
            const double x1 = x0 + h, x2 = x0 + 2.0 * h;
            const double M1 = (3.0 / (2.0 * h * h)) * (v0 - 2.0 * v1 + v2);
            if (x <= x1) {
                const double A = x1 - x, B = x - x0;
                return M1 * B * B * B / (6.0 * h) + (v0 / h) * A + (v1 / h - M1 * h / 6.0) * B;
            }
            const double A = x2 - x, B = x - x1;
            return M1 * A * A * A / (6.0 * h) + (v1 / h - M1 * h / 6.0) * A + (v2 / h) * B;
            };
        auto basis = [spline3](double t, int j) -> double {
            return spline3(t, j == 0 ? 1.0 : 0.0, j == 1 ? 1.0 : 0.0, j == 2 ? 1.0 : 0.0);
            };

        std::vector<double> c; double rss = 0.0;
        if (!detail::linear_ls(y, lo, hi, 3, basis, c, rss) || c.size() < 3)
            return flat();

        r.rss = rss;
        r.params = { c[0], c[1], c[2] };   // the three knot values
        const double v0 = c[0], v1 = c[1], v2 = c[2];
        r.f = [spline3, v0, v1, v2](double t) { return spline3(t, v0, v1, v2); };
        return r;
    }

    inline FitResult selectAnchorModel(const std::vector<double>& y, int lo, int hi,
        FitMode mode = FitMode::Auto) {
        const int n = hi - lo + 1;
        if (n < 5) {
            // Too few samples for any meaningful fit; return a flat line.
            FitResult fallback;
            fallback.rss = 0.0;
            fallback.nparams = 1;
            fallback.type = FitType::FLAT;
            double mean = 0.0; int cnt = 0;
            for (int i = lo; i <= hi; ++i)
                if (!std::isnan(y[i])) { mean += y[i]; ++cnt; }
            if (cnt > 0) mean /= cnt;
            fallback.params = { mean };
            fallback.f = [=](double) { return mean; };
            return fallback;
        }

        // FORCED model: the operator picked one on the radio; return it directly
        // so the anchor is placed from that model (no BIC contest).
        switch (mode) {
        case FitMode::Linear:      return fitPiecewiseLinear(y, lo, hi);
        case FitMode::Sigmoid:     return fitSigmoid(y, lo, hi, fitPiecewiseLinear(y, lo, hi));
        case FitMode::FracPoly:    return fitFractionalPolynomial(y, lo, hi);   // full fit when forced
        case FitMode::CubicSpline: return fitCubicSpline(y, lo, hi);
        case FitMode::Auto:
        default:                   break;
        }

        // AUTO: the cheap models (piecewise, sigmoid, cubic spline) are always
        // fit and compared by BIC. The FRACTIONAL polynomial is expensive (a
        // power-pair search) and this runs on every glyph snapshot of every
        // panel on every alignment switch -- fitting it unconditionally made
        // alignment crawl. So it is only escalated to when the cheap winner
        // fits POORLY (the original gate). The focus panel still fits and draws
        // all four candidates (that path passes candOut and runs only on click).
        FitResult pw = fitPiecewiseLinear(y, lo, hi);
        FitResult sig = fitSigmoid(y, lo, hi, pw);
        FitResult sp = fitCubicSpline(y, lo, hi);

        FitResult best = pw; double bestBic = bic(pw.rss, n, pw.nparams);
        auto consider = [&](const FitResult& r) {
            const double b = bic(r.rss, n, r.nparams);
            if (b < bestBic) { best = r; bestBic = b; }
            };
        consider(sig);
        consider(sp);

        if (best.rss > POOR_FIT_FACTOR * n) {
            FitResult fp = fitFractionalPolynomial(y, lo, hi, POOR_FIT_FACTOR * n);
            consider(fp);
        }
        return best;
    }

    // =========================================================================
    // Anchor placement: fractional crossing
    // =========================================================================

    /**
     * @brief  Find where the fitted curve crosses a target amplitude level.
     *
     *         For onsets (Q-onset, P-onset): f ~ 0.10 (10% of the way from
     *         baseline B toward extremum E).
     *         For offsets (S-end, T-end): f ~ 0.90 (90% of the way from
     *         extremum back toward baseline, i.e. mostly recovered).
     *
     * @param fit   The selected model's FitResult.
     * @param lo    Left edge of the fit window (sample index).
     * @param hi    Right edge of the fit window (sample index).
     * @param B     Baseline level (e.g. PQ isoelectric, ST segment median).
     * @param E     Extremum level (e.g. Q trough, T peak).
     * @param f     Fraction in [0,1] defining the target level L = B + f*(E-B).
     * @return      Sub-sample position of the first crossing, or window
     *              midpoint as fallback.
     */
    inline double anchorAtFraction(const FitResult& fit, int lo, int hi,
        double B, double E, double f)
    {
        const double L = B + f * (E - B);

        // Scan fitted curve for the first sign-change crossing of level L.
        for (int t = lo; t < hi; ++t) {
            const double y0 = fit.f(static_cast<double>(t));
            const double y1 = fit.f(static_cast<double>(t + 1));
            const double d0 = y0 - L;
            const double d1 = y1 - L;

            // Exact hit.
            if (std::abs(d0) < 1e-12) return static_cast<double>(t);

            // Sign change => linear interpolation for sub-sample position.
            if (d0 * d1 < 0.0) {
                const double frac = (L - y0) / (y1 - y0);
                return t + frac;
            }
        }

        // Fallback: no crossing found.
        return (lo + hi) / 2.0;
    }

    // =========================================================================
    // Convenience: anchor for a peak (argmax or argmin of the fitted curve)
    // =========================================================================
    // Not currently used (R/S/T-peak stay on existing detectors), but
    // included per spec Step 3 for completeness.

    inline double anchorAtPeak(const FitResult& fit, int lo, int hi, bool findMin = false) {
        double bestVal = findMin ? std::numeric_limits<double>::infinity()
            : -std::numeric_limits<double>::infinity();
        double bestT = (lo + hi) / 2.0;
        for (int t = lo; t <= hi; ++t) {
            const double v = fit.f(static_cast<double>(t));
            if (findMin ? (v < bestVal) : (v > bestVal)) {
                bestVal = v; bestT = static_cast<double>(t);
            }
        }
        return bestT;
    }

} // namespace anchor_fit