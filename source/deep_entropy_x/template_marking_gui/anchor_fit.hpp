/**
 * @file   anchor_fit.hpp
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

namespace anchor_fit {

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
    enum class FitType { LINEAR, SIGMOID, FRACTIONAL, FLAT };

    struct FitResult {
        double rss = std::numeric_limits<double>::infinity();
        int    nparams = 0;
        FitType type = FitType::FLAT;
        std::function<double(double)> f;   // evaluator: f(sample_index) -> fitted value
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

        const int n = hi - lo + 1;
        if (n < 5) { result.rss = std::numeric_limits<double>::infinity(); return result; }

        // Collect non-NaN samples.
        std::vector<int> idx;
        idx.reserve(n);
        for (int i = lo; i <= hi; ++i)
            if (!std::isnan(y[i])) idx.push_back(i);
        const int npts = static_cast<int>(idx.size());
        if (npts < 5) { result.rss = std::numeric_limits<double>::infinity(); return result; }

        // Initialize from piecewise-linear.
        const double yLo = pwInit.f(static_cast<double>(lo));
        const double yHi = pwInit.f(static_cast<double>(hi));
        double a = yHi - yLo;
        double c = yLo;
        double t0 = 0.0;
        // Find breakpoint: evaluate pw at each sample, the breakpoint is
        // where the two segments meet -- approximate as midpoint of the window
        // scaled by the piecewise-linear shape.
        {
            double maxSlope = 0.0;
            t0 = (lo + hi) / 2.0;
            for (int i = lo + 1; i < hi; ++i) {
                const double s = std::abs(pwInit.f(i + 0.5) - pwInit.f(i - 0.5));
                if (s > maxSlope) { maxSlope = s; t0 = i; }
            }
        }
        double k = 4.0 / std::max(1, hi - lo);
        if (a < 0) k = -k;   // descending sigmoid

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

    inline FitResult fitFractionalPolynomial(const std::vector<double>& y, int lo, int hi) {
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

        for (int pi = 0; pi < npow; ++pi) {
            for (int pj = pi; pj < npow; ++pj) {
                const double p1 = powers[pi], p2 = powers[pj];

                auto basis = [&](double t, int j) -> double {
                    const double ts = (t - lo) + 1.0;
                    if (j == 0) return 1.0;
                    if (j == 1) return fp_basis(p1, ts);
                    return fp_basis(p2, ts);
                    };

                std::vector<double> coeffs;
                double rss;
                if (!detail::linear_ls(y, lo, hi, 3, basis, coeffs, rss)) continue;

                if (rss < best.rss) {
                    best.rss = rss;
                    const double c0 = coeffs[0], c1 = coeffs[1], c2 = coeffs[2];
                    const double pp1 = p1, pp2 = p2;
                    const int flo = lo;
                    best.f = [=](double t) -> double {
                        const double ts = (t - flo) + 1.0;
                        const double b1 = (pp1 == 0.0) ? std::log(ts) : std::pow(ts, pp1);
                        const double b2 = (pp2 == 0.0) ? std::log(ts) : std::pow(ts, pp2);
                        return c0 + c1 * b1 + c2 * b2;
                        };
                }
            }
        }
        return best;
    }

    // =========================================================================
    // Model selection (tiered, per spec Section 4.2 Step 2)
    // =========================================================================

    inline FitResult selectAnchorModel(const std::vector<double>& y, int lo, int hi) {
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
            fallback.f = [=](double) { return mean; };
            return fallback;
        }

        // Always fit both simple models.
        FitResult pw = fitPiecewiseLinear(y, lo, hi);
        FitResult sig = fitSigmoid(y, lo, hi, pw);

        const double bic_pw = bic(pw.rss, n, pw.nparams);
        const double bic_sig = bic(sig.rss, n, sig.nparams);
        FitResult best = (bic_pw <= bic_sig) ? pw : sig;

        // Escalate to fractional polynomial only if the winner fits poorly.
        // Per spec: POOR = 0.05 * n (tune from clean-template residuals).
        const double POOR = POOR_FIT_FACTOR * n;
        if (best.rss > POOR) {
            FitResult fp = fitFractionalPolynomial(y, lo, hi);
            if (bic(fp.rss, n, fp.nparams) < bic(best.rss, n, best.nparams))
                best = fp;
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