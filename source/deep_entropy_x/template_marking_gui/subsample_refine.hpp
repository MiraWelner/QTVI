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
#include "anchor_fit.hpp"

namespace subsample_refine {

    inline constexpr int kWindowHalfWidth = 7;     // 15-point window = seed +- 7
    inline constexpr double kResidualGuardFrac = 0.10;   // 10% of peak amplitude

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

    // ---------------------------------------------------------------------
    // Symmetric extrema: Gaussian-weighted quadratic, vertex = extremum.
    // ---------------------------------------------------------------------
    inline double symmetricExtremum(const std::vector<double>& signal, int seed, double sigma) {
        const WindowSamples ws = gatherWindow(signal, seed, sigma);
        if (!ws.ok) return static_cast<double>(seed);

        // Weighted least squares for y = a t^2 + b t + c.
        double S0 = 0, S1 = 0, S2 = 0, S3 = 0, S4 = 0, Y0 = 0, Y1 = 0, Y2 = 0;
        for (size_t i = 0; i < ws.t.size(); ++i) {
            const double wt = ws.w[i], ti = ws.t[i], ti2 = ti * ti;
            S0 += wt; S1 += wt * ti; S2 += wt * ti2; S3 += wt * ti2 * ti; S4 += wt * ti2 * ti2;
            Y0 += wt * ws.y[i]; Y1 += wt * ti * ws.y[i]; Y2 += wt * ti2 * ws.y[i];
        }
        std::vector<double> sol;
        const bool ok = solveLinear({ {S4,S3,S2},{S3,S2,S1},{S2,S1,S0} }, { Y2,Y1,Y0 }, sol);
        if (!ok || std::fabs(sol[0]) < 1e-12) return fivePointParabolaExtremum(signal, seed);

        std::vector<double> fitted(ws.t.size());
        for (size_t i = 0; i < ws.t.size(); ++i)
            fitted[i] = sol[0] * ws.t[i] * ws.t[i] + sol[1] * ws.t[i] + sol[2];
        if (weightedRmsResidualFrac(ws, fitted) > kResidualGuardFrac)
            return fivePointParabolaExtremum(signal, seed);

        const double tVertex = -sol[1] / (2.0 * sol[0]);
        return static_cast<double>(seed) + std::clamp(tVertex, (double)-kWindowHalfWidth, (double)kWindowHalfWidth);
    }

    // ---------------------------------------------------------------------
    // Asymmetric extrema: cubic fit on Gaussian-weighted samples, solve
    // dy/dt = 3a t^2 + 2b t + c = 0 analytically; pick the root inside the
    // window closest to t=0 (the seed).
    // ---------------------------------------------------------------------
    inline double asymmetricExtremum(const std::vector<double>& signal, int seed, double sigma) {
        const WindowSamples ws = gatherWindow(signal, seed, sigma);
        if (!ws.ok) return static_cast<double>(seed);

        // Weighted least squares for y = a t^3 + b t^2 + c t + d.
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
        if (!solveLinear(A, { Y[3],Y[2],Y[1],Y[0] }, sol))
            return fivePointParabolaExtremum(signal, seed);

        std::vector<double> fitted(ws.t.size());
        for (size_t i = 0; i < ws.t.size(); ++i) {
            const double ti = ws.t[i];
            fitted[i] = sol[0] * ti * ti * ti + sol[1] * ti * ti + sol[2] * ti + sol[3];
        }
        if (weightedRmsResidualFrac(ws, fitted) > kResidualGuardFrac)
            return fivePointParabolaExtremum(signal, seed);

        // dy/dt = 3a t^2 + 2b t + c = 0
        const double a = sol[0], b = sol[1], c = sol[2];
        if (std::fabs(a) < 1e-12) {
            if (std::fabs(b) < 1e-12) return fivePointParabolaExtremum(signal, seed);
            const double t = -c / (2.0 * b);
            return static_cast<double>(seed) + std::clamp(t, (double)-kWindowHalfWidth, (double)kWindowHalfWidth);
        }
        const double disc = 4.0 * b * b - 12.0 * a * c;
        if (disc < 0.0) return fivePointParabolaExtremum(signal, seed);
        const double sq = std::sqrt(disc);
        const double r1 = (-2.0 * b + sq) / (6.0 * a);
        const double r2 = (-2.0 * b - sq) / (6.0 * a);
        // Pick whichever root is inside the window and closer to t=0.
        const bool r1In = r1 >= -kWindowHalfWidth && r1 <= kWindowHalfWidth;
        const bool r2In = r2 >= -kWindowHalfWidth && r2 <= kWindowHalfWidth;
        double tBest;
        if (r1In && r2In) tBest = (std::fabs(r1) < std::fabs(r2)) ? r1 : r2;
        else if (r1In) tBest = r1;
        else if (r2In) tBest = r2;
        else return fivePointParabolaExtremum(signal, seed);
        return static_cast<double>(seed) + tBest;
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
    // Transition onsets/offsets: locally upsample a 40-sample window from
    // its native rate to 4x via cubic interpolation, then fit-and-select
    // (Section 4.2 machinery, Phase A) on the upsampled window, returning a
    // sub-sample position in the ORIGINAL sample-rate coordinate.
    // ---------------------------------------------------------------------
    inline double transitionAnchor(const std::vector<double>& signal, int seed,
        double fraction, int windowSamples = 40,
        double externalBaseline = std::numeric_limits<double>::quiet_NaN(),
        double boundLo = -1.0, double boundHi = -1.0) {
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

        auto fit = anchor_fit::selectAnchorModel(up, 0, nOut - 1);
        const double anchorUp = anchor_fit::anchorAtFraction(fit, 0, nOut - 1, B, E, fraction);
        (void)seedUp;
        // Back to original-rate coordinate: anchorUp is a position on the
        // 4x grid starting at local[0] == signal[lo].
        return static_cast<double>(lo) + anchorUp / upsampleFactor;
    }


    inline int cubicSplineNotch(const std::vector<double>&, int lo, int, int*) { return lo; }
}  // namespace subsample_refine