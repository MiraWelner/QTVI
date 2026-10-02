// ppg_pipeline.cpp -- implementations for ppg_pipeline.hpp.

#include "ppg_pipeline.hpp"
#include "fiducial_marker_finding/ppg_derivative.hpp"   // the one Savitzky-Golay bank
#include "fiducial_marker_finding/ppg_dicrotic.hpp"    // extremaOfFirstDerivative

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
    constexpr double kPi = 3.14159265358979323846;   // M_PI is not portable (MSVC needs _USE_MATH_DEFINES)

    // -----------------------------------------------------------------
    // Natural cubic spline through knots (xs, ys), xs strictly increasing,
    // sampled at every integer index in [0, N). Tails outside the knot span
    // are left NaN (callers that need full coverage anchor the endpoints).
    // -----------------------------------------------------------------
    std::vector<double> naturalCubicSplineImpl(const std::vector<int>& xs,
        const std::vector<double>& ys, int N) {
        const int n = static_cast<int>(xs.size());
        std::vector<double> out(N, kNaN);
        if (n == 0 || N <= 0) return out;
        if (n == 1) { std::fill(out.begin(), out.end(), ys[0]); return out; }
        if (n == 2) {
            const double m = (ys[1] - ys[0]) / static_cast<double>(xs[1] - xs[0]);
            for (int i = 0; i < N; ++i) out[i] = ys[0] + m * (i - xs[0]);
            for (int i = 0; i < xs.front(); ++i) out[i] = kNaN;
            for (int i = xs.back() + 1; i < N; ++i) out[i] = kNaN;
            return out;
        }
        std::vector<double> h(n - 1), alpha(n, 0.0), l(n), mu(n), z(n),
            c(n, 0.0), b(n - 1, 0.0), dcoef(n - 1, 0.0);
        for (int i = 0; i < n - 1; ++i) h[i] = static_cast<double>(xs[i + 1] - xs[i]);
        for (int i = 1; i < n - 1; ++i)
            alpha[i] = 3.0 * ((ys[i + 1] - ys[i]) / h[i] - (ys[i] - ys[i - 1]) / h[i - 1]);
        l[0] = 1.0; mu[0] = 0.0; z[0] = 0.0;                 // natural BC
        for (int i = 1; i < n - 1; ++i) {
            l[i] = 2.0 * (xs[i + 1] - xs[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }
        l[n - 1] = 1.0; z[n - 1] = 0.0; c[n - 1] = 0.0;      // natural BC
        for (int j = n - 2; j >= 0; --j) {
            c[j] = z[j] - mu[j] * c[j + 1];
            b[j] = (ys[j + 1] - ys[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
            dcoef[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
        }
        for (int i = 0; i < N; ++i) {
            int seg = static_cast<int>(
                std::upper_bound(xs.begin(), xs.end(), i) - xs.begin()) - 1;
            seg = std::clamp(seg, 0, n - 2);
            const double dx = i - xs[seg];
            out[i] = ys[seg] + b[seg] * dx + c[seg] * dx * dx + dcoef[seg] * dx * dx * dx;
        }
        for (int i = 0; i < xs.front(); ++i)    out[i] = kNaN;
        for (int i = xs.back() + 1; i < N; ++i) out[i] = kNaN;
        return out;
    }

    // -----------------------------------------------------------------
    // 2nd-order Butterworth low-pass (bilinear transform, Q = 1/sqrt(2))
    // and a zero-phase forward-backward apply (filtfilt), with odd
    // reflection padding to damp edge transients.
    // -----------------------------------------------------------------
    struct Biquad { double b0, b1, b2, a1, a2; };

    Biquad butter2Lowpass(double fc, double fs) {
        const double K = std::tan(kPi * fc / fs);
        const double Q = 0.70710678118654752440;   // Butterworth
        const double norm = 1.0 / (1.0 + K / Q + K * K);
        Biquad bq;
        bq.b0 = K * K * norm;
        bq.b1 = 2.0 * bq.b0;
        bq.b2 = bq.b0;
        bq.a1 = 2.0 * (K * K - 1.0) * norm;
        bq.a2 = (1.0 - K / Q + K * K) * norm;
        return bq;
    }

    std::vector<double> biquadForward(const std::vector<double>& x, const Biquad& q) {
        const int N = static_cast<int>(x.size());
        std::vector<double> y(N, 0.0);
        double z1 = 0.0, z2 = 0.0;               // direct form II transposed
        for (int i = 0; i < N; ++i) {
            const double in = x[i];
            const double out = q.b0 * in + z1;
            z1 = q.b1 * in - q.a1 * out + z2;
            z2 = q.b2 * in - q.a2 * out;
            y[i] = out;
        }
        return y;
    }

    std::vector<double> filtfilt(const std::vector<double>& x, const Biquad& q, double fs, double fc) {
        const int N = static_cast<int>(x.size());
        if (N < 2) return x;
        // Pad ~ one cutoff period each side, capped to the signal length.
        int pad = static_cast<int>(std::lround(fs / std::max(fc, 1e-6)));
        pad = std::clamp(pad, 3, N - 1);

        std::vector<double> ext(N + 2 * pad);
        for (int k = 0; k < pad; ++k) ext[k] = 2.0 * x[0] - x[pad - k];             // odd reflection
        for (int i = 0; i < N; ++i)   ext[pad + i] = x[i];
        for (int k = 0; k < pad; ++k) ext[pad + N + k] = 2.0 * x[N - 1] - x[N - 2 - k];

        std::vector<double> f = biquadForward(ext, q);
        std::reverse(f.begin(), f.end());
        f = biquadForward(f, q);
        std::reverse(f.begin(), f.end());

        return std::vector<double>(f.begin() + pad, f.begin() + pad + N);
    }

    // -----------------------------------------------------------------
    // Interior local extrema, NaN-skipping (used by the IEM sift).
    // -----------------------------------------------------------------
    // findExtrema REMOVED. It took 3-point extrema of the SIGNAL; E-5.3
    // requires extrema of the FIRST DERIVATIVE, located by sign changes of the
    // SECOND. Those are different and larger knot sets, and the knot set is
    // what the envelopes -- and therefore the residual the notch lives in --
    // are built from. ppg_dicrotic::extremaOfFirstDerivative is now the single
    // definition; see the note there on why the sign of d2 must not be used.

    std::vector<double> envelopeThrough(const std::vector<double>& s,
        const std::vector<int>& ext, int i0, int i1) {
        std::vector<int> xs; std::vector<double> ys;
        xs.push_back(i0); ys.push_back(s[i0]);
        for (int e : ext) if (e > i0 && e < i1) { xs.push_back(e); ys.push_back(s[e]); }
        if (xs.back() != i1) { xs.push_back(i1); ys.push_back(s[i1]); }
        return naturalCubicSplineImpl(xs, ys, static_cast<int>(s.size()));
    }

} // anonymous namespace

namespace ppg_pipeline {

    // E-5.2: Savitzky-Golay coefficient count nc, a WINDOW length, odd, about
    // one half-width (0.1 s) of the shortest cardiac feature. 25 at 256 Hz,
    // 49 at 500 Hz. The authors name non-adaptive SG parameters as their first
    // limitation and require nc be updated for a different sampling rate.
    static int sgNc(double fs) {
        int n = (fs > 0.0) ? static_cast<int>(std::lround(0.098 * fs)) : 25;
        if (n % 2 == 0) ++n;
        return std::max(n, 7);
    }

    // E-5.3's `splineThrough`. One definition, shared by the IEM's envelopes
    // and E-1's DC baseline.
    std::vector<double> naturalCubicSpline(const std::vector<int>& xs,
        const std::vector<double>& ys, int N) {
        return naturalCubicSplineImpl(xs, ys, N);
    }

    IemEnvelope iemEnvelope(const std::vector<double>& pulse, double fs,
        int maxIter, double beta) {
        const int N = static_cast<int>(pulse.size());
        IemEnvelope R;
        R.nonStationary.assign(N, kNaN);
        R.stationary.assign(N, 0.0);
        if (N < 5 || !(fs > 0.0)) return R;

        int i0 = 0;     while (i0 < N && std::isnan(pulse[i0])) ++i0;
        int i1 = N - 1; while (i1 >= 0 && std::isnan(pulse[i1])) --i1;
        if (i1 - i0 < 4) return R;

        std::vector<double> y = pulse;          // current UNSMOOTHED iterate
        double prevEnergy = 0.0;                // R_0 initialised to zero

        for (int it = 0; it < maxIter; ++it) {
            // Envelopes are built on the SMOOTHED signal, and the knots are
            // extrema of its FIRST derivative located by sign changes of the
            // SECOND (E-5.3). ppg_deriv supplies the Savitzky-Golay bank --
            // d0 is the smoothed signal, d1/d2 the derivatives.
            // nc per E-5.2: 25 at 256 Hz, 49 at 500 Hz. ppg_deriv takes a
            // HALF-width and its own default (12) would give nc = 25 at
            // 500 Hz -- half what E-5.2 requires.
            const int h = (sgNc(fs) - 1) / 2;
            const ppg_deriv::DerivBank D = ppg_deriv::buildDerivatives(y, fs, h, h, 4);
            if (D.d1.empty() || D.d2.empty()) { R.converged = true; break; }

            std::vector<int> mx, mn;
            ppg_dicrotic::extremaOfFirstDerivative(D.d1, D.d2, mx, mn);
            if ((int)mx.size() < 2 || (int)mn.size() < 2) { R.converged = true; break; }

            // Knot VALUES are the smoothed signal at the knot positions.
            const std::vector<double> up = envelopeThrough(D.d0, mx, i0, i1);
            const std::vector<double> lo = envelopeThrough(D.d0, mn, i0, i1);

            // ...but subtraction is from the UNSMOOTHED current iterate.
            // Equation (3).
            std::vector<double> r(N, kNaN);
            double energy = 0.0;
            int n = 0;
            for (int i = i0; i <= i1; ++i) {
                if (std::isnan(up[i]) || std::isnan(lo[i]) || std::isnan(y[i])) continue;
                const double m = 0.5 * (up[i] + lo[i]);
                r[i] = y[i] - m;
                R.stationary[i] += m;
                energy += r[i] * r[i];
                ++n;
            }
            if (n == 0) { R.converged = true; break; }
            energy /= n;                        // E{R_i^2(n)}

            R.upper = up; R.lower = lo;
            y.swap(r);                          // iterate on the residual
            R.iterations = it + 1;

            // STC_i = |E{R_{i-1}^2} - E{R_i^2}| < beta
            if (std::fabs(prevEnergy - energy) < beta) { R.converged = true; break; }
            prevEnergy = energy;
        }

        R.nonStationary = y;
        R.ok = true;
        return R;
    }

    std::vector<double> dcEnvelope(const std::vector<double>& ppg,
        const std::vector<int>& troughs, double fs) {
        const int N = static_cast<int>(ppg.size());
        std::vector<double> env(N, kNaN);
        if (N < 3 || fs <= 0.0) return env;

        // Knots at the finite, in-range troughs (sorted, de-duplicated x).
        std::vector<int> xs; std::vector<double> ys;
        std::vector<int> tr = troughs;
        std::sort(tr.begin(), tr.end());
        for (int t : tr) {
            if (t < 0 || t >= N || std::isnan(ppg[t])) continue;
            if (!xs.empty() && t == xs.back()) continue;
            xs.push_back(t); ys.push_back(ppg[t]);
        }
        if (xs.empty()) return env;

        // Anchor the record ends by holding the outermost trough values flat,
        // so the baseline is defined across the whole recording.
        if (xs.front() != 0) { xs.insert(xs.begin(), 0);   ys.insert(ys.begin(), ys.front()); }
        if (xs.back() != N - 1) { xs.push_back(N - 1);        ys.push_back(ys.back()); }

        std::vector<double> spline = naturalCubicSplineImpl(xs, ys, N);

        // 2nd-order Butterworth LPF at 0.1 Hz, zero-phase.
        const double fc = 0.1;
        if (fs <= 2.0 * fc) return spline;   // cutoff not resolvable at this rate
        const Biquad q = butter2Lowpass(fc, fs);
        return filtfilt(spline, q, fs, fc);
    }

    std::vector<double> perfusionIndex(const std::vector<double>& sys,
        const std::vector<double>& dc) {
        const size_t n = std::min(sys.size(), dc.size());
        std::vector<double> pi(sys.size(), kNaN);

        // 1st-percentile floor over the finite DC values.
        std::vector<double> s;
        s.reserve(n);
        for (size_t t = 0; t < n; ++t) if (std::isfinite(dc[t])) s.push_back(dc[t]);
        if (s.empty()) return pi;
        std::sort(s.begin(), s.end());
        const double floorVal = s[static_cast<size_t>(0.01 * s.size())];

        for (size_t t = 0; t < n; ++t) {
            if (!std::isfinite(sys[t]) || !std::isfinite(dc[t])) continue;
            const double safe = std::max(dc[t], floorVal);
            if (!(safe > 0.0)) continue;                      // non-positive baseline: drop
            const double val = (sys[t] - dc[t]) / safe * 100.0;
            pi[t] = (val < 0.1) ? kNaN : val;                 // exclude sub-0.1% beats
        }
        return pi;
    }

} // namespace ppg_pipeline