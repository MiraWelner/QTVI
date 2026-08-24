// ppg_pipeline.cpp -- implementations for ppg_pipeline.hpp.

#include "ppg_pipeline.hpp"

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
    std::vector<double> naturalCubicSpline(const std::vector<int>& xs,
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
    void findExtrema(const std::vector<double>& s,
        std::vector<int>& maxs, std::vector<int>& mins) {
        maxs.clear(); mins.clear();
        const int N = static_cast<int>(s.size());
        for (int i = 1; i < N - 1; ++i) {
            if (std::isnan(s[i - 1]) || std::isnan(s[i]) || std::isnan(s[i + 1])) continue;
            if (s[i] >= s[i - 1] && s[i] > s[i + 1]) maxs.push_back(i);
            else if (s[i] <= s[i - 1] && s[i] < s[i + 1]) mins.push_back(i);
        }
    }

    std::vector<double> envelopeThrough(const std::vector<double>& s,
        const std::vector<int>& ext, int i0, int i1) {
        std::vector<int> xs; std::vector<double> ys;
        xs.push_back(i0); ys.push_back(s[i0]);
        for (int e : ext) if (e > i0 && e < i1) { xs.push_back(e); ys.push_back(s[e]); }
        if (xs.back() != i1) { xs.push_back(i1); ys.push_back(s[i1]); }
        return naturalCubicSpline(xs, ys, static_cast<int>(s.size()));
    }

} // anonymous namespace

namespace ppg_pipeline {

    IemEnvelope iemEnvelope(const std::vector<double>& pulse, int maxIter, double sdThresh) {
        const int N = static_cast<int>(pulse.size());
        IemEnvelope R;
        R.envelope.assign(N, kNaN);
        R.imf.assign(N, kNaN);
        if (N < 5) return R;

        int i0 = 0;     while (i0 < N && std::isnan(pulse[i0])) ++i0;
        int i1 = N - 1; while (i1 >= 0 && std::isnan(pulse[i1])) --i1;
        if (i1 - i0 < 4) return R;

        std::vector<double> h = pulse;
        std::vector<int> mx, mn;

        for (int it = 0; it < maxIter; ++it) {
            findExtrema(h, mx, mn);
            if (mx.empty() || mn.empty()) break;   // no oscillation left -> h is the trend

            const std::vector<double> up = envelopeThrough(h, mx, i0, i1);
            const std::vector<double> lo = envelopeThrough(h, mn, i0, i1);

            std::vector<double> hn(N, kNaN);
            double num = 0.0, den = 0.0;
            for (int i = i0; i <= i1; ++i) {
                if (std::isnan(up[i]) || std::isnan(lo[i]) || std::isnan(h[i])) continue;
                const double m = 0.5 * (up[i] + lo[i]);
                hn[i] = h[i] - m;
                num += (h[i] - hn[i]) * (h[i] - hn[i]);
                den += h[i] * h[i];
            }
            R.upper = up; R.lower = lo;
            h.swap(hn);
            R.iterations = it + 1;
            if (den > 0.0 && num / den < sdThresh) break;
        }

        R.imf = h;
        for (int i = 0; i < N; ++i)
            R.envelope[i] = (!std::isnan(pulse[i]) && !std::isnan(h[i])) ? pulse[i] - h[i] : kNaN;
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

        std::vector<double> spline = naturalCubicSpline(xs, ys, N);

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