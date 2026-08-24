// ppg_dicrotic.cpp -- implementations for ppg_dicrotic.hpp (E-5).

#include "ppg_dicrotic.hpp"
#include "ppg_derivative.hpp"   // ppg_deriv::buildDerivatives, DerivBank

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

    // -----------------------------------------------------------------
    // Phase-1 stand-ins (replace with the repo's Task N / Task I-3 routines).
    // -----------------------------------------------------------------

    // Natural cubic spline through knots (xs, ys), evaluated on [0, N). Ends are
    // clamped by holding the first/last knot value flat to the array bounds, which
    // prevents the end-swing the spec warns about. (Phase 1 Task N mirrors the end
    // extrema instead; behavior is equivalent for the IEM's purposes.)
    std::vector<double> splineThrough(const std::vector<int>& xs_in,
        const std::vector<double>& ys_in, int N) {
        std::vector<int> xs; std::vector<double> ys;
        if (!xs_in.empty() && xs_in.front() != 0) { xs.push_back(0); ys.push_back(ys_in.front()); }
        for (size_t i = 0; i < xs_in.size(); ++i) {
            if (!xs.empty() && xs.back() == xs_in[i]) continue;
            xs.push_back(xs_in[i]); ys.push_back(ys_in[i]);
        }
        if (!xs.empty() && xs.back() != N - 1) { xs.push_back(N - 1); ys.push_back(ys_in.back()); }

        const int n = static_cast<int>(xs.size());
        std::vector<double> out(N, kNaN);
        if (n == 0) return out;
        if (n == 1) { std::fill(out.begin(), out.end(), ys[0]); return out; }
        if (n == 2) {
            const double m = (ys[1] - ys[0]) / double(xs[1] - xs[0]);
            for (int i = 0; i < N; ++i) out[i] = ys[0] + m * (i - xs[0]);
            return out;
        }
        std::vector<double> h(n - 1), al(n, 0.0), l(n), mu(n), z(n), c(n, 0.0),
            b(n - 1, 0.0), d(n - 1, 0.0);
        for (int i = 0; i < n - 1; ++i) h[i] = double(xs[i + 1] - xs[i]);
        for (int i = 1; i < n - 1; ++i)
            al[i] = 3.0 * ((ys[i + 1] - ys[i]) / h[i] - (ys[i] - ys[i - 1]) / h[i - 1]);
        l[0] = 1.0; mu[0] = 0.0; z[0] = 0.0;
        for (int i = 1; i < n - 1; ++i) {
            l[i] = 2.0 * (xs[i + 1] - xs[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (al[i] - h[i - 1] * z[i - 1]) / l[i];
        }
        l[n - 1] = 1.0; z[n - 1] = 0.0; c[n - 1] = 0.0;
        for (int j = n - 2; j >= 0; --j) {
            c[j] = z[j] - mu[j] * c[j + 1];
            b[j] = (ys[j + 1] - ys[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
            d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
        }
        for (int i = 0; i < N; ++i) {
            int seg = int(std::upper_bound(xs.begin(), xs.end(), i) - xs.begin()) - 1;
            seg = std::clamp(seg, 0, n - 2);
            const double dx = i - xs[seg];
            out[i] = ys[seg] + b[seg] * dx + c[seg] * dx * dx + d[seg] * dx * dx * dx;
        }
        return out;
    }

    // Sub-sample refinement of an extremum by parabolic vertex. (Phase 1 Task I-3
    // fits an asymmetric profile of width sigma; sigma is accepted but unused here.)
    double refineAsymmetric(const std::vector<double>& s, int idx, double /*sigma*/) {
        if (idx <= 0 || idx >= (int)s.size() - 1) return double(idx);
        const double a = s[idx - 1], b = s[idx], c = s[idx + 1];
        const double den = a - 2.0 * b + c;
        if (std::fabs(den) < 1e-300) return double(idx);
        double off = 0.5 * (a - c) / den;
        off = std::clamp(off, -1.0, 1.0);   // a sub-sample shift never exceeds one sample
        return idx + off;
    }

} // anonymous namespace

namespace ppg_dicrotic {

    IemResult iterativeEnvelopeMean(const std::vector<double>& x, int maxIter, double tol) {
        IemResult R;
        std::vector<double> cur = x;
        const int N = static_cast<int>(x.size());
        double x0 = 0.0; for (double v : x) x0 += v * v; x0 = std::sqrt(x0) + 1e-12;
        for (int k = 0; k < maxIter; ++k) {
            std::vector<int> mx, mn;
            for (int i = 1; i < N - 1; ++i) {
                if (cur[i] >= cur[i - 1] && cur[i] > cur[i + 1]) mx.push_back(i);
                if (cur[i] <= cur[i - 1] && cur[i] < cur[i + 1]) mn.push_back(i);
            }
            // Fewer than two extrema of either sign: no further oscillation to
            // separate -- the iterate is already the profile.
            if ((int)mx.size() < 2 || (int)mn.size() < 2) { R.converged = true; break; }
            auto vals = [&](const std::vector<int>& idx) {
                std::vector<double> v; v.reserve(idx.size());
                for (int i : idx) v.push_back(cur[i]);
                return v;
                };
            std::vector<double> up = splineThrough(mx, vals(mx), N);
            std::vector<double> lo = splineThrough(mn, vals(mn), N);
            std::vector<double> mean(N);
            double delta = 0.0;
            for (int i = 0; i < N; ++i) {
                mean[i] = 0.5 * (up[i] + lo[i]);
                delta += (mean[i] - cur[i]) * (mean[i] - cur[i]);
            }
            cur = mean;
            R.iterations = k + 1;
            if (std::sqrt(delta) / x0 < tol) { R.converged = true; break; }
        }
        R.profile = cur;
        R.residual.resize(N);
        for (int i = 0; i < N; ++i) R.residual[i] = x[i] - R.profile[i];
        return R;
    }

    std::vector<double> enhanceDescendingLimb(const std::vector<double>& pulse,
        double fs, int sysPeak, double gain) {
        if (gain <= 0.0) return pulse;
        const int N = static_cast<int>(pulse.size());
        std::vector<double> best(N, 0.0);
        for (double widthMs : {20.0, 30.0, 40.0, 50.0}) {
            const double sig = (widthMs * 1e-3 * fs) / 2.0;   // samples
            const int h = std::max(2, (int)std::ceil(3.0 * sig));
            std::vector<double> k(2 * h + 1);
            for (int j = -h; j <= h; ++j) {
                const double t = j / sig;
                k[j + h] = (t * t - 1.0) * std::exp(-0.5 * t * t) / (sig * sig);  // d2/dt2 Gaussian
            }
            for (int i = 0; i < N; ++i) {
                double s = 0.0;
                for (int j = -h; j <= h; ++j)
                    s += k[j + h] * pulse[std::clamp(i + j, 0, N - 1)];
                s *= sig * sig;                                // scale normalization
                best[i] = std::max(best[i], s);
            }
        }
        const double amp = *std::max_element(pulse.begin(), pulse.end())
            - *std::min_element(pulse.begin(), pulse.end());
        const double bmax = *std::max_element(best.begin(), best.end()) + 1e-12;
        std::vector<double> out = pulse;
        for (int i = std::max(0, sysPeak); i < N; ++i)
            out[i] += gain * amp * best[i] / bmax;
        return out;
    }

    DnResult windkesselDn(const std::vector<double>& pulse, double fs,
        int sysPeak, int winLo, int winHi) {
        const int N = static_cast<int>(pulse.size());
        // 1. Fit tau on the late diastolic decay (last 40%) by log-linear regression.
        const int t0 = (int)(0.60 * N), t1 = N - 1;
        const double pmin = *std::min_element(pulse.begin(), pulse.end());
        double sx = 0, sy = 0, sxx = 0, sxy = 0; int n = 0;
        for (int i = t0; i <= t1; ++i) {
            const double y = pulse[i] - pmin + 1e-6;
            if (y <= 0) continue;
            const double lx = (double)i / fs, ly = std::log(y);
            sx += lx; sy += ly; sxx += lx * lx; sxy += lx * ly; ++n;
        }
        if (n < 8) return {};
        const double slope = (n * sxy - sx * sy) / (n * sxx - sx * sx + 1e-12);
        if (slope >= 0.0) return {};                    // no decay: model does not apply
        const double tau = -1.0 / slope;                // seconds

        // 2. Flow-proportional signal g(t) = dP/dt + P(t)/tau.
        const ppg_deriv::DerivBank D = ppg_deriv::buildDerivatives(pulse, fs);
        std::vector<double> g(N);
        for (int i = 0; i < N; ++i) g[i] = D.d1[i] + (pulse[i] - pmin) / tau;

        // 3. First downward zero-crossing of g inside the window is valve closure.
        for (int i = std::max(winLo, sysPeak + 1); i < std::min(winHi, N); ++i)
            if (g[i - 1] > 0.0 && g[i] <= 0.0) {
                DnResult r; r.index = i; r.tier = DnResult::WINDKESSEL; r.confidence = 0.5;
                return r;
            }
        // 4. No clean crossing: take the minimum of |g| if it is a clear trough.
        int best = -1; double bv = 1e300;
        for (int i = std::max(winLo, sysPeak + 1); i < std::min(winHi, N); ++i)
            if (std::fabs(g[i]) < bv) { bv = std::fabs(g[i]); best = i; }
        if (best > 0) { DnResult r; r.index = best; r.tier = DnResult::WINDKESSEL; r.confidence = 0.25; return r; }
        return {};
    }

    DnResult detectDicroticNotch(const std::vector<double>& pulse, double fs,
        int sysPeak, double rrSeconds, const PpgConfig& cfg) {
        const int N = static_cast<int>(pulse.size());
        const int winLo = sysPeak + (int)(cfg.dnWindowLoMs * 1e-3 * fs);
        const int winHi = std::min(N, sysPeak + (int)(cfg.dnWindowHiRrFrac * rrSeconds * fs));
        if (winHi - winLo < 5) return {};

        // ---- Tier 1: IEM ----
        const std::vector<double> enh =
            enhanceDescendingLimb(pulse, fs, sysPeak, cfg.dnEnhanceGain);
        const IemResult iem = iterativeEnvelopeMean(enh, cfg.iemMaxIter, cfg.iemTol);
        const double amp = *std::max_element(pulse.begin(), pulse.end())
            - *std::min_element(pulse.begin(), pulse.end()) + 1e-12;

        int bestIdx = -1; double bestProm = 0.0;
        const std::vector<double>& r = iem.residual;
        for (int i = winLo + 1; i < winHi - 1; ++i) {
            if (!(r[i] <= r[i - 1] && r[i] < r[i + 1])) continue;   // local min
            double lmax = r[i], rmax = r[i];
            for (int j = i; j >= winLo; --j) lmax = std::max(lmax, r[j]);
            for (int j = i; j < winHi; ++j) rmax = std::max(rmax, r[j]);
            const double prom = std::min(lmax, rmax) - r[i];        // prominence
            if (prom > bestProm) { bestProm = prom; bestIdx = i; }
        }
        if (bestIdx > 0 && bestProm / amp >= cfg.dnMinProminence) {
            DnResult res; res.index = bestIdx; res.tier = DnResult::IEM;
            res.confidence = std::min(1.0, bestProm / (amp * cfg.dnMinProminence * 5.0));
            res.subSample = refineAsymmetric(iem.residual, bestIdx, 10.0);
            return res;
        }

        // ---- Tier 2: Windkessel ----
        DnResult wk = windkesselDn(pulse, fs, sysPeak, winLo, winHi);
        if (wk.index > 0) { wk.subSample = refineAsymmetric(pulse, wk.index, 10.0); return wk; }

        // ---- Tier 3: absent ----
        return {};   // tier == ABSENT
    }

} // namespace ppg_dicrotic