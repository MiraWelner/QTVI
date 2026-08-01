// ============================================================================
// File: FilterUtils.hpp
// Digital filtering utilities
// Includes: filter, filtfilt, butter (lowpass/highpass/bandpass), conv, medfilt1
// ============================================================================
#pragma once

#include <vector>
#include <string>
#include <complex>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

// FilterUtils uses the project-wide `vector` / `string` aliases (no std::
// prefix). Historically these were provided by a project-level force-include
// hitting every translation unit; a couple of TUs miss it (e.g. the caller
// that surfaced C7568 here), so pull them in directly to make this header
// stand on its own without breaking the callers that already have them.
using std::vector;
using std::string;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace filter_detail {

    using cd = std::complex<double>;

    // Polynomial from complex roots � matches MATLAB poly(roots)
    inline vector<double> poly_from_roots(const vector<cd>& roots) {
        int n = (int)roots.size();
        vector<cd> c(n + 1, cd(0.0, 0.0));
        c[0] = cd(1.0, 0.0);
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j >= 1; j--)
                c[j] -= roots[i] * c[j - 1];
        vector<double> result(n + 1);
        for (int i = 0; i <= n; i++)
            result[i] = c[i].real();
        return result;
    }

    // Filter with initial conditions (Direct Form II Transposed)
    inline vector<double> filter_ic(const vector<double>& b, const vector<double>& a,
        const vector<double>& x, vector<double> zi) {
        if (x.empty() || a.empty()) return {};
        size_t n = x.size();
        size_t nb = b.size();
        size_t na = a.size();
        size_t n_order = std::max(nb, na);
        zi.resize(n_order, 0.0);

        vector<double> y(n);
        double a0 = a[0];

        for (size_t i = 0; i < n; ++i) {
            y[i] = b[0] * x[i] / a0 + zi[0];
            for (size_t j = 1; j < n_order; ++j) {
                double bj = (j < nb) ? b[j] : 0.0;
                double aj = (j < na) ? a[j] : 0.0;
                if (j < n_order - 1)
                    zi[j - 1] = zi[j] + (bj * x[i] - aj * y[i]) / a0;
                else
                    zi[j - 1] = (bj * x[i] - aj * y[i]) / a0;
            }
        }
        return y;
    }

    // Compute initial conditions for filtfilt
    inline vector<double> compute_filtfilt_ic(const vector<double>& b, const vector<double>& a) {
        int n_order = (int)std::max(b.size(), a.size());
        int m = n_order - 1;
        if (m <= 0) return {};

        double a0 = a[0];
        vector<double> bn(n_order, 0.0), an(n_order, 0.0);
        for (int i = 0; i < (int)b.size() && i < n_order; i++) bn[i] = b[i] / a0;
        for (int i = 0; i < (int)a.size() && i < n_order; i++) an[i] = a[i] / a0;

        vector<vector<double>> mat(m, vector<double>(m, 0.0));
        vector<double> rhs(m, 0.0);

        mat[0][0] = 1.0 + an[1];
        for (int i = 1; i < m; i++)
            mat[i][0] = (i + 1 < n_order) ? an[i + 1] : 0.0;
        for (int i = 1; i < m; i++)
            mat[i][i] = 1.0;
        for (int i = 0; i < m - 1; i++)
            mat[i][i + 1] = -1.0;

        for (int i = 0; i < m; i++) {
            double bi = (i + 1 < (int)b.size()) ? b[i + 1] / a0 : 0.0;
            double ai = (i + 1 < (int)a.size()) ? a[i + 1] / a0 : 0.0;
            rhs[i] = bi - bn[0] * ai;
        }

        // Gaussian elimination with partial pivoting
        for (int col = 0; col < m; col++) {
            int pivot = col;
            double maxval = std::fabs(mat[col][col]);
            for (int row = col + 1; row < m; row++) {
                if (std::fabs(mat[row][col]) > maxval) {
                    maxval = std::fabs(mat[row][col]);
                    pivot = row;
                }
            }
            if (pivot != col) { std::swap(mat[col], mat[pivot]); std::swap(rhs[col], rhs[pivot]); }
            if (std::fabs(mat[col][col]) < 1e-15) continue;
            for (int row = col + 1; row < m; row++) {
                double factor = mat[row][col] / mat[col][col];
                for (int j = col; j < m; j++) mat[row][j] -= factor * mat[col][j];
                rhs[row] -= factor * rhs[col];
            }
        }

        vector<double> zi(m, 0.0);
        for (int i = m - 1; i >= 0; i--) {
            double s = rhs[i];
            for (int j = i + 1; j < m; j++) s -= mat[i][j] * zi[j];
            if (std::fabs(mat[i][i]) > 1e-15) zi[i] = s / mat[i][i];
        }
        return zi;
    }

} // namespace filter_detail

// ============================================================================
// Direct Form II Transposed IIR filter � matches y = filter(b, a, x)
// ============================================================================
inline vector<double> filter(const vector<double>& b, const vector<double>& a, const vector<double>& x) {
    if (x.empty() || a.empty()) return {};
    size_t n = x.size();
    size_t nb = b.size();
    size_t na = a.size();
    size_t n_order = std::max(nb, na);

    // Fast path: pure FIR with uniform b. Result = (b[0]) * causal sliding sum
    // of x over the last nb samples. This is the case for rpeakdetect step 5
    // (filter(d_kernel, {1.0}, sqr) with d_kernel = vector(W, 1.0)). The general
    // direct-form-II-transposed loop above would do O(N*W) work; the sliding
    // sum is O(N). Results match the original within FP reorder noise.
    if (na == 1 && nb >= 2 && a[0] != 0.0) {
        bool uniform = true;
        for (size_t j = 1; j < nb; ++j) {
            if (b[j] != b[0]) { uniform = false; break; }
        }
        if (uniform) {
            const double c = b[0] / a[0];
            vector<double> y(n);
            double s = 0.0;
            for (size_t i = 0; i < n; ++i) {
                s += x[i];
                if (i >= nb) s -= x[i - nb];
                y[i] = c * s;
            }
            return y;
        }
    }

    vector<double> y(n);
    vector<double> z(n_order, 0.0);
    double a0 = a[0];

    for (size_t i = 0; i < n; ++i) {
        y[i] = b[0] * x[i] / a0 + z[0];
        for (size_t j = 1; j < n_order; ++j) {
            double bj = (j < nb) ? b[j] : 0.0;
            double aj = (j < na) ? a[j] : 0.0;
            if (j < n_order - 1)
                z[j - 1] = z[j] + (bj * x[i] - aj * y[i]) / a0;
            else
                z[j - 1] = (bj * x[i] - aj * y[i]) / a0;
        }
    }
    return y;
}

// ============================================================================
// Zero-phase filtering (mirror padding + initial conditions)
// ============================================================================
inline vector<double> filtfilt(const vector<double>& b, const vector<double>& a, const vector<double>& x) {
    if (x.empty() || a.empty() || b.empty()) return x;

    int nfact = 3 * (std::max((int)b.size(), (int)a.size()) - 1);
    if ((int)x.size() <= nfact)
        return filter(b, a, x);

    auto zi_unit = filter_detail::compute_filtfilt_ic(b, a);
    int m = (int)zi_unit.size();

    vector<double> padded;
    padded.reserve(x.size() + 2 * nfact);
    for (int i = nfact; i > 0; --i)
        padded.push_back(2.0 * x[0] - x[i]);
    padded.insert(padded.end(), x.begin(), x.end());
    for (int i = 1; i <= nfact; ++i)
        padded.push_back(2.0 * x.back() - x[x.size() - 1 - i]);

    vector<double> zi_fwd(m);
    for (int i = 0; i < m; i++) zi_fwd[i] = zi_unit[i] * padded[0];
    vector<double> y = filter_detail::filter_ic(b, a, padded, zi_fwd);

    std::reverse(y.begin(), y.end());

    vector<double> zi_bwd(m);
    for (int i = 0; i < m; i++) zi_bwd[i] = zi_unit[i] * y[0];
    y = filter_detail::filter_ic(b, a, y, zi_bwd);

    std::reverse(y.begin(), y.end());
    if ((int)y.size() < 2 * nfact) return x;
    return vector<double>(y.begin() + nfact, y.end() - nfact);
}

// ============================================================================
// Butterworth lowpass/highpass design (arbitrary order N, ZPK pipeline).
//
// Replaces an earlier direct-3rd-order formula. That version silently ignored
// N (always producing 3rd-order) AND had a sign bug in the HP denominator
// that put its poles outside the unit circle at high frequencies -- the
// "highpass" was really an unstable filter that exponentially amplified
// content well above cutoff. This ZPK version follows the same design flow
// as the bandpass overload below:
//
//   1) Place N analog Butterworth prototype poles on the left-half unit
//      circle (zeros are all at infinity).
//   2) Prewarp the digital cutoff frequency for the bilinear transform.
//   3) Apply the s-domain LP->LP or LP->HP transformation (scale by cutoff,
//      or invert-and-scale).
//   4) Bilinear-transform to z-domain: (2*fs + s) / (2*fs - s).
//   5) Expand pole/zero polynomials, normalize DC gain (LP) or Nyquist
//      gain (HP) to unity.
//
// Wn is the digital cutoff normalized to Nyquist (0 < Wn < 1), same
// convention as MATLAB's butter().
// ============================================================================
inline void butter(int N, double Wn, const string& type, vector<double>& b, vector<double>& a) {
    using cd = std::complex<double>;
    if (N < 1) N = 1;

    const double fs = 2.0;
    const double u = 2.0 * fs * std::tan(M_PI * Wn / fs);   // prewarped cutoff

    // (1) Analog Butterworth prototype poles: unit-circle, left half plane.
    vector<cd> proto_poles;
    proto_poles.reserve(N);
    for (int k = 1; k <= N; ++k) {
        const double angle = M_PI * (2.0 * k + N - 1.0) / (2.0 * N);
        proto_poles.push_back(cd(std::cos(angle), std::sin(angle)));
    }

    // (2)+(3) LP->LP: multiply each proto pole by u.  LP->HP: divide u by
    // each proto pole (and place N zeros at s=0).
    vector<cd> s_poles, s_zeros;
    s_poles.reserve(N);
    if (type == "low") {
        for (const auto& p : proto_poles) s_poles.push_back(p * u);
        // No finite zeros; bilinear will add N zeros at z=-1.
    }
    else if (type == "high") {
        for (const auto& p : proto_poles) s_poles.push_back(u / p);
        s_zeros.assign(N, cd(0.0, 0.0));
    }
    else {
        // Unknown type: leave b/a empty so callers can detect the mistake.
        b.clear(); a.clear();
        return;
    }

    // (4) Bilinear s -> z:  z = (2*fs + s) / (2*fs - s).
    vector<cd> z_poles, z_zeros;
    z_poles.reserve(s_poles.size());
    z_zeros.reserve(std::max<size_t>(s_zeros.size(), s_poles.size()));
    for (const auto& s : s_poles)
        z_poles.push_back((2.0 * fs + s) / (2.0 * fs - s));
    for (const auto& s : s_zeros)
        z_zeros.push_back((2.0 * fs + s) / (2.0 * fs - s));
    // Missing finite zeros get placed at z=-1 (the LP case: N zeros there).
    while (z_zeros.size() < z_poles.size())
        z_zeros.push_back(cd(-1.0, 0.0));

    // (5) Expand to polynomial coefficients, then normalize DC (LP) or
    // Nyquist (HP) gain to unity so the passband peaks at 0 dB.
    a = filter_detail::poly_from_roots(z_poles);
    b = filter_detail::poly_from_roots(z_zeros);

    // Normalize a[0] to 1 first (poly_from_roots produces monic leading term
    // already, but keep this defensive against any FP roundoff).
    const double a0 = a[0];
    for (auto& v : a) v /= a0;
    for (auto& v : b) v /= a0;

    // Passband gain at z = 1 (LP) or z = -1 (HP).
    auto polyEval = [](const vector<double>& p, double z) {
        double zp = 1.0, s = 0.0;
        for (size_t i = 0; i < p.size(); ++i) { s += p[p.size() - 1 - i] * zp; zp *= z; }
        return s;
        };
    const double z_eval = (type == "low") ? 1.0 : -1.0;
    const double num = polyEval(b, z_eval);
    const double den = polyEval(a, z_eval);
    const double gain = num / den;
    if (std::fabs(gain) > 1e-15)
        for (auto& v : b) v /= gain;
}

// ============================================================================
// Butterworth bandpass design (ZPK pipeline)
// ============================================================================
inline void butter(int N, const vector<double>& Wn, vector<double>& b, vector<double>& a) {
    if (Wn.size() < 2) return;
    using cd = std::complex<double>;

    const double fs = 2.0;
    double u1 = 2.0 * fs * std::tan(M_PI * Wn[0] / fs);
    double u2 = 2.0 * fs * std::tan(M_PI * Wn[1] / fs);
    double BW = u2 - u1;
    double W0 = std::sqrt(u1 * u2);

    vector<cd> proto_poles;
    for (int k = 1; k <= N; k++) {
        double angle = M_PI * (2.0 * k + N - 1.0) / (2.0 * N);
        proto_poles.push_back(cd(std::cos(angle), std::sin(angle)));
    }

    vector<cd> bp_poles;
    for (const auto& p : proto_poles) {
        cd pBW2 = p * BW / 2.0;
        cd sq = std::sqrt(pBW2 * pBW2 - W0 * W0);
        bp_poles.push_back(pBW2 + sq);
        bp_poles.push_back(pBW2 - sq);
    }

    vector<cd> bp_zeros(N, cd(0.0, 0.0));
    double bp_gain = std::pow(BW, N);

    cd num_prod(1.0, 0.0);
    for (const auto& z : bp_zeros) num_prod *= (2.0 * fs - z);
    cd den_prod(1.0, 0.0);
    for (const auto& p : bp_poles) den_prod *= (2.0 * fs - p);
    double d_gain = bp_gain * (num_prod / den_prod).real();

    vector<cd> d_poles, d_zeros;
    for (const auto& p : bp_poles)
        d_poles.push_back((1.0 + p / (2.0 * fs)) / (1.0 - p / (2.0 * fs)));
    for (const auto& z : bp_zeros)
        d_zeros.push_back((1.0 + z / (2.0 * fs)) / (1.0 - z / (2.0 * fs)));
    while (d_zeros.size() < d_poles.size())
        d_zeros.push_back(cd(-1.0, 0.0));

    a = filter_detail::poly_from_roots(d_poles);
    b = filter_detail::poly_from_roots(d_zeros);

    double a0 = a[0];
    for (auto& v : b) v *= d_gain / a0;
    for (auto& v : a) v /= a0;
}

// ============================================================================
// Full convolution: output length = a.size() + b.size() - 1
// Includes a fast path for uniform kernels (all elements equal): a uniform
// kernel of width W reduces convolution to a scaled sliding sum, O(N) instead
// of O(N*W). Results agree with the naive form to within floating-point
// reorder noise (~1e-15 relative).
// ============================================================================
inline vector<double> conv(const vector<double>& a, const vector<double>& b) {
    if (a.empty() || b.empty()) return {};

    // Fast path: uniform kernel detection. Cheap to check, huge payoff
    // for the moving-average kernels in pan_tompkin and rpeakdetect.
    auto try_uniform = [](const vector<double>& sig, const vector<double>& kernel,
        vector<double>& out) -> bool {
            if (kernel.size() < 2) return false;
            const double c = kernel[0];
            for (size_t i = 1; i < kernel.size(); ++i) {
                if (kernel[i] != c) return false;
            }
            const size_t N = sig.size();
            const size_t W = kernel.size();
            const size_t outN = N + W - 1;
            out.assign(outN, 0.0);
            // Sliding sum: at output index k, accumulate sig[max(0,k-W+1)..min(N-1,k)]
            double s = 0.0;
            for (size_t k = 0; k < outN; ++k) {
                if (k < N) s += sig[k];
                if (k >= W) s -= sig[k - W];
                out[k] = c * s;
            }
            return true;
        };

    vector<double> res;
    if (try_uniform(a, b, res)) return res;
    if (try_uniform(b, a, res)) return res;

    // General path.
    res.assign(a.size() + b.size() - 1, 0.0);
    for (size_t i = 0; i < a.size(); ++i)
        for (size_t j = 0; j < b.size(); ++j)
            res[i + j] += a[i] * b[j];
    return res;
}

// ============================================================================
// 1D median filter
// ============================================================================
inline std::vector<double> medfilt1(const std::vector<double>& x, int n) {
    if (x.empty()) return {};
    if (n <= 1) return x;

    size_t len = x.size();
    std::vector<double> y(len);

    int n_minus = n / 2;
    int n_plus = (n % 2 == 0) ? (n / 2 - 1) : (n / 2);

    std::vector<double> window(n);   // <-- allocated ONCE, outside the loop

    for (int i = 0; i < (int)len; ++i) {
        int k = 0;
        for (int j = i - n_minus; j <= i + n_plus; ++j, ++k) {
            window[k] = (j >= 0 && j < (int)len) ? x[j] : 0.0;
        }
        std::sort(window.begin(), window.end());
        if (n % 2 != 0)
            y[i] = window[n / 2];
        else
            y[i] = (window[n / 2 - 1] + window[n / 2]) / 2.0;
    }
    return y;
}

// ============================================================================
// detrend: subtract a fitted trend from x. mode="constant" removes the mean;
// mode="linear" (default) fits a least-squares line and subtracts it. Matches
// MATLAB detrend(x) / detrend(x, 'constant'). NaN-aware: NaN samples are
// excluded from the fit, and are returned as NaN in the output.
// ============================================================================
inline std::vector<double> detrend(const std::vector<double>& x, const std::string& mode = "linear") {
    const size_t n = x.size();
    std::vector<double> y(n, std::numeric_limits<double>::quiet_NaN());
    if (n == 0) return y;

    if (mode == "constant") {
        double s = 0.0; size_t cnt = 0;
        for (double v : x) if (!std::isnan(v)) { s += v; ++cnt; }
        const double m = (cnt > 0) ? s / cnt : 0.0;
        for (size_t i = 0; i < n; ++i)
            if (!std::isnan(x[i])) y[i] = x[i] - m;
        return y;
    }

    // Linear least-squares fit: y = m*i + c over non-NaN samples.
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    size_t cnt = 0;
    for (size_t i = 0; i < n; ++i) {
        if (std::isnan(x[i])) continue;
        const double xi = static_cast<double>(i);
        sx += xi; sy += x[i]; sxx += xi * xi; sxy += xi * x[i]; ++cnt;
    }
    if (cnt < 2) {
        // Degenerate: fall back to constant detrend (or copy through if empty).
        const double m = (cnt == 1) ? sy : 0.0;
        for (size_t i = 0; i < n; ++i)
            if (!std::isnan(x[i])) y[i] = x[i] - m;
        return y;
    }
    const double dcnt = static_cast<double>(cnt);
    const double denom = dcnt * sxx - sx * sx;
    const double slope = (std::fabs(denom) > 1e-15) ? (dcnt * sxy - sx * sy) / denom : 0.0;
    const double intercept = (sy - slope * sx) / dcnt;
    for (size_t i = 0; i < n; ++i)
        if (!std::isnan(x[i])) y[i] = x[i] - (slope * static_cast<double>(i) + intercept);
    return y;
}

// ============================================================================
// notch_filter: remove a narrow band around notch_hz from x.
//
// Implemented as (x - narrow_bandpass(x)) using the ZPK bandpass overload,
// applied zero-phase with filtfilt. Q sets the notch width: bandwidth =
// notch_hz / Q. Q=30 gives ~2 Hz total width at 60 Hz -- tight enough to
// kill powerline without smearing neighboring content. Order N applies to
// the internal bandpass; N=4 is a solid default.
//
// If notch_hz <= 0 or fs <= 0, returns the input unchanged (this is what
// callers gated on notch_filter_hz should see when the setting is disabled).
// ============================================================================
inline std::vector<double> notch_filter(const std::vector<double>& x, double notch_hz,
    double fs, double Q = 30.0, int N = 4) {
    if (x.empty() || notch_hz <= 0.0 || fs <= 0.0 || Q <= 0.0) return x;
    const double nyq = fs / 2.0;
    const double bw = notch_hz / Q;
    const double lo = notch_hz - 0.5 * bw;
    const double hi = notch_hz + 0.5 * bw;
    if (lo <= 0.0 || hi >= nyq) return x;   // out of band; no-op

    // NaN-aware: filtfilt is IIR + zero-phase, so a single NaN in the input
    // contaminates every output sample. Real templates commonly have NaN
    // tails (column-median-over-variable-length beats), so plain filtfilt
    // would return all-NaN for almost every bin. Replace NaNs with
    // linearly-interpolated placeholders for the filter run and restore
    // them at the end so the caller sees exactly the same NaN mask.
    std::vector<bool> was_nan(x.size(), false);
    std::vector<double> xi = x;
    size_t total_valid = 0;
    for (size_t i = 0; i < xi.size(); ++i) {
        if (std::isnan(xi[i])) was_nan[i] = true;
        else ++total_valid;
    }
    if (total_valid == 0) return x;   // nothing to filter; every sample was NaN
    if (total_valid < xi.size()) {
        // Forward pass: any leading NaNs get the first valid sample.
        double last = std::numeric_limits<double>::quiet_NaN();
        for (size_t i = 0; i < xi.size(); ++i)
            if (!was_nan[i]) { last = xi[i]; break; }
        for (size_t i = 0; i < xi.size(); ++i) {
            if (was_nan[i]) xi[i] = last;
            else last = xi[i];
        }
        // Backward pass: any remaining leading-run (never had a valid sample
        // to their left originally, but the forward pass already fixed
        // them). This second pass handles interior NaN runs by linearly
        // interpolating between the valid neighbors on either side, so we
        // don't inject a flat step that ripples through the filter.
        for (size_t i = 0; i < xi.size(); ++i) {
            if (!was_nan[i]) continue;
            // Find end of this NaN run.
            size_t j = i;
            while (j < xi.size() && was_nan[j]) ++j;
            // xi[i-1] is the last valid before, xi[j] is the first valid
            // after (if any). Linear interp between them for the whole run.
            if (i > 0 && j < xi.size()) {
                const double a = xi[i - 1];
                const double b = xi[j];
                const double denom = static_cast<double>(j - (i - 1));
                for (size_t k = i; k < j; ++k)
                    xi[k] = a + (b - a) * static_cast<double>(k - (i - 1)) / denom;
            }
            // else: trailing NaN run past the last valid sample; the forward
            // pass already filled these with the last valid value, which
            // is the best we can do without extrapolation.
            i = j;
        }
    }

    std::vector<double> b, a;
    std::vector<double> Wn = { lo / nyq, hi / nyq };
    butter(N, Wn, b, a);
    const auto bp = filtfilt(b, a, xi);

    std::vector<double> y(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        if (was_nan[i]) y[i] = std::numeric_limits<double>::quiet_NaN();
        else            y[i] = xi[i] - bp[i];
    }
    return y;
}

// ============================================================================
// waveform_highpass: high-pass a waveform at cutoff_hz, applied zero-phase
// via filtfilt. Spec: disabled by default (returns input unchanged if
// cutoff_hz <= 0), 0.5 Hz when enabled. Applied before template averaging.
// Order default 3 matches the LP/HP default in the rest of this codebase.
// ============================================================================
inline std::vector<double> waveform_highpass(const std::vector<double>& x,
    double cutoff_hz, double fs, int N = 3) {
    if (x.empty() || cutoff_hz <= 0.0 || fs <= 0.0) return x;
    const double Wn = cutoff_hz / (fs / 2.0);
    if (Wn <= 0.0 || Wn >= 1.0) return x;

    // NaN-aware, same reasoning as notch_filter above: IIR + zero-phase
    // means one NaN contaminates the whole output. Interpolate NaN runs
    // linearly for the filter, restore the NaN mask at the end.
    std::vector<bool> was_nan(x.size(), false);
    std::vector<double> xi = x;
    size_t total_valid = 0;
    for (size_t i = 0; i < xi.size(); ++i) {
        if (std::isnan(xi[i])) was_nan[i] = true;
        else ++total_valid;
    }
    if (total_valid == 0) return x;
    if (total_valid < xi.size()) {
        double last = std::numeric_limits<double>::quiet_NaN();
        for (size_t i = 0; i < xi.size(); ++i)
            if (!was_nan[i]) { last = xi[i]; break; }
        for (size_t i = 0; i < xi.size(); ++i) {
            if (was_nan[i]) xi[i] = last;
            else last = xi[i];
        }
        for (size_t i = 0; i < xi.size(); ++i) {
            if (!was_nan[i]) continue;
            size_t j = i;
            while (j < xi.size() && was_nan[j]) ++j;
            if (i > 0 && j < xi.size()) {
                const double a = xi[i - 1];
                const double b = xi[j];
                const double denom = static_cast<double>(j - (i - 1));
                for (size_t k = i; k < j; ++k)
                    xi[k] = a + (b - a) * static_cast<double>(k - (i - 1)) / denom;
            }
            i = j;
        }
    }

    std::vector<double> b, a;
    butter(N, Wn, "high", b, a);
    auto y = filtfilt(b, a, xi);
    for (size_t i = 0; i < x.size(); ++i)
        if (was_nan[i]) y[i] = std::numeric_limits<double>::quiet_NaN();
    return y;
}