// ============================================================================
// File: FilterUtils.hpp
// Digital filtering utilities
// Includes: filter, filtfilt, butter (lowpass/highpass/bandpass), conv, medfilt1
// ============================================================================
#pragma once

#include "SignalProcessingTypes.h"
#include <complex>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace filter_detail {

    using cd = std::complex<double>;

    // Polynomial from complex roots — matches MATLAB poly(roots)
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
// Direct Form II Transposed IIR filter — matches y = filter(b, a, x)
// ============================================================================
inline vector<double> filter(const vector<double>& b, const vector<double>& a, const vector<double>& x) {
    if (x.empty() || a.empty()) return {};
    size_t n = x.size();
    size_t nb = b.size();
    size_t na = a.size();
    size_t n_order = std::max(nb, na);

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
// Butterworth lowpass/highpass design (3rd order)
// ============================================================================
inline void butter(int N, double Wn, const string& type, vector<double>& b, vector<double>& a) {
    double alpha = std::tan(M_PI * Wn / 2.0);
    double a2 = alpha * alpha;
    double a3 = alpha * a2;

    if (type == "low") {
        double d = 1.0 + 2.0 * alpha + 2.0 * a2 + a3;
        b = { a3 / d, 3.0 * a3 / d, 3.0 * a3 / d, a3 / d };
        a = { 1.0,
             (3.0 * a3 + 2.0 * a2 - 2.0 * alpha - 3.0) / d,
             (3.0 * a3 - 2.0 * a2 - 2.0 * alpha + 3.0) / d,
             (a3 - 2.0 * a2 + 2.0 * alpha - 1.0) / d };
    }
    else if (type == "high") {
        double d = 1.0 + 2.0 * alpha + 2.0 * a2 + a3;
        b = { 1.0 / d, -3.0 / d, 3.0 / d, -1.0 / d };
        a = { 1.0,
             (3.0 + 2.0 * alpha - 2.0 * a2 - 3.0 * a3) / d,
             (3.0 - 2.0 * alpha - 2.0 * a2 + 3.0 * a3) / d,
             (1.0 - 2.0 * alpha + 2.0 * a2 - a3) / d };
    }
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
// ============================================================================
inline vector<double> conv(const vector<double>& a, const vector<double>& b) {
    if (a.empty() || b.empty()) return {};
    vector<double> res(a.size() + b.size() - 1, 0.0);
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

    for (int i = 0; i < (int)len; ++i) {
        std::vector<double> window;
        window.reserve(n);
        for (int j = i - n_minus; j <= i + n_plus; ++j) {
            if (j >= 0 && j < (int)len)
                window.push_back(x[j]);
            else
                window.push_back(0.0);
        }
        std::sort(window.begin(), window.end());
        if (n % 2 != 0)
            y[i] = window[n / 2];
        else
            y[i] = (window[n / 2 - 1] + window[n / 2]) / 2.0;
    }
    return y;
}