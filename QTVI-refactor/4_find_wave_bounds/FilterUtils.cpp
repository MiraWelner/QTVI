// ============================================================================
// File: FilterUtils.cpp
// Comprehensive implementation of MATLAB-equivalent DSP functions
// ============================================================================
#include "FilterUtils.h"
#include <cmath>
#include <complex>
#include <algorithm>
#include <vector>
#include <set>

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef std::complex<double> cd;

// ============================================================================
// Helper: polynomial from complex roots
// Matches: poly(roots) in MATLAB
// ============================================================================
static vector<double> poly_from_roots(const vector<cd>& roots) {
    int n = (int)roots.size();
    vector<cd> c(n + 1, cd(0.0, 0.0));
    c[0] = cd(1.0, 0.0);
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j >= 1; j--) {
            c[j] -= roots[i] * c[j - 1];
        }
    }
    vector<double> result(n + 1);
    for (int i = 0; i <= n; i++) {
        result[i] = c[i].real();
    }
    return result;
}

// ============================================================================
// 1. Digital IIR Filter (Direct Form II Transposed)
// Matches: y = filter(b, a, x)
// ============================================================================
vector<double> filter(const vector<double>& b, const vector<double>& a, const vector<double>& x) {
    if (x.empty() || a.empty()) return {};
    size_t n = x.size();
    size_t nb = b.size();
    size_t na = a.size();
    size_t n_order = max(nb, na);

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
// Helper: filter with initial conditions (Direct Form II Transposed)
// Matches: y = filter(b, a, x, zi)
// ============================================================================
static vector<double> filter_ic(const vector<double>& b, const vector<double>& a,
    const vector<double>& x, vector<double> zi) {
    if (x.empty() || a.empty()) return {};
    size_t n = x.size();
    size_t nb = b.size();
    size_t na = a.size();
    size_t n_order = max(nb, na);
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

// ============================================================================
// Helper: Compute initial conditions for filtfilt
// Matches MATLAB's internal filtic-like computation:
//   solve (I - A) * zi = B  where A and B are built from filter coeffs
// ============================================================================
static vector<double> compute_filtfilt_ic(const vector<double>& b, const vector<double>& a) {
    int n_order = (int)max(b.size(), a.size());
    int m = n_order - 1;
    if (m <= 0) return {};

    // Normalize coefficients
    double a0 = a[0];
    vector<double> bn(n_order, 0.0), an(n_order, 0.0);
    for (int i = 0; i < (int)b.size() && i < n_order; i++) bn[i] = b[i] / a0;
    for (int i = 0; i < (int)a.size() && i < n_order; i++) an[i] = a[i] / a0;

    // Build matrix (I - A) where A is the companion-form matrix
    // and vector B = b[1:] - a[1:]*b[0]
    // This solves for zi such that the filter output is constant
    // when the input is constant (the DC steady-state condition).
    //
    // MATLAB equivalent: 
    //   rows = [1:m-1; 1:m-1];  cols = [2:m; 1:m-1];
    //   vals = [ones(1,m-1); -ones(1,m-1)];
    //   rhs = b(2:end) - b(1)*a(2:end);
    //   zi = (eye(m) + diag(a(2:m),1)) \ rhs;

    // Simple Gaussian elimination for small system (m x m)
    // Build the system: mat * zi = rhs
    vector<vector<double>> mat(m, vector<double>(m, 0.0));
    vector<double> rhs(m, 0.0);

    // First row
    mat[0][0] = 1.0 + an[1];
    for (int j = 1; j < m; j++) {
        mat[0][j] = an[j + 1];
    }
    rhs[0] = bn[1] - an[1] * bn[0];

    // Rows 2 to m-1
    for (int i = 1; i < m; i++) {
        mat[i][i] = 1.0;
        mat[i][i - 1] = -1.0;
        if (i + 1 < n_order) {
            rhs[i] = bn[i + 1] - an[i + 1] * bn[0];
        }
        for (int j = i + 1; j < m; j++) {
            // These are zero already
        }
        // Add a[j+1] contributions to first equation style
        // Actually the MATLAB system is:
        // (I + A') zi = B where A' has a(2..m) in first row, -1 on subdiagonal
    }

    // Actually let me implement this more carefully matching MATLAB exactly.
    // MATLAB builds: 
    //   A = eye(m) + [a(2) a(3) ... a(m) 0; -eye(m-1) zeros(m-1,1)]  (wrong)
    // The actual MATLAB code from filtfilt.m is:
    //   zi = ( eye(nfilt-1) - [-a(2:nfilt)'; eye(nfilt-2) zeros(nfilt-2,1)] ) \ 
    //        ( b(2:nfilt)' - b(1)*a(2:nfilt)' );

    // Build matrix matching MATLAB's sparse construction in filtfilt.m:
    //   rows = [1:m, 2:m, 1:m-1];
    //   cols = [ones(1,m), 2:m, 2:m];
    //   vals = [1+a(2), a(3:nfilt), ones(1,m-1), -ones(1,m-1)];
    //
    // This gives:
    //   First column:      [1+a[1], a[2], a[3], ..., a[m]]
    //   Diagonal (i>=1):   1
    //   Superdiagonal:     -1
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            mat[i][j] = 0.0;
        }
    }

    // First column: a coefficients
    mat[0][0] = 1.0 + an[1];
    for (int i = 1; i < m; i++) {
        mat[i][0] = (i + 1 < n_order) ? an[i + 1] : 0.0;
    }

    // Diagonal (rows 1..m-1): 1
    for (int i = 1; i < m; i++) {
        mat[i][i] = 1.0;
    }

    // Superdiagonal (rows 0..m-2): -1
    for (int i = 0; i < m - 1; i++) {
        mat[i][i + 1] = -1.0;
    }

    // RHS: b(2:nfilt) - b(1)*a(2:nfilt)
    for (int i = 0; i < m; i++) {
        double bi = (i + 1 < (int)b.size()) ? b[i + 1] / a0 : 0.0;
        double ai = (i + 1 < (int)a.size()) ? a[i + 1] / a0 : 0.0;
        rhs[i] = bi - bn[0] * ai;
    }

    // Gaussian elimination with partial pivoting
    for (int col = 0; col < m; col++) {
        // Find pivot
        int pivot = col;
        double maxval = fabs(mat[col][col]);
        for (int row = col + 1; row < m; row++) {
            if (fabs(mat[row][col]) > maxval) {
                maxval = fabs(mat[row][col]);
                pivot = row;
            }
        }
        if (pivot != col) {
            swap(mat[col], mat[pivot]);
            swap(rhs[col], rhs[pivot]);
        }
        if (fabs(mat[col][col]) < 1e-15) continue;

        for (int row = col + 1; row < m; row++) {
            double factor = mat[row][col] / mat[col][col];
            for (int j = col; j < m; j++) {
                mat[row][j] -= factor * mat[col][j];
            }
            rhs[row] -= factor * rhs[col];
        }
    }

    // Back substitution
    vector<double> zi(m, 0.0);
    for (int i = m - 1; i >= 0; i--) {
        double sum = rhs[i];
        for (int j = i + 1; j < m; j++) {
            sum -= mat[i][j] * zi[j];
        }
        if (fabs(mat[i][i]) > 1e-15)
            zi[i] = sum / mat[i][i];
    }

    return zi;
}

// ============================================================================
// 2. Zero-Phase Filtering (Mirror Padding + Initial Conditions)
// Matches: y = filtfilt(b, a, x)
// ============================================================================
vector<double> filtfilt(const vector<double>& b, const vector<double>& a, const vector<double>& x) {
    if (x.empty()) return x;
    if (a.empty() || b.empty()) return x;

    int nfact = 3 * (max((int)b.size(), (int)a.size()) - 1);

    // CRITICAL: If signal is too short for padding, just do regular filtering
    if ((int)x.size() <= nfact) {
        return filter(b, a, x);
    }

    // Compute initial condition unit vector
    vector<double> zi_unit = compute_filtfilt_ic(b, a);
    int m = (int)zi_unit.size();

    // Reflective padding
    vector<double> padded;
    padded.reserve(x.size() + 2 * nfact);

    for (int i = nfact; i > 0; --i) {
        padded.push_back(2.0 * x[0] - x[i]);
    }
    padded.insert(padded.end(), x.begin(), x.end());
    for (int i = 1; i <= nfact; ++i) {
        padded.push_back(2.0 * x.back() - x[x.size() - 1 - i]);
    }

    // Forward pass with initial conditions scaled by first padded sample
    vector<double> zi_fwd(m);
    for (int i = 0; i < m; i++) zi_fwd[i] = zi_unit[i] * padded[0];
    vector<double> y = filter_ic(b, a, padded, zi_fwd);

    // Reverse
    std::reverse(y.begin(), y.end());

    // Backward pass with initial conditions scaled by first sample of reversed
    vector<double> zi_bwd(m);
    for (int i = 0; i < m; i++) zi_bwd[i] = zi_unit[i] * y[0];
    y = filter_ic(b, a, y, zi_bwd);

    // Reverse back
    std::reverse(y.begin(), y.end());

    // Extract middle
    if ((int)y.size() < 2 * nfact) return x;
    return vector<double>(y.begin() + nfact, y.end() - nfact);
}

// ============================================================================
// 3a. Butterworth Lowpass/Highpass Design (3rd Order)
// Matches: [b, a] = butter(3, Wn, 'low') or butter(3, Wn, 'high')
// ============================================================================
void butter(int N, double Wn, const string& type, vector<double>& b, vector<double>& a) {
    double alpha = tan(M_PI * Wn / 2.0);
    double a2 = alpha * alpha;
    double a3 = alpha * a2;

    if (type == "low") {
        double denom = 1.0 + 2.0 * alpha + 2.0 * a2 + a3;
        b = { a3 / denom, 3.0 * a3 / denom, 3.0 * a3 / denom, a3 / denom };
        a = { 1.0,
             (3.0 * a3 + 2.0 * a2 - 2.0 * alpha - 3.0) / denom,
             (3.0 * a3 - 2.0 * a2 - 2.0 * alpha + 3.0) / denom,
             (a3 - 2.0 * a2 + 2.0 * alpha - 1.0) / denom };
    }
    else if (type == "high") {
        double denom = 1.0 + 2.0 * alpha + 2.0 * a2 + a3;
        b = { 1.0 / denom, -3.0 / denom, 3.0 / denom, -1.0 / denom };
        a = { 1.0,
             (3.0 + 2.0 * alpha - 2.0 * a2 - 3.0 * a3) / denom,
             (3.0 - 2.0 * alpha - 2.0 * a2 + 3.0 * a3) / denom,
             (1.0 - 2.0 * alpha + 2.0 * a2 - a3) / denom };
    }
}

// ============================================================================
// 3b. Butterworth Bandpass Design (Proper ZPK implementation)
// Matches: [b, a] = butter(N, [Wn1 Wn2])
//
// Algorithm: buttap -> lp2bp -> bilinear -> zp2tf
// This is a faithful replication of MATLAB's internal pipeline.
// ============================================================================
void butter(int N, const vector<double>& Wn, vector<double>& b, vector<double>& a) {
    if (Wn.size() < 2) return;

    const double fs = 2.0;  // MATLAB convention

    // Pre-warp cutoff frequencies to analog domain
    // MATLAB: u = 2*fs*tan(pi*Wn/fs)
    double u1 = 2.0 * fs * tan(M_PI * Wn[0] / fs);
    double u2 = 2.0 * fs * tan(M_PI * Wn[1] / fs);
    double BW = u2 - u1;
    double W0 = sqrt(u1 * u2);

    // ---- buttap(N): Analog Butterworth prototype ----
    // Poles on unit circle in left-half s-plane, no zeros, gain = 1
    vector<cd> proto_poles;
    for (int k = 1; k <= N; k++) {
        double angle = M_PI * (2.0 * k + N - 1.0) / (2.0 * N);
        proto_poles.push_back(cd(cos(angle), sin(angle)));
    }

    // ---- lp2bp: Lowpass to Bandpass transformation ----
    // Each prototype pole p maps to two bandpass poles:
    //   s = p*BW/2 +/- sqrt((p*BW/2)^2 - W0^2)
    // N zeros are added at s = 0
    // Gain becomes k * BW^N
    vector<cd> bp_poles;
    for (const auto& p : proto_poles) {
        cd pBW2 = p * BW / 2.0;
        cd sq = sqrt(pBW2 * pBW2 - W0 * W0);
        bp_poles.push_back(pBW2 + sq);
        bp_poles.push_back(pBW2 - sq);
    }

    vector<cd> bp_zeros(N, cd(0.0, 0.0));
    double bp_gain = pow(BW, N);

    // ---- bilinear: Analog to Digital ----
    // MATLAB: kd = k * real(prod(2*fs - z) / prod(2*fs - p))
    cd num_prod(1.0, 0.0);
    for (const auto& z : bp_zeros) {
        num_prod *= (2.0 * fs - z);
    }
    cd den_prod(1.0, 0.0);
    for (const auto& p : bp_poles) {
        den_prod *= (2.0 * fs - p);
    }
    double d_gain = bp_gain * (num_prod / den_prod).real();

    // Map poles: zd = (1 + s/(2*fs)) / (1 - s/(2*fs))
    vector<cd> d_poles, d_zeros;
    for (const auto& p : bp_poles) {
        d_poles.push_back((1.0 + p / (2.0 * fs)) / (1.0 - p / (2.0 * fs)));
    }
    for (const auto& z : bp_zeros) {
        d_zeros.push_back((1.0 + z / (2.0 * fs)) / (1.0 - z / (2.0 * fs)));
    }

    // Extra zeros at z = -1 (from s = infinity)
    while (d_zeros.size() < d_poles.size()) {
        d_zeros.push_back(cd(-1.0, 0.0));
    }

    // ---- zp2tf: Zeros/Poles/Gain to Transfer Function ----
    a = poly_from_roots(d_poles);
    b = poly_from_roots(d_zeros);

    // Apply gain and normalize so a[0] = 1
    double a0 = a[0];
    for (auto& v : b) v *= d_gain / a0;
    for (auto& v : a) v /= a0;
}

// ============================================================================
// 4. Convolution
// Matches: y = conv(a, b)
// ============================================================================
vector<double> conv(const vector<double>& a, const vector<double>& b) {
    if (a.empty() || b.empty()) return {};
    vector<double> res(a.size() + b.size() - 1, 0.0);
    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < b.size(); ++j) {
            res[i + j] += a[i] * b[j];
        }
    }
    return res;
}

// ============================================================================
// 5. Median Filter
// Matches: y = medfilt1(x, window_size)
// ============================================================================
std::vector<double> medfilt1(const std::vector<double>& x, int n) {
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
            if (j >= 0 && j < (int)len) {
                window.push_back(x[j]);
            }
            else {
                window.push_back(0.0);
            }
        }

        std::sort(window.begin(), window.end());

        if (n % 2 != 0) {
            y[i] = window[n / 2];
        }
        else {
            y[i] = (window[n / 2 - 1] + window[n / 2]) / 2.0;
        }
    }
    return y;
}