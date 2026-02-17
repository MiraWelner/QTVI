// ============================================================================
// File: FilterUtils.cpp
// Comprehensive implementation of MATLAB-equivalent DSP functions
// ============================================================================
#include "FilterUtils.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <set>

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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
    vector<double> z(n_order, 0.0); // Filter state
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
// 2. Zero-Phase Filtering (Mirror Padding)
// Matches: y = filtfilt(b, a, x)
// ============================================================================
vector<double> filtfilt(const vector<double>& b, const vector<double>& a, const vector<double>& x) {
    if (x.empty()) return x;
    int nfact = 3 * (max((int)b.size(), (int)a.size()) - 1);
    if (x.size() <= (size_t)nfact) return filter(b, a, x);

    // Reflective padding to match MATLAB's transient suppression
    vector<double> padded;
    padded.reserve(x.size() + 2 * nfact);
    for (int i = nfact; i > 0; --i) padded.push_back(2.0 * x[0] - x[i]);
    padded.insert(padded.end(), x.begin(), x.end());
    for (int i = 1; i <= nfact; ++i) padded.push_back(2.0 * x.back() - x[x.size() - 1 - i]);

    // Forward and backward filtering
    vector<double> y = filter(b, a, padded);
    reverse(y.begin(), y.end());
    y = filter(b, a, y);
    reverse(y.begin(), y.end());

    // Extract original signal from padding
    return vector<double>(y.begin() + nfact, y.end() - nfact);
}

// ============================================================================
// 3. Butterworth Filter Design (3rd Order Bilinear Transform)
// Matches: [b, a] = butter(3, Wn, 'type')
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

// Bandpass design (Cascade)
void butter(int N, const vector<double>& Wn, vector<double>& b, vector<double>& a) {
    vector<double> b_high, a_high, b_low, a_low;
    butter(N, Wn[0], "high", b_high, a_high);
    butter(N, Wn[1], "low", b_low, a_low);
    b = conv(b_high, b_low);
    a = conv(a_high, a_low);
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
vector<double> medfilt1(const vector<double>& x, int window_size) {
    if (x.empty()) return {};
    if (window_size <= 1) return x;

    vector<double> y(x.size());
    int half = window_size / 2;

    for (int i = 0; i < (int)x.size(); ++i) {
        vector<double> window;
        for (int j = -half; j <= half; ++j) {
            int idx = i + j;
            // MATLAB uses zero-padding for medfilt1 boundaries
            window.push_back((idx < 0 || idx >= (int)x.size()) ? 0.0 : x[idx]);
        }
        // Use nth_element for efficient median calculation O(n)
        nth_element(window.begin(), window.begin() + window.size() / 2, window.end());
        y[i] = window[window.size() / 2];
    }
    return y;
}
