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
    if (a.empty() || b.empty()) return x;

    int nfact = 3 * (max((int)b.size(), (int)a.size()) - 1);

    // CRITICAL: If signal is too short for padding, just do regular filtering
    if (x.size() <= (size_t)nfact) {
        return filter(b, a, x);
    }

    vector<double> padded;
    padded.reserve(x.size() + 2 * nfact);

    // Forward padding (reflective)
    for (int i = nfact; i > 0; --i) {
        padded.push_back(2.0 * x[0] - x[i]);
    }

    padded.insert(padded.end(), x.begin(), x.end());

    // Backward padding (reflective)
    for (int i = 1; i <= nfact; ++i) {
        padded.push_back(2.0 * x.back() - x[x.size() - 1 - i]);
    }

    vector<double> y = filter(b, a, padded);
    std::reverse(y.begin(), y.end());
    y = filter(b, a, y);
    std::reverse(y.begin(), y.end());

    // Safely extract the middle part
    if (y.size() < (size_t)(2 * nfact)) return x;
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
std::vector<double> medfilt1(const std::vector<double>& x, int n) {
    if (x.empty()) return {};
    if (n <= 1) return x;

    size_t len = x.size();
    std::vector<double> y(len);

    // MATLAB centered window logic:
    // If n=10, window is x(k-5 : k+4)
    // If n=11, window is x(k-5 : k+5)
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
                // MATLAB Rule: Assume 0 beyond endpoints
                window.push_back(0.0);
            }
        }

        std::sort(window.begin(), window.end());

        if (n % 2 != 0) {
            // Odd: Middle element
            y[i] = window[n / 2];
        }
        else {
            // Even: Average of two middle elements (Rule for medfilt1)
            y[i] = (window[n / 2 - 1] + window[n / 2]) / 2.0;
        }
    }
    return y;
}