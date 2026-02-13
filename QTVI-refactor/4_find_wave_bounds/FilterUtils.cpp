
// ============================================================================
// File: FilterUtils.cpp
// ============================================================================
#define _USE_MATH_DEFINES 
#include <set>
#include "FilterUtils.h"

// Helper function for Butterworth filter design (simplified)
void butter(int N, double Wn, const string& type, vector<double>& b, vector<double>& a) {
    // Simplified Butterworth filter implementation
    // For a complete implementation, consider using a DSP library
    // This is a basic approximation for low and high pass filters

    if (type == "low") {
        // Simplified low-pass Butterworth
        double c = 1.0 / std::tan(M_PI * Wn);
        double c2 = c * c;
        double sqrt2 = std::sqrt(2.0);

        if (N == 3) {
            // 3rd order Butterworth
            b.resize(4);
            a.resize(4);

            double norm = 1.0 / (1.0 + sqrt2 * c + c2);
            b[0] = norm;
            b[1] = 3 * norm;
            b[2] = 3 * norm;
            b[3] = norm;

            a[0] = 1.0;
            a[1] = (2.0 - 2.0 * c2) * norm;
            a[2] = (1.0 - sqrt2 * c + c2) * norm;
            a[3] = 0.0;
        }
    }
    else if (type == "high") {
        // Simplified high-pass Butterworth
        double c = std::tan(M_PI * Wn);
        double c2 = c * c;
        double sqrt2 = std::sqrt(2.0);

        if (N == 3) {
            // 3rd order Butterworth
            b.resize(4);
            a.resize(4);

            double norm = 1.0 / (1.0 + sqrt2 * c + c2);
            b[0] = norm;
            b[1] = -3 * norm;
            b[2] = 3 * norm;
            b[3] = -norm;

            a[0] = 1.0;
            a[1] = (2.0 * c2 - 2.0) * norm;
            a[2] = (1.0 - sqrt2 * c + c2) * norm;
            a[3] = 0.0;
        }
    }
}

// Butterworth bandpass filter
void butter(int N, const vector<double>& Wn, vector<double>& b, vector<double>& a) {
    // Simplified bandpass Butterworth filter
    // For production use, consider using a dedicated DSP library
    if (Wn.size() != 2) {
        throw std::invalid_argument("Wn must have 2 elements for bandpass filter");
    }

    double w1 = Wn[0];
    double w2 = Wn[1];
    double w0 = std::sqrt(w1 * w2);
    double bw = w2 - w1;

    // Simplified 3rd order bandpass
    if (N == 3) {
        b.resize(7);
        a.resize(7);

        // This is a simplified version - use a DSP library for accurate implementation
        double norm = 1.0;
        b[0] = bw * bw * bw * norm;
        b[1] = 0;
        b[2] = -3 * bw * bw * bw * norm;
        b[3] = 0;
        b[4] = 3 * bw * bw * bw * norm;
        b[5] = 0;
        b[6] = -bw * bw * bw * norm;

        a[0] = 1.0;
        a[1] = 0;
        a[2] = 0;
        a[3] = 0;
        a[4] = 0;
        a[5] = 0;
        a[6] = 0;
    }
}

// Apply IIR filter
vector<double> filter(const vector<double>& b, const vector<double>& a, const vector<double>& x) {
    size_t n = x.size();
    size_t nb = b.size();
    size_t na = a.size();

    vector<double> y(n, 0.0);
    vector<double> z(std::max(na, nb) - 1, 0.0);

    for (size_t i = 0; i < n; ++i) {
        // Calculate output
        y[i] = b[0] * x[i] + z[0];

        // Update state
        for (size_t j = 1; j < z.size(); ++j) {
            z[j - 1] = z[j];
            if (j < nb) z[j - 1] += b[j] * x[i];
            if (j < na) z[j - 1] -= a[j] * y[i];
        }

        // Last state element
        if (z.size() > 0) {
            z.back() = 0;
            if (nb > z.size()) z.back() += b[z.size()] * x[i];
            if (na > z.size()) z.back() -= a[z.size()] * y[i];
        }
    }

    return y;
}

// Zero-phase forward and reverse filtering
vector<double> filtfilt(const vector<double>& b, const vector<double>& a, const vector<double>& x) {
    // Forward filter
    vector<double> y = filter(b, a, x);

    // Reverse the result
    std::reverse(y.begin(), y.end());

    // Filter again
    y = filter(b, a, y);

    // Reverse back
    std::reverse(y.begin(), y.end());

    return y;
}

// Convolution
vector<double> conv(const vector<double>& a, const vector<double>& b) {
    size_t na = a.size();
    size_t nb = b.size();
    size_t nc = na + nb - 1;

    vector<double> c(nc, 0.0);

    for (size_t i = 0; i < na; ++i) {
        for (size_t j = 0; j < nb; ++j) {
            c[i + j] += a[i] * b[j];
        }
    }

    return c;
}

// Median filter
vector<double> medfilt1(const vector<double>& x, int window_size) {
    size_t n = x.size();
    if (n == 0) return {};
    vector<double> y(n);
    int half = window_size / 2;

    std::multiset<double> window;

    // Fill initial window with padding (or zeros) to match Matlab behavior
    for (int i = -half; i <= half; ++i) {
        window.insert((i < 0) ? 0.0 : x[std::min((size_t)i, n - 1)]);
    }

    auto get_median = [&]() {
        auto it = window.begin();
        std::advance(it, window.size() / 2);
        return *it;
        };

    y[0] = get_median();

    for (size_t i = 1; i < n; ++i) {
        // Remove the element that fell off the left side
        int left_idx = (int)i - half - 1;
        window.erase(window.find(left_idx < 0 ? 0.0 : x[std::min((size_t)left_idx, n - 1)]));

        // Add the new element appearing on the right side
        int right_idx = (int)i + half;
        window.insert(right_idx >= n ? 0.0 : x[right_idx]);

        y[i] = get_median();
    }

    return y;
}
