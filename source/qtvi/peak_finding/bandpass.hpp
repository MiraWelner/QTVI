/**
 * @file   bandpass.hpp
 * @brief  Butterworth filter via cascaded biquads with filtfilt
 */
#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/*
* These are the biquad coefficeints
* Second order biquad coefficients are used in digital biquad filters, 
* which are recursive filters characterized by two poles and two zeros.
*/
struct Biquad {
    double b0, b1, b2;
    double a1, a2;  // a0 is always 1.0 (normalized)
};

namespace bandpass_detail {

    inline std::vector<double> butterworth_poles(int order) {
        std::vector<double> angles;
        for (int k = 0; k < order; k++)
            angles.push_back(M_PI * (2.0 * k + order + 1.0) / (2.0 * order));
        return angles;
    }

} // namespace bandpass_detail

// ============================================================================
// Lowpass design via bilinear transform - cascaded second-order sections
// ============================================================================
inline std::vector<Biquad> butterworth_lowpass(int order, double cutoff_hz, double sample_rate) {
    if (order < 1) throw std::invalid_argument("Order must be >= 1");
    if (cutoff_hz <= 0 || cutoff_hz >= sample_rate / 2.0)
        throw std::invalid_argument("Cutoff must be in (0, Nyquist)");

    double wc = std::tan(M_PI * cutoff_hz / sample_rate);
    auto angles = bandpass_detail::butterworth_poles(order);
    std::vector<Biquad> sections;

    int i = 0;
    while (i < order) {
        Biquad bq;
        if (i + 1 < order) {
            double re = std::cos(angles[i]);
            double A = 1.0, B = -2.0 * re * wc, C = wc * wc;
            double a0 = A + B + C;
            double a1_coeff = -2.0 * A + 2.0 * C;
            double a2_coeff = A - B + C;
            double b0 = C, b1_val = 2.0 * C, b2 = C;

            bq.b0 = b0 / a0; bq.b1 = b1_val / a0; bq.b2 = b2 / a0;
            bq.a1 = a1_coeff / a0; bq.a2 = a2_coeff / a0;
            sections.push_back(bq);
            i += 2;
        }
        else {
            double a0 = 1.0 + wc;
            bq.b0 = wc / a0; bq.b1 = wc / a0; bq.b2 = 0.0;
            bq.a1 = (-1.0 + wc) / a0; bq.a2 = 0.0;
            sections.push_back(bq);
            i += 1;
        }
    }
    return sections;
}

// ============================================================================
// Highpass design via bilinear transform
// ============================================================================
inline std::vector<Biquad> butterworth_highpass(int order, double cutoff_hz, double sample_rate) {
    if (order < 1) throw std::invalid_argument("Order must be >= 1");
    if (cutoff_hz <= 0 || cutoff_hz >= sample_rate / 2.0)
        throw std::invalid_argument("Cutoff must be in (0, Nyquist)");

    double wc = std::tan(M_PI * cutoff_hz / sample_rate);
    auto angles = bandpass_detail::butterworth_poles(order);
    std::vector<Biquad> sections;

    int i = 0;
    while (i < order) {
        Biquad bq;
        if (i + 1 < order) {
            double re = std::cos(angles[i]);
            double A = 1.0, B = -2.0 * re * wc, C = wc * wc;
            double a0 = A + B + C;
            double a1_coeff = -2.0 * A + 2.0 * C;
            double a2_coeff = A - B + C;
            double b0 = A, b1_val = -2.0 * A, b2 = A;

            bq.b0 = b0 / a0; bq.b1 = b1_val / a0; bq.b2 = b2 / a0;
            bq.a1 = a1_coeff / a0; bq.a2 = a2_coeff / a0;
            sections.push_back(bq);
            i += 2;
        }
        else {
            double a0 = 1.0 + wc;
            bq.b0 = 1.0 / a0; bq.b1 = -1.0 / a0; bq.b2 = 0.0;
            bq.a1 = (-1.0 + wc) / a0; bq.a2 = 0.0;
            sections.push_back(bq);
            i += 1;
        }
    }
    return sections;
}

// ============================================================================
// Apply single biquad — Direct Form II Transposed
// ============================================================================
inline std::vector<double> apply_biquad(const Biquad& bq, const std::vector<double>& x) {
    size_t n = x.size();
    std::vector<double> y(n);
    double z1 = 0.0, z2 = 0.0;
    for (size_t i = 0; i < n; i++) {
        double in = x[i];
        double out = bq.b0 * in + z1;
        z1 = bq.b1 * in - bq.a1 * out + z2;
        z2 = bq.b2 * in - bq.a2 * out;
        y[i] = out;
    }
    return y;
}

// ============================================================================
// Apply cascaded second-order sections forward
// ============================================================================
inline std::vector<double> apply_sos(const std::vector<Biquad>& sos, const std::vector<double>& x) {
    std::vector<double> y = x;
    for (const auto& bq : sos)
        y = apply_biquad(bq, y);
    return y;
}

// ============================================================================
// filtfilt - zero-phase filtering (forward + reverse) for biquad cascade
// ============================================================================
inline std::vector<double> filtfilt(const std::vector<Biquad>& sos, const std::vector<double>& x) {
    if (x.size() < 4) return x;

    size_t pad_len = 3 * sos.size();
    if (pad_len >= x.size()) pad_len = x.size() - 1;

    std::vector<double> padded(pad_len + x.size() + pad_len);

    for (size_t i = 0; i < pad_len; i++)
        padded[i] = 2.0 * x[0] - x[pad_len - i];
    for (size_t i = 0; i < x.size(); i++)
        padded[pad_len + i] = x[i];
    for (size_t i = 0; i < pad_len; i++)
        padded[pad_len + x.size() + i] = 2.0 * x.back() - x[x.size() - 2 - i];

    std::vector<double> y = apply_sos(sos, padded);
    std::reverse(y.begin(), y.end());
    y = apply_sos(sos, y);
    std::reverse(y.begin(), y.end());

    std::vector<double> result(x.size());
    for (size_t i = 0; i < x.size(); i++)
        result[i] = y[pad_len + i];
    return result;
}

// ============================================================================
// Convenience: bandpass via cascaded HP + LP with filtfilt
// ============================================================================
inline std::vector<double> bandpass_filtfilt(
    int order, double low_hz, double high_hz, double sample_rate,
    const std::vector<double>& x)
{
    auto hp = butterworth_highpass(order, low_hz, sample_rate);
    auto lp = butterworth_lowpass(order, high_hz, sample_rate);
    auto y = filtfilt(hp, x);
    return filtfilt(lp, y);
}