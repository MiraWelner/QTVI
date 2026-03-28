#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>

// ============================================================================
// Second-order section (biquad) coefficients
// ============================================================================
struct Biquad {
    double b0, b1, b2;
    double a1, a2;  // a0 is always 1.0 (normalized)
};

// ============================================================================
// Butterworth filter — supports arbitrary order via cascaded biquads
//
// Usage:
//   auto coeffs = butterworth_lowpass(order, cutoff_hz, sample_rate);
//   auto filtered = filtfilt(coeffs, signal);
//
//   // Bandpass = highpass then lowpass
//   auto hp = butterworth_highpass(order, low_hz, sample_rate);
//   auto lp = butterworth_lowpass(order, high_hz, sample_rate);
//   auto filtered = filtfilt(lp, filtfilt(hp, signal));
// ============================================================================

// Design functions — return vector of biquad sections
std::vector<Biquad> butterworth_lowpass(int order, double cutoff_hz, double sample_rate);
std::vector<Biquad> butterworth_highpass(int order, double cutoff_hz, double sample_rate);

// Apply a single biquad section forward (Direct Form II Transposed)
std::vector<double> apply_biquad(const Biquad& bq, const std::vector<double>& x);

// Apply cascaded biquad sections forward
std::vector<double> apply_sos(const std::vector<Biquad>& sos, const std::vector<double>& x);

// Zero-phase filtering (forward + reverse) — equivalent to MATLAB filtfilt
std::vector<double> filtfilt(const std::vector<Biquad>& sos, const std::vector<double>& x);

// Convenience: bandpass via cascaded highpass + lowpass, with filtfilt
std::vector<double> bandpass_filtfilt(
    int order,
    double low_hz,
    double high_hz,
    double sample_rate,
    const std::vector<double>& x);