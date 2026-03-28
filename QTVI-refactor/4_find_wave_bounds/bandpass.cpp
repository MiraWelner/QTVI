
#include "bandpass.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================================================================
// Internal: compute analog pole angles for Butterworth of given order
// ============================================================================
static std::vector<double> butterworth_poles(int order) {
    std::vector<double> angles;
    for (int k = 0; k < order; k++) {
        angles.push_back(M_PI * (2.0 * k + order + 1.0) / (2.0 * order));
    }
    return angles;
}

// ============================================================================
// Lowpass design via bilinear transform
// Returns cascaded second-order sections
// ============================================================================
std::vector<Biquad> butterworth_lowpass(int order, double cutoff_hz, double sample_rate) {
    if (order < 1) throw std::invalid_argument("Order must be >= 1");
    if (cutoff_hz <= 0 || cutoff_hz >= sample_rate / 2.0)
        throw std::invalid_argument("Cutoff must be in (0, Nyquist)");

    // Pre-warp
    double wc = tan(M_PI * cutoff_hz / sample_rate);

    auto angles = butterworth_poles(order);
    std::vector<Biquad> sections;

    int i = 0;
    while (i < order) {
        Biquad bq;

        if (i + 1 < order) {
            // Conjugate pair -> second-order section
            double re = cos(angles[i]);
            // For Butterworth, conjugate pairs share the same real part
            // analog prototype: s^2 - 2*re*s + 1
            // bilinear with pre-warp:
            double a = wc * wc;
            double b_coeff = -2.0 * re * wc;
            double c = 1.0;

            double norm = a + b_coeff + c;
            // bilinear: s = (2/T)(z-1)/(z+1), T=1 for normalized
            // Direct computation for 2nd order section:

            double d = 1.0 + b_coeff + a;  // = c + b + a with c=1
            // Actually let's do this properly with bilinear transform

            // Analog second order: H(s) = wc^2 / (s^2 - 2*re*wc*s + wc^2)
            // Bilinear: s = (1 - z^-1) / (1 + z^-1)  (pre-warped wc already in cutoff)

            double A = 1.0;
            double B = -2.0 * re * wc;
            double C = wc * wc;

            // Numerator of analog = C (constant)
            // Denominator = A*s^2 + B*s + C

            // After bilinear s = (1-z^-1)/(1+z^-1):
            // s^2 = (1 - 2z^-1 + z^-2) / (1 + 2z^-1 + z^-2)
            // s   = (1 - z^-1) / (1 + z^-1)

            // Denominator in z:
            // A*(1 - 2z^-1 + z^-2) + B*(1 - z^-1)*(1 + z^-1)... no
            // multiply through by (1+z^-1)^2:
            // A(1-z^-1)^2 + B(1-z^-1)(1+z^-1) + C(1+z^-1)^2

            double a0 = A + B + C;
            double a1_coeff = -2.0 * A + 2.0 * C;
            double a2_coeff = A - B + C;

            // Numerator: C * (1 + z^-1)^2 = C*(1 + 2z^-1 + z^-2)
            double b0 = C;
            double b1_val = 2.0 * C;
            double b2 = C;

            // Normalize
            bq.b0 = b0 / a0;
            bq.b1 = b1_val / a0;
            bq.b2 = b2 / a0;
            bq.a1 = a1_coeff / a0;
            bq.a2 = a2_coeff / a0;

            sections.push_back(bq);
            i += 2;
        }
        else {
            // Odd order: single real pole
            // Analog: H(s) = wc / (s + wc)
            // Bilinear: s = (1 - z^-1)/(1 + z^-1)
            // H(z) = wc*(1 + z^-1) / ((1 - z^-1) + wc*(1 + z^-1))
            //       = wc*(1 + z^-1) / ((1+wc) + (-1+wc)*z^-1)

            double a0 = 1.0 + wc;
            bq.b0 = wc / a0;
            bq.b1 = wc / a0;
            bq.b2 = 0.0;
            bq.a1 = (-1.0 + wc) / a0;
            bq.a2 = 0.0;

            sections.push_back(bq);
            i += 1;
        }
    }

    return sections;
}

// ============================================================================
// Highpass design via bilinear transform
// ============================================================================
std::vector<Biquad> butterworth_highpass(int order, double cutoff_hz, double sample_rate) {
    if (order < 1) throw std::invalid_argument("Order must be >= 1");
    if (cutoff_hz <= 0 || cutoff_hz >= sample_rate / 2.0)
        throw std::invalid_argument("Cutoff must be in (0, Nyquist)");

    // Pre-warp
    double wc = tan(M_PI * cutoff_hz / sample_rate);

    auto angles = butterworth_poles(order);
    std::vector<Biquad> sections;

    int i = 0;
    while (i < order) {
        Biquad bq;

        if (i + 1 < order) {
            // Conjugate pair
            double re = cos(angles[i]);

            // Analog highpass prototype: H(s) = s^2 / (s^2 - 2*re*wc*s + wc^2)
            // Bilinear with s = (1 - z^-1)/(1 + z^-1):

            double A = 1.0;
            double B = -2.0 * re * wc;
            double C = wc * wc;

            // Denominator in z (same as lowpass):
            double a0 = A + B + C;
            double a1_coeff = -2.0 * A + 2.0 * C;
            double a2_coeff = A - B + C;

            // Numerator: A * (1 - z^-1)^2 = A*(1 - 2z^-1 + z^-2)
            double b0 = A;
            double b1_val = -2.0 * A;
            double b2 = A;

            bq.b0 = b0 / a0;
            bq.b1 = b1_val / a0;
            bq.b2 = b2 / a0;
            bq.a1 = a1_coeff / a0;
            bq.a2 = a2_coeff / a0;

            sections.push_back(bq);
            i += 2;
        }
        else {
            // Odd order: single real pole
            // Analog highpass: H(s) = s / (s + wc)
            // H(z) = (1-z^-1) / ((1+wc) + (-1+wc)*z^-1)

            double a0 = 1.0 + wc;
            bq.b0 = 1.0 / a0;
            bq.b1 = -1.0 / a0;
            bq.b2 = 0.0;
            bq.a1 = (-1.0 + wc) / a0;
            bq.a2 = 0.0;

            sections.push_back(bq);
            i += 1;
        }
    }

    return sections;
}

// ============================================================================
// Apply single biquad — Direct Form II Transposed
// ============================================================================
std::vector<double> apply_biquad(const Biquad& bq, const std::vector<double>& x) {
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
std::vector<double> apply_sos(const std::vector<Biquad>& sos, const std::vector<double>& x) {
    std::vector<double> y = x;
    for (const auto& bq : sos) {
        y = apply_biquad(bq, y);
    }
    return y;
}

// ============================================================================
// filtfilt — zero-phase filtering (forward + reverse)
// Includes edge padding to reduce transient artifacts
// ============================================================================
std::vector<double> filtfilt(const std::vector<Biquad>& sos, const std::vector<double>& x) {
    if (x.size() < 4) return x;

    // Pad length: 3x the number of biquad sections (like scipy)
    size_t pad_len = 3 * sos.size();
    if (pad_len >= x.size()) pad_len = x.size() - 1;

    // Edge-reflected padding to reduce startup transients
    std::vector<double> padded(pad_len + x.size() + pad_len);

    // Front padding: reflect around x[0]
    for (size_t i = 0; i < pad_len; i++) {
        padded[i] = 2.0 * x[0] - x[pad_len - i];
    }

    // Copy signal
    for (size_t i = 0; i < x.size(); i++) {
        padded[pad_len + i] = x[i];
    }

    // Back padding: reflect around x[end]
    for (size_t i = 0; i < pad_len; i++) {
        padded[pad_len + x.size() + i] = 2.0 * x.back() - x[x.size() - 2 - i];
    }

    // Forward pass
    std::vector<double> y = apply_sos(sos, padded);

    // Reverse
    std::reverse(y.begin(), y.end());

    // Reverse pass
    y = apply_sos(sos, y);

    // Reverse back
    std::reverse(y.begin(), y.end());

    // Strip padding
    std::vector<double> result(x.size());
    for (size_t i = 0; i < x.size(); i++) {
        result[i] = y[pad_len + i];
    }

    return result;
}

// ============================================================================
// Convenience: bandpass via cascaded HP + LP with filtfilt
// ============================================================================
std::vector<double> bandpass_filtfilt(
    int order,
    double low_hz,
    double high_hz,
    double sample_rate,
    const std::vector<double>& x)
{
    auto hp = butterworth_highpass(order, low_hz, sample_rate);
    auto lp = butterworth_lowpass(order, high_hz, sample_rate);

    auto y = filtfilt(hp, x);
    return filtfilt(lp, y);
}