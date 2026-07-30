#pragma once
/**
 * @file   notch_filter.hpp
 * @brief  Powerline notch (band-stop) filter, applied zero-phase. Sibling to
 *         FilterUtils.hpp (the polyphase resampler); same "no net delay"
 *         design goal, but via a forward-backward IIR pass instead of a
 *         linear-phase FIR, since a notch this narrow as an FIR would need a
 *         very long kernel.
 *
 *         Config wiring: cfg.notch_filter_hz (0/50/60) is the only exposed
 *         knob; Q is fixed here (see kNotchQ) since the config has no
 *         separate bandwidth field.
 */
#include <cmath>
#include <vector>
#include <algorithm>

namespace notch_filter {

    // Q of ~30 gives a narrow stop-band (roughly +-1 Hz around f0 at 50/60 Hz)
    // so it removes powerline hum without visibly touching nearby ECG content.
    inline constexpr double kNotchQ = 30.0;

    struct Biquad { double b0, b1, b2, a1, a2; };  // a0 normalized to 1

    // RBJ audio-EQ-cookbook notch design.
    inline Biquad design_notch(double sampleRateHz, double notchHz, double Q) {
        const double w0 = 2.0 * M_PI * notchHz / sampleRateHz;
        const double alpha = std::sin(w0) / (2.0 * Q);
        const double cosw0 = std::cos(w0);
        const double a0 = 1.0 + alpha;
        return Biquad{
            1.0 / a0, (-2.0 * cosw0) / a0, 1.0 / a0,
            (-2.0 * cosw0) / a0, (1.0 - alpha) / a0
        };
    }

    // One causal pass (Direct Form I).
    inline std::vector<double> apply_biquad(const std::vector<double>& x, const Biquad& c) {
        std::vector<double> y(x.size(), 0.0);
        double x1 = 0, x2 = 0, y1 = 0, y2 = 0;
        for (size_t n = 0; n < x.size(); ++n) {
            const double xn = x[n];
            const double yn = c.b0 * xn + c.b1 * x1 + c.b2 * x2 - c.a1 * y1 - c.a2 * y2;
            y[n] = yn;
            x2 = x1; x1 = xn; y2 = y1; y1 = yn;
        }
        return y;
    }

    // Zero-phase notch: forward pass, reverse the result, forward pass again,
    // reverse back. Two passes deepen the notch (and its skirts) but introduce
    // no net time shift, matching upsample()'s "zero net delay" property.
    // No-op (returns a copy of x) if the filter is disabled or would sit at or
    // past Nyquist for this signal's rate.
    inline std::vector<double> apply(const std::vector<double>& x,
        double sampleRateHz, int notchHz, double Q = kNotchQ) {
        if (notchHz <= 0 || sampleRateHz <= 0.0 ||
            (double)notchHz >= sampleRateHz / 2.0 || x.size() < 4)
            return x;
        const Biquad c = design_notch(sampleRateHz, (double)notchHz, Q);
        std::vector<double> fwd = apply_biquad(x, c);
        std::reverse(fwd.begin(), fwd.end());
        std::vector<double> back = apply_biquad(fwd, c);
        std::reverse(back.begin(), back.end());
        return back;
    }

}  // namespace notch_filter