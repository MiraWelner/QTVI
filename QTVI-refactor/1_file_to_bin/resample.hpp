/**
 * @file   resample.hpp
 * @brief  Linear-interpolation upsample to TARGET_RATE.
 *         For each output slot, the two nearest input samples are blended
 *         proportionally to the output slot's position between them. The
 *         result is a smooth, continuous waveform rather than a stepped one.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 */
#pragma once

#include <vector>
#include <cmath>
#include <cstddef>

static constexpr int TARGET_RATE = 1000;

inline std::vector<double> upsample(const std::vector<double>& input, double sourceRate) {
    /**
     * @brief  Resample to TARGET_RATE using linear interpolation between
     *         adjacent input samples.
     *
     *   time_of_output_m = m / TARGET_RATE     (seconds)
     *   fractional_input_idx = m * (sourceRate / TARGET_RATE)
     *   i0 = floor(frac), i1 = i0 + 1
     *   output[m] = input[i0] * (1 - f) + input[i1] * f   where f = frac - i0
     */
    if (input.empty()) return {};
    if (sourceRate <= 0.0) return input;
    if (std::abs(sourceRate - TARGET_RATE) < 0.01) return input;

    const size_t inLen = input.size();

    double outLenD = std::ceil(static_cast<double>(inLen) * TARGET_RATE / sourceRate);
    if (outLenD <= 0.0 || outLenD > 1e10) return input;  // sanity bail
    const size_t outLen = static_cast<size_t>(outLenD);

    std::vector<double> output(outLen);

    const double step = sourceRate / static_cast<double>(TARGET_RATE);

    for (size_t m = 0; m < outLen; ++m) {
        const double frac = static_cast<double>(m) * step;
        const size_t i0 = static_cast<size_t>(frac);

        if (i0 >= inLen - 1) {
            // At or past the last input sample: hold the final value
            // (nothing to interpolate towards).
            output[m] = input[inLen - 1];
            continue;
        }

        const double f = frac - static_cast<double>(i0);
        output[m] = input[i0] * (1.0 - f) + input[i0 + 1] * f;
    }
    return output;
}