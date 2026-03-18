/**
 * @file   resample.hpp
 * @brief  Polyphase rational resampling to 2000 Hz.
 *         Memory-efficient: no intermediate upsampled array is created.
 *         256->2000 (P=125, Q=16) and 500->2000 (P=4, Q=1).
 */
#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <string>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace resample_detail {

    inline int gcd(int a, int b) {
        a = std::abs(a); b = std::abs(b);
        while (b) { int t = b; b = a % b; a = t; }
        return a;
    }

    /**
     * @brief  Build P polyphase sub-filters, each normalized to sum to 1.0.
     *
     *         The prototype lowpass has cutoff = pi/max(P,Q) in the upsampled
     *         domain, which equals the original Nyquist — preserving all
     *         original bandwidth while preventing aliases.
     *
     * @param  P           Upsample factor
     * @param  Q           Downsample factor
     * @param  halfLobes   Number of sinc lobes per side (controls quality)
     * @return             Vector of P sub-filters
     */
    inline std::vector<std::vector<double>> buildPolyphaseBank(int P, int Q, int halfLobes) {
        int maxPQ = std::max(P, Q);
        int numTaps = 2 * halfLobes * maxPQ + 1;
        double fc = 1.0 / static_cast<double>(maxPQ);  // normalized cutoff
        int M = numTaps - 1;
        double halfM = M / 2.0;

        // Design windowed-sinc prototype
        std::vector<double> h(numTaps);
        for (int n = 0; n < numTaps; ++n) {
            double x = n - halfM;
            double sinc = (std::abs(x) < 1e-12) ? 1.0
                : std::sin(M_PI * fc * x) / (M_PI * x);
            // Blackman window — good stopband, gentle rolloff
            double w = 0.42 - 0.5 * std::cos(2.0 * M_PI * n / M)
                + 0.08 * std::cos(4.0 * M_PI * n / M);
            h[n] = sinc * w;
        }

        // Decompose into P polyphase sub-filters
        int subLen = (numTaps + P - 1) / P;
        std::vector<std::vector<double>> bank(P, std::vector<double>(subLen, 0.0));
        for (int i = 0; i < numTaps; ++i) {
            bank[i % P][i / P] = h[i];
        }

        // Normalize each sub-filter to sum to 1.0
        // This guarantees unity gain and prevents amplitude spikes
        for (int p = 0; p < P; ++p) {
            double s = 0.0;
            for (int k = 0; k < subLen; ++k) s += bank[p][k];
            if (std::abs(s) > 1e-15) {
                for (int k = 0; k < subLen; ++k) bank[p][k] /= s;
            }
        }

        return bank;
    }

    /**
     * @brief  Polyphase rational resampling by P/Q.
     *         Only computes output samples — never builds the upsampled signal.
     *         Memory usage = input + output + filter bank (tiny).
     */
    inline std::vector<double> polyphaseResample(const std::vector<double>& input, int P, int Q) {
        if (input.empty()) return {};
        if (P == 1 && Q == 1) return input;

        // Use more lobes for larger ratios to avoid artifacts
        int halfLobes = std::max(16, std::max(P, Q) / 2);
        auto bank = buildPolyphaseBank(P, Q, halfLobes);
        int subLen = static_cast<int>(bank[0].size());

        long long inLen = static_cast<long long>(input.size());
        long long outLen = static_cast<long long>(
            std::ceil(static_cast<double>(inLen) * P / Q));
        std::vector<double> output(outLen);

        // The filter center in the upsampled domain is at tap (numTaps-1)/2.
        // In input-sample units that's halfLobes * max(P,Q) / P.
        int maxPQ = std::max(P, Q);
        int filterCenter = halfLobes * maxPQ / P;

        for (long long m = 0; m < outLen; ++m) {
            long long upsampledIdx = m * Q;
            int phase = static_cast<int>(upsampledIdx % P);
            long long baseInput = upsampledIdx / P;

            double sum = 0.0;
            const auto& sub = bank[phase];
            for (int k = 0; k < subLen; ++k) {
                long long inIdx = baseInput - k + filterCenter;
                if (inIdx >= 0 && inIdx < inLen) {
                    sum += sub[k] * input[static_cast<size_t>(inIdx)];
                }
            }
            output[m] = sum;
        }

        return output;
    }

} // namespace resample_detail

/**
 * @brief  Resample a signal to 2000 Hz using polyphase filtering.
 *         Supports 256 Hz and 500 Hz source rates.
 *         Memory efficient — only input + output arrays are allocated.
 *
 * @param  input       Input samples
 * @param  sourceRate  Original sampling rate (256.0 or 500.0)
 * @return             Resampled signal at 2000 Hz
 */
inline std::vector<double> upsample(const std::vector<double>& input, double sourceRate) {
    if (input.empty()) return {};

    static constexpr double TARGET_RATE = 2000.0;
    if (std::abs(sourceRate - TARGET_RATE) < 0.01) return input;

    int rateNum = static_cast<int>(std::round(TARGET_RATE));
    int rateDen = static_cast<int>(std::round(sourceRate));
    int g = resample_detail::gcd(rateNum, rateDen);
    int P = rateNum / g;
    int Q = rateDen / g;

    if (P > 1000 || Q > 1000) {
        throw std::runtime_error("Resampling ratio " + std::to_string(P) + "/" +
            std::to_string(Q) + " too large — source rate: " +
            std::to_string(sourceRate));
    }

    return resample_detail::polyphaseResample(input, P, Q);
}