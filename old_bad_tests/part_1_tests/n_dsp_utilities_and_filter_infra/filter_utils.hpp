#pragma once
/**
 * @file  filter_utils.hpp
 * @brief Shared signal-processing helpers: a polyphase FIR resampler, a
 *        zero-phase (forward-backward) IIR lowpass, and linear detrending.
 *
 *        This header is the single source of truth for the resampler --
 *        file_to_bin.cpp includes it rather than carrying its own copy, so
 *        the unit tests exercise exactly the code that converts files.
 *
 * @note  polyphase_resample() spawns std::thread for the interior loop, so any
 *        target that includes this header and calls upsample() must link with
 *        -pthread on GCC/Clang.
 */
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
 // MSVC omits M_SQRT2 unless _USE_MATH_DEFINES is set, same as M_PI above.
#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

namespace filterutils {

    // ========================================================================
    // Polyphase FIR resampler
    //
    // Builds a filter bank of small filters (one per fractional phase between
    // input samples), then for each output sample dot-products the right
    // sub-filter against a window of inputs. Output is split into a leading
    // boundary (filter hangs off the left edge), a multi-threaded interior (no
    // bounds checks), and a trailing boundary.
    // ========================================================================

    inline int greatest_common_divisor(int a, int b) {
        a = std::abs(a); b = std::abs(b);
        while (b) { int t = b; b = a % b; a = t; }
        return a;
    }

    inline std::vector<std::vector<double>> buildPolyphaseBank(int P, int Q, int halfLobes) {
        int maxPQ = std::max(P, Q);
        int numTaps = 2 * halfLobes * maxPQ + 1;
        double fc = 1.0 / static_cast<double>(maxPQ);
        int M = numTaps - 1;
        double halfM = M / 2.0;
        std::vector<double> h(numTaps);
        for (int n = 0; n < numTaps; ++n) {
            double x = n - halfM;
            double sinc = (std::abs(x) < 1e-12) ? 1.0 : std::sin(M_PI * fc * x) / (M_PI * x);
            double w = 0.42 - 0.5 * std::cos(2.0 * M_PI * n / M) + 0.08 * std::cos(4.0 * M_PI * n / M);
            h[n] = sinc * w;
        }
        int subLen = (numTaps + P - 1) / P;
        std::vector<std::vector<double>> bank(P, std::vector<double>(subLen, 0.0));
        for (int i = 0; i < numTaps; ++i) bank[i % P][i / P] = h[i];
        for (int p = 0; p < P; ++p) {
            double s = 0.0;
            for (int k = 0; k < subLen; ++k) s += bank[p][k];
            if (std::abs(s) > 1e-15) for (int k = 0; k < subLen; ++k) bank[p][k] /= s;
        }
        return bank;
    }

    inline void processRangeBoundary(const double* inPtr, long long inLen, double* outPtr,
        long long mStart, long long mEnd, const std::vector<const double*>& bankPtrs,
        int subLen, int filterCenter, int P, int Q) {
        for (long long m = mStart; m < mEnd; ++m) {
            long long upsampledIdx = m * Q;
            int phase = static_cast<int>(upsampledIdx % P);
            long long baseInput = upsampledIdx / P;
            const double* sub = bankPtrs[phase];
            double sum = 0.0;
            for (int k = 0; k < subLen; ++k) {
                long long inIdx = baseInput - k + filterCenter;
                if (inIdx >= 0 && inIdx < inLen) sum += sub[k] * inPtr[inIdx];
            }
            outPtr[m] = sum;
        }
    }

    inline void processRangeInterior(const double* inPtr, double* outPtr, long long mStart,
        long long mEnd, const std::vector<const double*>& bankPtrs, int subLen,
        int filterCenter, int P, int Q) {
        for (long long m = mStart; m < mEnd; ++m) {
            long long upsampledIdx = m * Q;
            int phase = static_cast<int>(upsampledIdx % P);
            long long baseInput = upsampledIdx / P;
            const double* sub = bankPtrs[phase];
            const double* inBase = inPtr + baseInput + filterCenter;
            double sum = 0.0;
            for (int k = 0; k < subLen; ++k) sum += sub[k] * inBase[-k];
            outPtr[m] = sum;
        }
    }

    inline std::vector<double> polyphase_resample(const std::vector<double>& input, int P, int Q) {
        if (input.empty()) return {};
        if (P == 1 && Q == 1) return input;
        int halfLobes = std::max(16, std::max(P, Q) / 2);
        auto bank = buildPolyphaseBank(P, Q, halfLobes);
        int subLen = static_cast<int>(bank[0].size());
        long long inLen = static_cast<long long>(input.size());
        long long outLen = static_cast<long long>(std::ceil(static_cast<double>(inLen) * P / Q));
        std::vector<double> output(outLen);
        int maxPQ = std::max(P, Q);
        int filterCenter = halfLobes * maxPQ / P;
        const double* inPtr = input.data();
        double* outPtr = output.data();
        std::vector<const double*> bankPtrs(P);
        for (int p = 0; p < P; ++p) bankPtrs[p] = bank[p].data();
        long long safeStartM = 0, safeEndM = 0;
        for (long long m = 0; m < outLen; ++m) {
            long long baseInput = (m * Q) / P;
            if (baseInput - subLen + 1 + filterCenter >= 0 && baseInput + filterCenter < inLen) { safeStartM = m; break; }
        }
        for (long long m = outLen - 1; m >= safeStartM; --m) {
            long long baseInput = (m * Q) / P;
            if (baseInput - subLen + 1 + filterCenter >= 0 && baseInput + filterCenter < inLen) { safeEndM = m + 1; break; }
        }
        processRangeBoundary(inPtr, inLen, outPtr, 0, safeStartM, bankPtrs, subLen, filterCenter, P, Q);
        processRangeBoundary(inPtr, inLen, outPtr, safeEndM, outLen, bankPtrs, subLen, filterCenter, P, Q);
        long long interiorLen = safeEndM - safeStartM;
        if (interiorLen <= 0) return output;

        // Each output index m is written exactly once and reads only from the
        // read-only input and bank, so splitting the interior into chunks
        // cannot change results: output is bit-for-bit identical to the
        // single-threaded path.
        unsigned int numThreads = std::thread::hardware_concurrency();
        if (numThreads == 0) numThreads = 4;
        if (interiorLen < 100000) numThreads = 1;
        numThreads = std::min(numThreads, static_cast<unsigned int>(interiorLen));

        if (numThreads == 1) {
            processRangeInterior(inPtr, outPtr, safeStartM, safeEndM,
                bankPtrs, subLen, filterCenter, P, Q);
        }
        else {
            std::vector<std::thread> threads;
            threads.reserve(numThreads);
            long long chunkSize = interiorLen / numThreads;
            for (unsigned int t = 0; t < numThreads; ++t) {
                long long start = safeStartM + t * chunkSize;
                long long end = (t == numThreads - 1) ? safeEndM : start + chunkSize;
                threads.emplace_back(processRangeInterior, inPtr, outPtr, start, end,
                    std::cref(bankPtrs), subLen, filterCenter, P, Q);
            }
            for (auto& th : threads) th.join();
        }
        return output;
    }

    inline std::vector<double> upsample(const std::vector<double>& input, double sourceRate, double targetRate) {
        if (input.empty()) return {};
        if (sourceRate == targetRate) return input;
        int gcd = greatest_common_divisor((int)targetRate, (int)sourceRate);
        int P = (int)targetRate / gcd;
        int Q = (int)sourceRate / gcd;
        if (P > 1000 || Q > 1000) {
            throw std::runtime_error(
                "Resampling ratio " + std::to_string(P) + "/" +
                std::to_string(Q) + " too large - source rate: " +
                std::to_string(sourceRate));
        }
        return polyphase_resample(input, P, Q);
    }

    // Group delay of the interior path, in OUTPUT samples: the interior taps read
    // inBase[-k] = inPtr[baseInput+filterCenter-k], so the effective input-side
    // center offset is filterCenter input samples -> filterCenter*P/Q output samples.
    inline double group_delay_out_samples(double sourceRate, double targetRate) {
        int gcd = greatest_common_divisor((int)targetRate, (int)sourceRate);
        int P = (int)targetRate / gcd, Q = (int)sourceRate / gcd;
        int halfLobes = std::max(16, std::max(P, Q) / 2);
        int maxPQ = std::max(P, Q);
        int filterCenter = halfLobes * maxPQ / P;
        return (double)filterCenter * (double)P / (double)Q;
    }

    // ========================================================================
    // Zero-phase IIR lowpass filtering and linear detrending
    //
    // Unlike the resampler above, these introduce no delay: filtfilt's
    // forward-backward passes cancel the phase response, and detrend is not a
    // filter at all. Do NOT apply group_delay_out_samples() to their output --
    // that correction exists only for the causal FIR resampler path.
    // ========================================================================

    /// Second-order section: y = b[0]*x + b[1]*x[-1] + b[2]*x[-2]
    ///                         - a[1]*y[-1] - a[2]*y[-2]   (a[0] == 1)
    struct Biquad { double b[3], a[3]; };

    /**
     * @brief 2nd-order Butterworth lowpass via the bilinear transform.
     * @param fc  Cutoff in Hz. Must satisfy 0 < fc < fs/2; the prewarping
     *            tan(pi*fc/fs) diverges at Nyquist and the coefficients are
     *            meaningless outside that range (not checked here).
     * @param fs  Sample rate in Hz.
     */
    inline Biquad butterLP(double fc, double fs) {
        double w = std::tan(M_PI * fc / fs);
        double w2 = w * w, norm = 1.0 / (1.0 + M_SQRT2 * w + w2);
        Biquad bq;
        bq.b[0] = w2 * norm; bq.b[1] = 2.0 * bq.b[0]; bq.b[2] = bq.b[0];
        bq.a[0] = 1.0; bq.a[1] = 2.0 * (w2 - 1.0) * norm; bq.a[2] = (1.0 - M_SQRT2 * w + w2) * norm;
        return bq;
    }

    /**
     * @brief Zero-phase filtering (filtfilt): run the section forward, then
     *        backward, so the phase responses cancel. This doubles the
     *        effective order -- a butterLP section becomes 4th-order,
     *        24 dB/octave, with a squared magnitude response.
     *
     * @note Both passes start from zero state, so there is a startup transient
     *       at each end whose length scales with fs/fc. SciPy suppresses this
     *       with reflected padding and a steady-state initial condition; this
     *       version does not.
     */
    inline std::vector<double> filtfilt(const Biquad& bq, const std::vector<double>& x) {
        auto fwd = [&](const std::vector<double>& in) {
            std::vector<double> out(in.size());
            double z1 = 0, z2 = 0;
            for (size_t i = 0; i < in.size(); ++i) {
                double w = in[i] - bq.a[1] * z1 - bq.a[2] * z2;
                out[i] = bq.b[0] * w + bq.b[1] * z1 + bq.b[2] * z2;
                z2 = z1; z1 = w;
            }
            return out;
            };
        auto y = fwd(x);
        std::reverse(y.begin(), y.end());
        y = fwd(y);
        std::reverse(y.begin(), y.end());
        return y;
    }

    /**
     * @brief Linear detrend: least-squares fit of a line against sample index,
     *        subtracted from the signal. Returns {} for empty input.
     */
    inline std::vector<double> detrend(const std::vector<double>& x) {
        if (x.empty()) return {};
        long long n = static_cast<long long>(x.size());
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        for (long long i = 0; i < n; ++i) {
            double di = static_cast<double>(i);
            sx += di; sy += x[i]; sxx += di * di; sxy += di * x[i];
        }
        double dn = static_cast<double>(n);
        double den = dn * sxx - sx * sx;
        double m = (den != 0) ? (dn * sxy - sx * sy) / den : 0;
        double c = (sy - m * sx) / dn;
        std::vector<double> out(x.size());
        for (long long i = 0; i < n; ++i) out[i] = x[i] - (m * static_cast<double>(i) + c);
        return out;
    }

} // namespace filterutils