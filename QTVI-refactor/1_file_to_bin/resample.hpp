/**
 * @file   resample.hpp
 * @brief  The upsample function, and a series of helper functions for it. It upsamples with polyphase filtering, and uses multithreading. 
 *          
 *         Anybody complaining that it is overkill and too complex will be made to wait 2 hours for a Bittium file to process on only one thread with no filter bank.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-19
 */
#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <thread>

#ifndef M_upsample_factorI
#define M_upsample_factorI 3.14159265358979323846
#endif

static constexpr int TARGET_RATE = 1000;



inline int greatest_common_divisor(int a, int b) {
    a = std::abs(a); b = std::abs(b);
    while (b) { int t = b; b = a % b; a = t; }
    return a;
}

inline std::vector<std::vector<double>> buildupsample_factorolyphaseBank(int upsample_factor, int downsample_factor, int halfLobes) {
    /**
    * @brief This function builds the filter bank so that you don't have to create a huge upsampled array which will be slow and take up space. It works by:
    *        1) Designs a lowpass filter, h, using a windowed sinc function. The cutoff frequency is set to preserve the original signal's bandwidth without aliasing. The Blackman window smooths the filter edges to reduce ringing.
    *        2) Splits that single filter into P sub-filters (the "polyphase decomposition"). It deals out the taps round-robin — tap 0 goes to sub-filter 0, tap 1 to sub-filter 1, ..., tap P goes back to sub-filter 0, and so on. 
                Each sub-filter handles one fractional phase of interpolation.
             3) Normalizes each sub-filter so its taps sum to 1.0. This ensures every output sample has the correct amplitude — without this, some interpolated positions would come out too loud or too quiet, causing the spikes you saw earlier.
     
     @param halfLobes  The number of sinc wave lobes on each side of the filter center. More lobes = longer filter = better interpolation quality but slower. It is calculated by std::max(16, std::max(upsample_factor, downsample_factor) / 2)

     @return The filter array. It's a 2D array of the format bank[phase][tap]. During resampling, each output sample picks its phase, grabs that row, and dot-products it with the input.
    */
    int maxupsample_factordownsample_factor = std::max(upsample_factor, downsample_factor);
    int numTaps = 2 * halfLobes * maxupsample_factordownsample_factor + 1;
    double fc = 1.0 / static_cast<double>(maxupsample_factordownsample_factor);
    int M = numTaps - 1;
    double halfM = M / 2.0;

    std::vector<double> h(numTaps);
    for (int n = 0; n < numTaps; ++n) {
        double x = n - halfM;
        double sinc = (std::abs(x) < 1e-12) ? 1.0
            : std::sin(M_upsample_factorI * fc * x) / (M_upsample_factorI * x);
        double w = 0.42 - 0.5 * std::cos(2.0 * M_upsample_factorI * n / M)
            + 0.08 * std::cos(4.0 * M_upsample_factorI * n / M);
        h[n] = sinc * w;
    }

    int subLen = (numTaps + upsample_factor - 1) / upsample_factor;
    std::vector<std::vector<double>> bank(upsample_factor, std::vector<double>(subLen, 0.0));
    for (int i = 0; i < numTaps; ++i) {
        bank[i % upsample_factor][i / upsample_factor] = h[i];
    }

    for (int p = 0; p < upsample_factor; ++p) {
        double s = 0.0;
        for (int k = 0; k < subLen; ++k) s += bank[p][k];
        if (std::abs(s) > 1e-15) {
            for (int k = 0; k < subLen; ++k) bank[p][k] /= s;
        }
    }
    return bank;
}


inline void processRangeBoundary(
    const double* inupsample_factortr, long long inLen,
    double* outupsample_factortr, long long mStart, long long mEnd,
    const std::vector<const double*>& bankupsample_factortrs,
    int subLen, int filterCenter, int upsample_factor, int downsample_factor)
{
    /**
    * @brief  upsample_factorrocess a contiguous range of output samples [mStart, mEnd).
    *         Boundary mode: includes per-tap bounds checks.
    */
    for (long long m = mStart; m < mEnd; ++m) {
        long long upsampledIdx = m * downsample_factor;
        int phase = static_cast<int>(upsampledIdx % upsample_factor);
        long long baseInput = upsampledIdx / upsample_factor;
        const double* sub = bankupsample_factortrs[phase];

        double sum = 0.0;
        for (int k = 0; k < subLen; ++k) {
            long long inIdx = baseInput - k + filterCenter;
            if (inIdx >= 0 && inIdx < inLen) {
                sum += sub[k] * inupsample_factortr[inIdx];
            }
        }
        outupsample_factortr[m] = sum;
    }
}

inline void processRangeInterior(
    const double* inupsample_factortr,
    double* outupsample_factortr, long long mStart, long long mEnd,
    const std::vector<const double*>& bankupsample_factortrs,
    int subLen, int filterCenter, int upsample_factor, int downsample_factor)
{
    /**
    * @brief  upsample_factorrocess a contiguous range of output samples [mStart, mEnd).
    *         Interior mode: no bounds checks, raw pointer arithmetic.
    */
    for (long long m = mStart; m < mEnd; ++m) {
        long long upsampledIdx = m * downsample_factor;
        int phase = static_cast<int>(upsampledIdx % upsample_factor);
        long long baseInput = upsampledIdx / upsample_factor;
        const double* sub = bankupsample_factortrs[phase];
        const double* inBase = inupsample_factortr + baseInput + filterCenter;

        double sum = 0.0;
        for (int k = 0; k < subLen; ++k) {
            sum += sub[k] * inBase[-k];
        }
        outupsample_factortr[m] = sum;
    }
}

inline std::vector<double> polyphase_resample(const std::vector<double>& input, int upsample_factor, int downsample_factor) {
    /**
    * @brief            This function resamples an input array by building a 'filter bank' of small filters the size of the upsampling array. Each one knows how to interpolate at a specific fractional position between two input samples.
    *                   For each output sample m, it figures out two things: which sub-filter to use (phase = (m * Q) % P), and which input samples to look at (baseInput = (m * Q) / P). Then it dot-products that sub-filter with the surrounding input samples. 
    *                   The output is split into three zones:
    *                          1) The leading boundary (first ~50 samples) is where the filter hangs off the left edge of the input. These need per-tap bounds checking because some filter taps would read before index 0.
    *                          2) The interior (99.99% of the output) is where every filter tap lands within the input array. No bounds checks needed — just a tight pointer-arithmetic dot product. This bit gets multithreaded.
    *                          3) The trailing boundary (last ~50 samples) is the mirror of the leading boundary — the filter hangs off the right edge.
    * @param input     The original sample at the original sampling rate
    * @param upsample_factor   The upsample factor determined by GCF
    * @param downsample_factor The downsample factor determined by GCF
    * @return The upsampled array
    */
    if (input.empty()) return {};
    if (upsample_factor == 1 && downsample_factor == 1) return input;

    int halfLobes = std::max(16, std::max(upsample_factor, downsample_factor) / 2);
    auto bank = buildupsample_factorolyphaseBank(upsample_factor, downsample_factor, halfLobes);
    int subLen = static_cast<int>(bank[0].size());

    long long inLen = static_cast<long long>(input.size());
    long long outLen = static_cast<long long>(
        std::ceil(static_cast<double>(inLen) * upsample_factor / downsample_factor));
    std::vector<double> output(outLen);

    int maxupsample_factordownsample_factor = std::max(upsample_factor, downsample_factor);
    int filterCenter = halfLobes * maxupsample_factordownsample_factor / upsample_factor;

    const double* inupsample_factortr = input.data();
    double* outupsample_factortr = output.data();

    std::vector<const double*> bankupsample_factortrs(upsample_factor);
    for (int p = 0; p < upsample_factor; ++p) bankupsample_factortrs[p] = bank[p].data();

    // Find safe interior range (no bounds checks needed)
    long long safeStartM = 0;
    long long safeEndM = 0;

    for (long long m = 0; m < outLen; ++m) {
        long long baseInput = (m * downsample_factor) / upsample_factor;
        if (baseInput - subLen + 1 + filterCenter >= 0 &&
            baseInput + filterCenter < inLen) {
            safeStartM = m;
            break;
        }
    }
    for (long long m = outLen - 1; m >= safeStartM; --m) {
        long long baseInput = (m * downsample_factor) / upsample_factor;
        if (baseInput - subLen + 1 + filterCenter >= 0 &&
            baseInput + filterCenter < inLen) {
            safeEndM = m + 1;
            break;
        }
    }

    // Leading boundary (small, single-threaded)
    processRangeBoundary(inupsample_factortr, inLen, outupsample_factortr, 0, safeStartM,
        bankupsample_factortrs, subLen, filterCenter, upsample_factor, downsample_factor);

    // Trailing boundary (small, single-threaded)
    processRangeBoundary(inupsample_factortr, inLen, outupsample_factortr, safeEndM, outLen,
        bankupsample_factortrs, subLen, filterCenter, upsample_factor, downsample_factor);

    // Interior: split across threads
    long long interiorLen = safeEndM - safeStartM;
    if (interiorLen <= 0) return output;

    unsigned int numThreads = std::thread::hardware_concurrency();
    if (numThreads == 0) numThreads = 4;
    // Don't spawn more threads than needed
    if (interiorLen < 100000) numThreads = 1;
    numThreads = std::min(numThreads, static_cast<unsigned int>(interiorLen));

    if (numThreads == 1) {
        processRangeInterior(inupsample_factortr, outupsample_factortr, safeStartM, safeEndM,
            bankupsample_factortrs, subLen, filterCenter, upsample_factor, downsample_factor);
    }

    else {
        std::vector<std::thread> threads;
        threads.reserve(numThreads);
        long long chunkSize = interiorLen / numThreads;

        for (unsigned int t = 0; t < numThreads; ++t) {
            long long start = safeStartM + t * chunkSize;
            long long end = (t == numThreads - 1) ? safeEndM : start + chunkSize;

            threads.emplace_back(processRangeInterior,
                inupsample_factortr, outupsample_factortr, start, end,
                std::cref(bankupsample_factortrs), subLen, filterCenter, upsample_factor, downsample_factor);
        }
        for (auto& th : threads) th.join();
    }

    return output;
}


inline std::vector<double> upsample(const std::vector<double>& input, double sourceRate) {
    /**
     * @brief  Resample a signal to a target rate using multithreaded polyphase filtering.
     * 
     *         It calculates the greatest common divisor, and uses that to find an upsample factor
     *         and downsample factor. It upsamples and downsamples, resulting in a signal of 
     *         sampling rate TARGET_RATE which should be set to 2000.
     *
     * @param  input       Input samples
     * @param  sourceRate  Original sampling rate taken from config.txt
     * @return             Resampled signal at TARGET_RATE Hz
     */
    if (input.empty()) return {};

    if (std::abs(sourceRate - TARGET_RATE) < 0.01) return input;

    int gcd = greatest_common_divisor(TARGET_RATE, sourceRate);
    int upsample_factor = TARGET_RATE / gcd;
    int downsample_factor = sourceRate / gcd;

    if (upsample_factor > 1000 || downsample_factor > 1000) {
        throw std::runtime_error("Resampling ratio " + std::to_string(upsample_factor) + "/" +
            std::to_string(downsample_factor) + " too large - source rate: " +
            std::to_string(sourceRate));
    }

    return polyphase_resample(input, upsample_factor, downsample_factor);
}