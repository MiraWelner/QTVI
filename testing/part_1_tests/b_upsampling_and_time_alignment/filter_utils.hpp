#pragma once
// Polyphase FIR resampler extracted verbatim from file_to_bin.cpp so it can be
// unit-tested and shared. Behaviour is identical to the in-place version.
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace filterutils {

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
        processRangeInterior(inPtr, outPtr, safeStartM, safeEndM, bankPtrs, subLen, filterCenter, P, Q);
        return output;
    }

    inline std::vector<double> upsample(const std::vector<double>& input, double sourceRate, double targetRate) {
        if (input.empty()) return {};
        if (sourceRate == targetRate) return input;
        int gcd = greatest_common_divisor((int)targetRate, (int)sourceRate);
        int P = (int)targetRate / gcd;
        int Q = (int)sourceRate / gcd;
        if (P > 1000 || Q > 1000)
            throw std::runtime_error("Resampling ratio too large");
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

} // namespace filterutils