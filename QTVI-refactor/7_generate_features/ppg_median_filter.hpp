#pragma once

#include <vector>
#include <algorithm>
#include <cmath>

/**
 * @file ppg_median_filter.hpp
 * @brief Median-based baseline wander removal for PPG signals.
 *
 * Port of PPGmedianfilter.m (Khaustov, Incart, 2008; modified Da Poian 2017).
 */

namespace ppg {

    /**
     * @brief Remove baseline wander via windowed median-like filter.
     * @param x      Input PPG signal.
     * @param freq   Sampling frequency of the input signal.
     * @param freqd  Desired decimation frequency (default 125).
     * @return Baseline-corrected signal (same length as input).
     */
    inline std::vector<double> ppg_median_filter(
        const std::vector<double>& x,
        double freq,
        double freqd = 125.0)
    {
        const int n = static_cast<int>(x.size());
        if (n == 0) return {};

        int factor = std::max(1, static_cast<int>(std::floor(freq / freqd)));

        // Decimate
        std::vector<double> y;
        for (int i = 0; i < n; i += factor)
            y.push_back(x[i]);

        const int ny = static_cast<int>(y.size());
        int halfwin = static_cast<int>(std::floor(1.3 * freqd / 2.0));

        // Windowed trimmed mean (37.5% – 62.5% of sorted window)
        std::vector<double> m(ny);
        for (int i = 0; i < ny; ++i) {
            int lo = std::max(0, i - halfwin);
            int hi = std::min(ny - 1, i + halfwin);

            std::vector<double> ss(y.begin() + lo, y.begin() + hi + 1);
            std::sort(ss.begin(), ss.end());

            int slen = static_cast<int>(ss.size());
            int i_lo = static_cast<int>(std::floor(0.375 * slen));
            int i_hi = static_cast<int>(std::floor(0.625 * slen));
            i_lo = std::max(0, i_lo);
            i_hi = std::min(slen - 1, i_hi);

            double sum = 0.0;
            int cnt = 0;
            for (int k = i_lo; k <= i_hi; ++k) { sum += ss[k]; ++cnt; }
            m[i] = (cnt > 0) ? sum / cnt : 0.0;
        }

        // Upsample baseline estimate back to original rate
        std::vector<double> res(n);
        for (int i = 0; i < n; ++i) {
            int idx = std::min(static_cast<int>((i) / factor), ny - 1);
            res[i] = m[idx];
        }

        // Subtract baseline
        std::vector<double> xfilt(n);
        for (int i = 0; i < n; ++i)
            xfilt[i] = x[i] - res[i];

        return xfilt;
    }

} // namespace ppg