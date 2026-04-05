
#pragma once

#include "common.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

/**
 * @file dicrotic_notch.hpp
 * @brief Detect the dicrotic notch in a PPG beat waveform.
 *
 * Port of dumbDicrotic.m.
 * Uses a shear-transform approach to locate the notch between the
 * systolic peak and the end of the beat.
 */

namespace ppg {

    namespace detail {

        /**
         * @brief Apply a shear transform to a segment.
         * @param x_indices  Sample indices of the segment.
         * @param y_values   Corresponding amplitude values.
         * @return Shear-transformed signal.
         */
        inline std::vector<double> shear_transform(
            const std::vector<int>& x_indices,
            const std::vector<double>& y_values)
        {
            const int n = static_cast<int>(y_values.size());
            if (n < 2) return y_values;

            double x0 = x_indices.front(), xn = x_indices.back();
            double y0 = y_values.front(), yn = y_values.back();
            double slope = (yn - y0) / (xn - x0);

            std::vector<double> result(n);
            for (int i = 0; i < n; ++i) {
                double baseline = y0 + slope * (x_indices[i] - x0);
                result[i] = y_values[i] - baseline;
            }
            return result;
        }

        /**
         * @brief Check if max orthogonal distance exceeds threshold.
         * @param norm_time  Normalised time axis [0,1].
         * @param norm_line  Normalised baseline.
         * @param norm_press Normalised pressure.
         * @param threshold  Distance threshold.
         * @return True if max distance exceeds threshold.
         */
        inline bool orthogonal_dist_thresh(
            const std::vector<double>& norm_time,
            const std::vector<double>& norm_line,
            const std::vector<double>& norm_press,
            double threshold)
        {
            for (size_t i = 0; i < norm_time.size(); ++i) {
                double d = std::abs(norm_press[i] - norm_line[i]);
                if (d > threshold) return true;
            }
            return false;
        }

    } // namespace detail

    /**
     * @brief Detect the dicrotic notch index in a single PPG beat.
     * @param beat_raw    PPG beat waveform (valley to valley).
     * @param sp_ratio    Optional systolic-peak-to-end ratio estimate.
     *                    Pass NaN or omit for automatic estimation.
     * @param smooth_win  Smoothing window width in samples (default 117 ≈ 59ms at 2000 Hz).
     * @return 0-based index of the dicrotic notch within the beat, or -1 if not found.
     */
    inline int dicrotic_notch(
        const std::vector<double>& beat_raw,
        double sp_ratio = kNaN,
        int smooth_win = 117)
    {
        if (beat_raw.size() < 4) return -1;

        // Smooth (window width should be ~59ms; at 2000 Hz that's ~117 samples)
        std::vector<double> beat = fast_smooth(beat_raw, smooth_win);
        if (beat.empty()) return -1;

        const int n = static_cast<int>(beat.size());

        // Find systolic peak
        auto it_max = std::max_element(beat.begin(), beat.end());
        int pmax = static_cast<int>(it_max - beat.begin());
        int pmin = n - 1;

        if (pmax == pmin || pmin - pmax < 2) return -1;

        // Estimate regions
        int p_min_dpdt_region, init_EP;
        if (!std::isnan(sp_ratio)) {
            int endlen = n - 1 - pmax;
            int notch_estimate = pmax + static_cast<int>(std::round(endlen * sp_ratio));
            p_min_dpdt_region = pmax + static_cast<int>(std::round((notch_estimate - pmax) * (2.0 / 3.0)));
            if (p_min_dpdt_region < pmax + 2)
                p_min_dpdt_region = pmax + static_cast<int>(std::round((pmin - pmax) / 3.0));
            init_EP = notch_estimate + static_cast<int>(std::round((n - 1 - notch_estimate) * 0.5));
            if (init_EP <= p_min_dpdt_region)
                init_EP = pmax + static_cast<int>(std::round((pmin - pmax) * 3.0 / 4.0));
        }
        else {
            p_min_dpdt_region = pmax + static_cast<int>(std::round((pmin - pmax) / 3.0));
            init_EP = pmax + static_cast<int>(std::round((pmin - pmax) * 3.0 / 4.0));
        }

        p_min_dpdt_region = std::min(p_min_dpdt_region, n - 1);
        init_EP = std::min(init_EP, n - 1);

        // Min derivative in [pmax, p_min_dpdt_region]
        std::vector<double> slice(beat.begin() + pmax, beat.begin() + p_min_dpdt_region + 1);
        auto d_slice = diff(slice);
        if (d_slice.empty()) return -1;
        int p_min_dpdt = static_cast<int>(std::min_element(d_slice.begin(), d_slice.end()) - d_slice.begin());
        p_min_dpdt += pmax;

        // Half-amplitude threshold
        double p_half = beat[pmax] - (beat[pmax] - beat[p_min_dpdt]) / 2.0;

        // Find candidate start point (SP) closest to p_half
        std::vector<std::pair<double, int>> candidates;
        for (int k = pmax; k <= p_min_dpdt_region && k < n; ++k)
            candidates.push_back({ std::abs(beat[k] - p_half), k });
        std::sort(candidates.begin(), candidates.end());

        int SP = pmax;
        for (auto& [dist, idx] : candidates) {
            // Check shear-transform criterion
            if (idx >= init_EP || idx >= n) continue;
            std::vector<int>    xi;
            std::vector<double> yi;
            for (int k = idx; k <= init_EP && k < n; ++k) {
                xi.push_back(k);
                yi.push_back(beat[k]);
            }
            auto transform = detail::shear_transform(xi, yi);
            int below = 0;
            for (size_t m = 0; m < transform.size(); ++m)
                if (transform[m] < yi[m]) ++below;
            if (static_cast<double>(below) / transform.size() < 0.5) {
                SP = idx;
                break;
            }
        }

        // Shrink EP
        int EP = init_EP;
        while (EP > SP + 1) {
            std::vector<double> seg(beat.begin() + SP, beat.begin() + EP + 1);
            int seg_len = static_cast<int>(seg.size());
            if (seg_len < 2) break;

            double slope_val = (seg.back() - seg.front()) / (seg_len - 1);

            // Normalise
            auto norm = [](const std::vector<double>& v) {
                double mn = *std::min_element(v.begin(), v.end());
                double mx = *std::max_element(v.begin(), v.end());
                double rng = mx - mn;
                std::vector<double> r(v.size());
                for (size_t i = 0; i < v.size(); ++i)
                    r[i] = (rng > 1e-15) ? (v[i] - mn) / rng : 0.0;
                return r;
                };

            std::vector<double> shearline(seg_len);
            double m = seg.back() / seg.front(); // simplified slope
            for (int k = 0; k < seg_len; ++k)
                shearline[k] = slope_val * k;

            auto norm_line = norm(shearline);
            auto norm_press = norm(seg);
            auto norm_time = linspace(0.0, 1.0, seg_len);

            if (detail::orthogonal_dist_thresh(norm_time, norm_line, norm_press, 0.3))
                --EP;
            else
                break;
        }

        if (SP >= EP || SP >= n || EP >= n) return -1;

        // Shear transform on final [SP, EP] to find minimum
        std::vector<int>    xi_final;
        std::vector<double> yi_final;
        for (int k = SP; k <= EP && k < n; ++k) {
            xi_final.push_back(k);
            yi_final.push_back(beat[k]);
        }
        auto transform = detail::shear_transform(xi_final, yi_final);
        int min_shear_local = static_cast<int>(
            std::min_element(transform.begin(), transform.end()) - transform.begin());
        int min_shear = SP + min_shear_local;

        // Find start of diastolic relaxation (local max after min_shear)
        int start_diastolic = min_shear;
        double dmax = beat[min_shear];
        for (int k = min_shear; k <= pmin && k < n; ++k) {
            if (beat[k] > dmax) { dmax = beat[k]; start_diastolic = k; }
        }

        // Dicrotic notch = min between min_shear and start_diastolic
        if (min_shear >= start_diastolic) return -1;
        int notch = min_shear;
        double nmin = beat[min_shear];
        for (int k = min_shear; k <= start_diastolic && k < n; ++k) {
            if (beat[k] < nmin) { nmin = beat[k]; notch = k; }
        }

        return notch;
    }

} // namespace ppg