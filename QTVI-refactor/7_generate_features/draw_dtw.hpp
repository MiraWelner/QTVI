#pragma once

#include "common.hpp"
#include <vector>
#include <cmath>
#include <set>
#include <algorithm>

/**
 * @file draw_dtw.hpp
 * @brief Reconstruct a time-warped version of signal y2 aligned to y1 via DTW path.
 *
 * Port of draw_dtw.m.
 */

namespace ppg {

    /**
     * @brief Result of the DTW drawing/warping operation.
     */
    struct DrawDtwResult {
        std::vector<double> y_modify;   ///< Blended signal (y1 adjusted toward y2 if close).
        std::vector<double> y2_modify;  ///< Warped version of y2 resampled onto y1's time base.
        double r;                       ///< Normalised distance ratio.
    };

    /**
     * @brief Warp signal y2 onto y1's timebase using DTW alignment path.
     * @param y1    Template signal.
     * @param pla1  PLA breakpoint indices for y1 (0-based).
     * @param p     DTW path row indices (into pla1, 0-based).
     * @param y2    Beat signal.
     * @param pla2  PLA breakpoint indices for y2 (0-based).
     * @param q     DTW path column indices (into pla2, 0-based).
     * @return DrawDtwResult with warped signals and distance ratio.
     */
    inline DrawDtwResult draw_dtw(
        const std::vector<double>& y1, const std::vector<int>& pla1, const std::vector<int>& p,
        const std::vector<double>& y2, const std::vector<int>& pla2, const std::vector<int>& q)
    {
        const int l = static_cast<int>(p.size());
        if (l < 2) return { y1, y2, 1.0 };

        // Skip initial repeated entries
        int i_start = 0;
        while (i_start + 1 < l && p[i_start + 1] == p[0]) ++i_start;
        int j_start = 0;
        while (j_start + 1 < l && q[j_start + 1] == q[0]) ++j_start;
        int point = std::max(i_start, j_start) + 1;

        int la = pla1[p[0]];
        int lb = pla2[q[0]];
        int i_cur = 0, j_cur = 0;

        std::vector<double> out_x, out_y;

        while (point < l) {
            // Advance past repeated path entries
            while (point + 1 < l &&
                (p[point + 1] == p[point] || q[point + 1] == q[point]))
                ++point;

            if (point < l) {
                int la_old = la, lb_old = lb;
                la = pla1[p[point]];
                lb = pla2[q[point]];

                double intvb = (la != la_old)
                    ? static_cast<double>(lb - lb_old) / static_cast<double>(la - la_old) / 10.0
                    : 0.0;

                // Build xx: samples along y1's time axis
                std::vector<double> xx;
                if (intvb != 0.0) {
                    for (double v = pla1[p[i_cur]]; v <= pla1[p[point]] + 1e-9; v += std::abs(intvb))
                        xx.push_back(v);
                }
                if (xx.empty()) xx.push_back(static_cast<double>(pla1[p[point]]));

                // Extract yy: slice of y2
                int y2_start = pla2[q[j_cur]];
                int y2_end = pla2[q[point]];
                if (y2_start > y2_end) std::swap(y2_start, y2_end);
                std::vector<double> yy;
                for (int idx = y2_start; idx <= y2_end && idx < static_cast<int>(y2.size()); ++idx)
                    yy.push_back(y2[idx]);

                if (yy.empty()) { ++point; continue; }

                // Resample yy onto xx via spline
                std::vector<double> y_interp;
                if (xx.size() <= 1) {
                    y_interp.push_back(yy.back());
                }
                else {
                    std::vector<double> x1_coords = linspace(xx.front(), xx.back(), static_cast<int>(yy.size()));
                    y_interp = interp1_spline(x1_coords, yy, xx);
                }

                out_x.insert(out_x.end(), xx.begin(), xx.end());
                out_y.insert(out_y.end(), y_interp.begin(), y_interp.end());
                i_cur = point;
                j_cur = point;
            }
            ++point;
        }

        // Remove duplicates (keep last value at each x)
        std::vector<double> ux, uy;
        if (!out_x.empty()) {
            for (size_t k = 0; k < out_x.size(); ++k) {
                if (k + 1 < out_x.size() && out_x[k] == out_x[k + 1]) continue;
                ux.push_back(out_x[k]);
                uy.push_back(out_y[k]);
            }
        }

        // Resample to match y1 length
        std::vector<double> y2_modify;
        if (ux.size() >= 2 && y1.size() >= 2) {
            // Deduplicate ux for valid interpolation
            std::vector<double> ux_u, uy_u;
            for (size_t k = 0; k < ux.size(); ++k) {
                if (ux_u.empty() || ux[k] != ux_u.back()) {
                    ux_u.push_back(ux[k]);
                    uy_u.push_back(uy[k]);
                }
            }
            if (ux_u.size() >= 2) {
                auto xq = linspace(ux_u.front(), ux_u.back(), static_cast<int>(y1.size()));
                y2_modify = interp1_spline(ux_u, uy_u, xq);
            }
            else {
                y2_modify.assign(y1.size(), uy_u.empty() ? 0.0 : uy_u[0]);
            }
        }
        else {
            y2_modify.assign(y1.size(), 0.0);
        }

        // Replace NaN with 0
        for (auto& v : y2_modify)
            if (std::isnan(v)) v = 0.0;

        // Compute normalised distance
        double sum_sq = 0.0;
        for (size_t k = 0; k < y1.size() && k < y2_modify.size(); ++k) {
            double d = y2_modify[k] - y1[k];
            sum_sq += d * d;
        }
        double dist_val = std::sqrt(sum_sq);

        double mean_y1 = 0.0;
        for (double v : y1) mean_y1 += v * v;
        mean_y1 = std::sqrt(mean_y1);

        double r = (mean_y1 > 1e-15) ? dist_val / mean_y1 : 1.0;

        // Blend
        std::vector<double> y_modify = y1;
        if (r < 0.001) {
            for (size_t k = 0; k < y_modify.size(); ++k)
                y_modify[k] = y1[k] * 0.9 + y2_modify[k] * 0.1;
        }

        return { std::move(y_modify), std::move(y2_modify), r };
    }

} // namespace ppg