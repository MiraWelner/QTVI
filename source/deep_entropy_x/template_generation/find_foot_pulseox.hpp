/**
 * @file   find_foot_pulseox.hpp
 * @brief  Find the foot (onset) of a pulse oximetry waveform using the
 *         intersecting tangent / max derivative rotation method.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-26
 */
#pragma once

#include "template_structs.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

 // Simple local-maxima finder (no min-distance constraint).
 // Only used here for derivative peaks - not the full R-peak findpeaks.
static inline void fp_findpeaks_simple(const std::vector<double>& data,
    std::vector<double>& pks,
    std::vector<size_t>& locs) {
    pks.clear();
    locs.clear();
    if (data.size() < 3) return;
    for (size_t i = 1; i < data.size() - 1; ++i) {
        if (!std::isnan(data[i]) && data[i] > data[i - 1] && data[i] >= data[i + 1]) {
            size_t j = i;
            while (j < data.size() - 1 && data[j] == data[j + 1]) ++j;
            if (j < data.size() - 1 && data[j] > data[j + 1]) {
                pks.push_back(data[i]);
                locs.push_back(i + (j - i) / 2);
                i = j;
            }
        }
    }
}

inline FootResult find_foot_pulseox(const std::vector<std::vector<double>>& data) {
    FootResult res;
    size_t nrows = data.size();
    res.val.resize(nrows, 0.0);
    res.idx.resize(nrows, 0);

#pragma omp parallel for schedule(dynamic)
    for (int r = 0; r < static_cast<int>(nrows); ++r) {
        const auto& row = data[r];

        // Effective length (strip trailing NaNs)
        size_t len = row.size();
        while (len > 0 && std::isnan(row[len - 1])) --len;

        if (len < 4) {
            res.idx[r] = 0;
            res.val[r] = 0.0;
            continue;
        }

        // data = data - max(data)
        double row_max = -std::numeric_limits<double>::infinity();
        for (size_t k = 0; k < len; ++k) {
            if (!std::isnan(row[k]) && row[k] > row_max) row_max = row[k];
        }
        std::vector<double> d(len);
        for (size_t k = 0; k < len; ++k) {
            d[k] = std::isnan(row[k]) ? std::numeric_limits<double>::quiet_NaN() : row[k] - row_max;
        }

        // Find max position
        size_t m_pos = 0;
        double m_val = -std::numeric_limits<double>::infinity();
        for (size_t k = 0; k < len; ++k) {
            if (!std::isnan(d[k]) && d[k] > m_val) { m_val = d[k]; m_pos = k; }
        }

        // Compute diff of shifted signal
        std::vector<double> dd(len > 0 ? len - 1 : 0);
        for (size_t k = 0; k + 1 < len; ++k) {
            dd[k] = d[k + 1] - d[k];
        }

        // Find local maxima on diff, keep only those at or before the max
        std::vector<double> pks;
        std::vector<size_t> locs;
        fp_findpeaks_simple(dd, pks, locs);

        std::vector<std::pair<double, size_t>> before_peaks;
        for (size_t k = 0; k < locs.size(); ++k) {
            if (locs[k] <= m_pos) {
                before_peaks.push_back({ pks[k], locs[k] });
            }
        }

        // Greedy ascending scan for greatest diff-peak
        size_t diff_peak_loc = 0;
        double diff_peak_val = 0.0;

        if (!before_peaks.empty()) {
            std::sort(before_peaks.begin(), before_peaks.end(),
                [](const auto& a, const auto& b) { return a.first > b.first; });

            double best_val = 0.0;
            int overshoot = 0;
            for (const auto& [amp, loc] : before_peaks) {
                if (amp > best_val) {
                    best_val = amp;
                    diff_peak_val = amp;
                    diff_peak_loc = loc;
                    overshoot = 0;
                }
                else {
                    ++overshoot;
                }
                if (overshoot > 2) break;
            }
        }

        if (diff_peak_loc == 0 && !dd.empty()) {
            for (size_t k = 0; k < dd.size(); ++k) {
                if (dd[k] > diff_peak_val) { diff_peak_val = dd[k]; diff_peak_loc = k; }
            }
        }

        if (diff_peak_loc == 0) {
            res.idx[r] = 0;
            res.val[r] = d[0];
            continue;
        }

        // Translate so begin = (0, 0)
        double p1y = -d[0];
        std::vector<double> moved(diff_peak_loc + 1);
        for (size_t k = 0; k <= diff_peak_loc && k < len; ++k) {
            moved[k] = d[k] + p1y;
        }

        double ex = static_cast<double>(diff_peak_loc);
        double ey = moved[diff_peak_loc];

        // Rotation angle (match MATLAB: abs(atand(ey/ex)))
        double theta = std::atan2(ey, ex) * 180.0 / M_PI;
        if (theta < 0) theta = -theta;

        double ct = std::cos(theta * M_PI / 180.0);
        double st = std::sin(theta * M_PI / 180.0);

        // Rotate and find max of x-component (= the foot)
        double best_rx = -std::numeric_limits<double>::infinity();
        size_t best_k = 0;
        for (size_t k = 0; k <= diff_peak_loc; ++k) {
            double x = static_cast<double>(k);
            double y = moved[k];
            double rx = ct * x + st * y;
            if (rx > best_rx) { best_rx = rx; best_k = k; }
        }

        res.idx[r] = best_k;
        res.val[r] = best_rx;
    }

    return res;
}