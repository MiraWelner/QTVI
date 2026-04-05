#pragma once

#include "common.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>

/**
 * @file find_foot_pulseox.hpp
 * @brief Detect the "foot" (onset) of a PPG pulse beat using
 *        tangent-intersection via rotation.
 *
 * Port of find_foot_pulseox.m.
*/

const double M_PI = 3.14159265358979;

namespace ppg {

    /**
     * @brief Result of foot detection.
     */
    struct FootResult {
        double value;  ///< Value of the rotated-max metric.
        int    index;  ///< 0-based index of the foot within the beat.
    };

    /**
     * @brief Find peaks in the first-difference of data that precede the max.
     * @param data  A single beat waveform (1-D vector).
     * @return Pair of vectors: (peak_amplitudes, peak_indices) — 0-based indices.
     */
    inline std::pair<std::vector<double>, std::vector<int>> find_diff_peaks(
        const std::vector<double>& data)
    {
        auto d = diff(data);
        std::vector<double> amps;
        std::vector<int>    locs;

        for (int i = 1; i + 1 < static_cast<int>(d.size()); ++i) {
            if (d[i] > d[i - 1] && d[i] > d[i + 1]) {
                amps.push_back(d[i]);
                locs.push_back(i);
            }
        }
        return { amps, locs };
    }

    /**
     * @brief Detect the foot (onset) index of a single PPG beat.
     * @param beat  PPG beat waveform (single pulse, valley to valley).
     * @return FootResult with the foot index and associated metric value.
     */
    inline FootResult find_foot_pulseox(const std::vector<double>& beat) {
        if (beat.empty() || beat.size() < 4)
            return { 0.0, 0 };

        const int n = static_cast<int>(beat.size());

        // Normalise: subtract max so waveform is ≤ 0
        double mx_val = *std::max_element(beat.begin(), beat.end());
        std::vector<double> data(n);
        for (int i = 0; i < n; ++i) data[i] = beat[i] - mx_val;

        // Find location of max in original
        int m_idx = static_cast<int>(std::max_element(beat.begin(), beat.end()) - beat.begin());

        // Find diff-peaks before the max
        auto [peak_amps, peak_locs] = find_diff_peaks(data);

        std::vector<double> pre_amps;
        std::vector<int>    pre_locs;
        for (size_t k = 0; k < peak_locs.size(); ++k) {
            if (peak_locs[k] <= m_idx) {
                pre_amps.push_back(peak_amps[k]);
                pre_locs.push_back(peak_locs[k]);
            }
        }

        // Find the dominant diff-peak (largest amplitude, allowing small overshoot)
        int best_loc = 0;
        double best_amp = 0.0;
        if (!pre_amps.empty()) {
            // Sort descending by amplitude
            std::vector<size_t> order(pre_amps.size());
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(),
                [&](size_t a, size_t b) { return pre_amps[a] > pre_amps[b]; });

            int overshoot = 0;
            for (size_t idx : order) {
                if (pre_amps[idx] > best_amp) {
                    best_amp = pre_amps[idx];
                    best_loc = pre_locs[idx];
                    overshoot = 0;
                }
                else {
                    ++overshoot;
                }
                if (overshoot > 2) break;
            }
        }

        if (best_loc == 0) {
            // Fallback: max of diff
            auto d = diff(data);
            if (!d.empty()) {
                best_loc = static_cast<int>(std::max_element(d.begin(), d.end()) - d.begin());
            }
        }

        int end_pt = best_loc;
        if (end_pt < 1) end_pt = 1;

        // Tangent-intersection via rotation
        // Move begin point to origin
        double begin_y = data[0];
        double end_x = static_cast<double>(end_pt);
        double end_y = data[end_pt] - begin_y;

        std::vector<double> moved(end_pt + 1);
        for (int i = 0; i <= end_pt; ++i) moved[i] = data[i] - begin_y;

        // Rotation angle
        double theta = std::atan2(end_x, end_y) * 180.0 / M_PI;
        if (theta < 0.0) theta = -theta;
        double rad = theta * M_PI / 180.0;
        double cos_t = std::cos(rad), sin_t = std::sin(rad);

        // Rotate the curve and find the max of the x-component
        double max_rot = -kInf;
        int foot_idx = 0;
        for (int i = 0; i <= end_pt; ++i) {
            double rx = cos_t * (i + 1) + (-sin_t) * moved[i];
            if (rx > max_rot) {
                max_rot = rx;
                foot_idx = i;
            }
        }

        return { max_rot, foot_idx };
    }

} // namespace ppg