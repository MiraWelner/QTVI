/**
 * @file   AlignWaves.hpp
 * @brief  Align matrix of waves to a specified alignment point.
 *         Port of AlignWaves.m (Author: Daniel Wendelken)
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-26
 */
#pragma once
#include <vector>
#include "TemplateTypes.hpp"

inline AlignWavesResult AlignWaves(const vector<vector<double>>& waves,
    const vector<size_t>& wave_alignment_points)
{
    AlignWavesResult res;
    size_t wave_count = waves.size();
    if (wave_count == 0 || wave_alignment_points.size() != wave_count) return res;

    // Length of each wave (ignoring trailing NaNs)
    vector<size_t> lengths(wave_count);
    for (size_t i = 0; i < wave_count; ++i) {
        size_t len = waves[i].size();
        while (len > 0 && std::isnan(waves[i][len - 1])) --len;
        lengths[i] = len;
    }

    size_t max_left = *std::max_element(wave_alignment_points.begin(),
        wave_alignment_points.end());

    size_t max_right = 0;
    for (size_t i = 0; i < wave_count; ++i) {
        size_t right = (lengths[i] > wave_alignment_points[i])
            ? lengths[i] - wave_alignment_points[i]
            : 0;
        if (right > max_right) max_right = right;
    }

    res.move_dist.resize(wave_count);
    for (size_t i = 0; i < wave_count; ++i) {
        res.move_dist[i] = static_cast<int>(max_left) -
            static_cast<int>(wave_alignment_points[i]);
    }

    size_t total_cols = max_left + max_right;
    res.alignedWaves.assign(wave_count, vector<double>(total_cols, NaN));

    for (size_t i = 0; i < wave_count; ++i) {
        size_t dst_start = static_cast<size_t>(std::max(0, res.move_dist[i]));
        for (size_t j = 0; j < lengths[i] && dst_start + j < total_cols; ++j) {
            res.alignedWaves[i][dst_start + j] = waves[i][j];
        }
    }
    return res;
}