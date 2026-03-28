/**
 * @file   WaveDiff.hpp
 * @brief  Compute squared-sum difference matrix between aligned wave templates.
 *         Port of WaveDiff.m
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-26
 */
#pragma once

#include "TemplateTypes.hpp"

inline vector<vector<double>> WaveDiff(const vector<vector<double>>& waves) {
    size_t n = waves.size();
    vector<vector<double>> diff_matrix(n, vector<double>(n, 0.0));

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < static_cast<int>(n); ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            double s = 0.0;
            size_t cols = std::max(waves[i].size(), waves[j].size());
            for (size_t k = 0; k < cols; ++k) {
                double vi = (k < waves[i].size()) ? waves[i][k] : NaN;
                double vj = (k < waves[j].size()) ? waves[j][k] : NaN;
                if (!std::isnan(vi) && !std::isnan(vj)) {
                    s += (vi - vj);
                }
            }
            diff_matrix[i][j] = s * s;
        }
    }
    return diff_matrix;
}