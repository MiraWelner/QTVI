#pragma once

#include "common.hpp"
#include <vector>
#include <cmath>
#include <limits>

/**
 * @file dp_dtw.hpp
 * @brief Dynamic programming to find min-cost DTW path through a cost matrix.
 *
 * Port of dp_dtw2.m (dpwe@ee.columbia.edu, 2003).
 */

namespace ppg {

    /**
     * @brief Result of the DP-based DTW path search.
     */
    struct DtwPathResult {
        std::vector<int> p;  ///< Row indices of the optimal path (0-based).
        std::vector<int> q;  ///< Column indices of the optimal path (0-based).
        Matrix D;            ///< Accumulated cost matrix.
    };

    /**
     * @brief Find minimum-cost path through cost matrix M using dynamic programming.
     * @param M    Cost matrix (rows × cols).
     * @param t_a  Time-axis weights for rows.
     * @param t_b  Time-axis weights for columns.
     * @return DtwPathResult with path indices (p, q) and accumulated cost matrix D.
     */
    inline DtwPathResult dp_dtw(
        const Matrix& M,
        const std::vector<double>& t_a,
        const std::vector<double>& t_b)
    {
        const int r = M.rows;
        const int c = M.cols;

        // Accumulated cost matrix with NaN border (use kInf as sentinel)
        Matrix D(r + 1, c + 1, kNaN);
        D(0, 0) = 0.0;
        for (int i = 0; i < r; ++i)
            for (int j = 0; j < c; ++j)
                D(i + 1, j + 1) = M(i, j);

        // Traceback matrix: 1 = diagonal, 2 = up, 3 = left
        Matrix phi(r, c, 0.0);

        for (int i = 0; i < r; ++i) {
            for (int j = 0; j < c; ++j) {
                double mij = M(i, j);
                double cost_diag = D(i, j) + mij * (t_a[i] + t_b[j]);
                double cost_up = D(i, j + 1) + mij * t_a[i];
                double cost_left = D(i + 1, j) + mij * t_b[j];

                // Handle NaN borders: treat NaN as infinity
                auto safe = [](double v) { return std::isnan(v) ? kInf : v; };
                double sd = safe(cost_diag), su = safe(cost_up), sl = safe(cost_left);

                double dmin;
                int tb;
                if (sd <= su && sd <= sl) { dmin = cost_diag; tb = 1; }
                else if (su <= sd && su <= sl) { dmin = cost_up;   tb = 2; }
                else { dmin = cost_left;  tb = 3; }

                D(i + 1, j + 1) = D(i + 1, j + 1) + dmin;
                phi(i, j) = static_cast<double>(tb);
            }
        }

        // Traceback from bottom-right
        int i = r - 1, j = c - 1;
        std::vector<int> p_path = { i };
        std::vector<int> q_path = { j };

        while (i > 0 && j > 0) {
            int tb = static_cast<int>(phi(i, j));
            if (tb == 1) { --i; --j; }
            else if (tb == 2) { --i; }
            else { --j; }
            p_path.push_back(i);
            q_path.push_back(j);
        }

        // Prepend start
        if (p_path.back() != 0 || q_path.back() != 0) {
            p_path.push_back(0);
            q_path.push_back(0);
        }

        std::reverse(p_path.begin(), p_path.end());
        std::reverse(q_path.begin(), q_path.end());

        // Strip border from D
        Matrix D_out(r, c);
        for (int ii = 0; ii < r; ++ii)
            for (int jj = 0; jj < c; ++jj)
                D_out(ii, jj) = D(ii + 1, jj + 1);

        return { std::move(p_path), std::move(q_path), std::move(D_out) };
    }

} // namespace ppg