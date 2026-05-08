/**
* @file   CombineTemplatesGraph.hpp
* @brief  Group similar templates using connected-component analysis on a
*         difference matrix.  Port of CombineTemplatesGraph.m
*
*         NOTE: the original MATLAB code computed connected components via
*         graph/conncomp, but then immediately overrode that with
*         bin_numbers = 1:N (each template gets its own bin).  We replicate
*         that final behaviour here.
*
* @author Mira Welner
* @email  MEW386@pitt.edu
* @date   2026-03-26
*/
#pragma once

#include "TemplateTypes.hpp"

inline CombineResult CombineTemplatesGraph(
    const vector<vector<double>>& aligned_templates,
    const vector<vector<double>>& /*diff_matrix*/,
    double /*threshold_percent*/)
{
    CombineResult res;
    size_t n = aligned_templates.size();
    if (n == 0) return res;

    // Each template is its own bin (matches final MATLAB behaviour)
    res.bin_numbers.resize(n);
    for (size_t i = 0; i < n; ++i) res.bin_numbers[i] = i;

    res.bin_templates.resize(n);
    res.foot_locations.resize(n, 0);

    for (size_t i = 0; i < n; ++i) {
        // Strip leading/trailing NaNs
        size_t first_valid = aligned_templates[i].size();
        size_t last_valid = 0;
        for (size_t c = 0; c < aligned_templates[i].size(); ++c) {
            if (!std::isnan(aligned_templates[i][c])) {
                if (c < first_valid) first_valid = c;
                last_valid = c;
            }
        }

        vector<double> tmpl;
        if (first_valid <= last_valid) {
            for (size_t c = first_valid; c <= last_valid; ++c) {
                tmpl.push_back(std::isnan(aligned_templates[i][c])
                    ? 0.0
                    : aligned_templates[i][c]);
            }
        }
        res.bin_templates[i] = tmpl;
        res.foot_locations[i] = 0;
    }

    return res;
}