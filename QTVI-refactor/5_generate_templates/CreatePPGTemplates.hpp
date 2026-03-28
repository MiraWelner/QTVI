/**
 * @file   CreatePPGTemplates.hpp
 * @brief  Create PPG templates for each bin using EnsembleTemplate.
 *         Port of CreatePPGTemplates.m
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-26
 */
#pragma once

#include "TemplateTypes.hpp"
#include "EnsembleTemplate.hpp"

inline vector<vector<double>> CreatePPGTemplates(
    const vector<output_binfile_data>& bins,
    double std_multiplier)
{
    size_t n = bins.size();
    vector<vector<double>> templates(n);

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < static_cast<int>(n); ++i) {
        if (bins[i].bad_segment || bins[i].ppgMinAmps.empty()) {
            templates[i] = {};
            continue;
        }
        try {
            templates[i] = EnsembleTemplate(
                bins[i].ppgSignal,
                bins[i].ppgMinAmps,
                std_multiplier,
                "ppg");
        }
        catch (...) {
            templates[i] = {};
        }
    }

    // Pad all to same length with NaN
    size_t max_len = 0;
    for (const auto& t : templates)
        if (t.size() > max_len) max_len = t.size();

    for (auto& t : templates)
        t.resize(max_len, NaN);

    return templates;
}