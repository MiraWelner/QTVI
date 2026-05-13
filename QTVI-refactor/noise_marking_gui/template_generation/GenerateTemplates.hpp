/**
 * @file   GenerateTemplates.hpp
 * @brief  Orchestrate the full template generation pipeline.
 *         Port of GenerateTemplates.m
 *
 *         Optimizations vs original:
 *           - Skip WaveDiff (O(n^2)) and CombineTemplatesGraph since
 *             combine assigns each template its own bin anyway.
 *           - Skip full PPG align/diff/combine pipeline; just use
 *             per-bin PPG templates directly after foot-stripping.
 *           - ECG gate on bad_segment, not PPG availability.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-26
 */
#pragma once

#include "TemplateTypes.hpp"
#include "CreatePPGTemplates.hpp"
#include "CreateEcgTemplates.hpp"
#include "find_foot_pulseox.hpp"
#include "AlignWaves.hpp"
#include <iostream>

inline vector<TemplateInfo> GenerateTemplates(const vector<output_binfile_data>& wave_data) {
    constexpr double std_multiplier = 2.5;

    size_t n = wave_data.size();

    // Check if any bin has PPG data
    bool has_ppg = false;
    for (size_t i = 0; i < n; ++i) {
        if (!wave_data[i].ppgSignal.empty() && !wave_data[i].ppgMinAmps.empty()) {
            has_ppg = true;
            break;
        }
    }

    // PPG templates
    vector<vector<double>> ppg_templates;
    vector<bool> template_good(n, false);

    if (has_ppg) {
        vector<vector<double>> template_matrix = CreatePPGTemplates(wave_data, std_multiplier);

        for (size_t i = 0; i < n; ++i) {
            for (double v : template_matrix[i]) {
                if (!std::isnan(v)) { template_good[i] = true; break; }
            }
        }

        // MATLAB pipeline: find_foot → AlignWaves → CombineTemplatesGraph.
        // The graph-combine step does the NaN-stripping that produces the
        // final per-bin PPG templates.
        FootResult feet = find_foot_pulseox(template_matrix);

        // AlignWaves: shifts each template so its foot lands at the global
        // max foot column, padding with NaN.
        AlignWavesResult aligned = AlignWaves(template_matrix, feet.idx);
        ppg_templates.resize(n);

#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < static_cast<int>(n); ++i) {
            if (!template_good[i]) {
                ppg_templates[i] = {};
                continue;
            }

            const auto& row = aligned.alignedWaves[i];
            vector<double> stripped;
            stripped.reserve(row.size());
            for (double v : row) {
                if (!std::isnan(v)) stripped.push_back(v);
            }
            ppg_templates[i] = std::move(stripped);
        }
    }

    // ECG templates (4 methods x 3 channels per bin)
    EcgTemplateResult ecg_res = CreateEcgTemplates(wave_data, std_multiplier);

    // Assemble TemplateInfo
    auto fill_channel = [](ChannelTemplates& dst, const EcgChannelResult& src, size_t i) {
        dst.ecgTemplate_raw = src.ecgTemplates_raw[i];
        dst.ecgTemplate_squared = src.ecgTemplates_squared[i];
        dst.ecgTemplate_absval = src.ecgTemplates_absval[i];
        dst.ecgTemplate_unfiltered = src.ecgTemplates_unfiltered[i];

        dst.alignment_point_raw = std::isnan(src.ppg_alignment_point_raw[i]) ? 0.0 : src.ppg_alignment_point_raw[i];
        dst.alignment_point_squared = std::isnan(src.ppg_alignment_point_squared[i]) ? 0.0 : src.ppg_alignment_point_squared[i];
        dst.alignment_point_absval = std::isnan(src.ppg_alignment_point_absval[i]) ? 0.0 : src.ppg_alignment_point_absval[i];
        dst.alignment_point_unfiltered = std::isnan(src.ppg_alignment_point_unfiltered[i]) ? 0.0 : src.ppg_alignment_point_unfiltered[i];

        dst.avg_r_expand_raw = src.avg_r_expand_raw[i];
        dst.avg_r_expand_squared = src.avg_r_expand_squared[i];
        dst.avg_r_expand_absval = src.avg_r_expand_absval[i];
        dst.avg_r_expand_unfiltered = src.avg_r_expand_unfiltered[i];
        };

    auto clear_channel = [](ChannelTemplates& dst) {
        dst.ecgTemplate_raw = {};
        dst.ecgTemplate_squared = {};
        dst.ecgTemplate_absval = {};
        dst.ecgTemplate_unfiltered = {};

        dst.alignment_point_raw = NaN;
        dst.alignment_point_squared = NaN;
        dst.alignment_point_absval = NaN;
        dst.alignment_point_unfiltered = NaN;

        dst.avg_r_expand_raw = 0.0;
        dst.avg_r_expand_squared = 0.0;
        dst.avg_r_expand_absval = 0.0;
        dst.avg_r_expand_unfiltered = 0.0;
        };

    vector<TemplateInfo> result(n);
    for (size_t i = 0; i < n; ++i) {
        auto& info = result[i];

        const bool ppg_template_good = has_ppg && template_good[i];

        if (ppg_template_good) {
            fill_channel(info.ch1, ecg_res.ch1, i);
            fill_channel(info.ch2, ecg_res.ch2, i);
            fill_channel(info.ch3, ecg_res.ch3, i);
            info.ppgTemplate = ppg_templates[i];

            if (i < ecg_res.ch1.kept_beats_raw.size()) {
                info.kept_beats_ch1_raw = std::move(ecg_res.ch1.kept_beats_raw[i]);
            }
        }
        else {
            clear_channel(info.ch1);
            clear_channel(info.ch2);
            clear_channel(info.ch3);
            info.kept_beats_ch1_raw.clear();
        }
    }
    return result;
}