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
        // Find feet and strip leading NaN before foot
        FootResult feet = find_foot_pulseox(template_matrix);

        ppg_templates.resize(n);

#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < static_cast<int>(n); ++i) {
            if (!template_good[i]) {
                ppg_templates[i] = {};
                continue;
            }

            const auto& tmpl = template_matrix[i];
            size_t foot = feet.idx[i];

            // Strip leading NaN before foot, trailing NaN after end
            size_t first_valid = tmpl.size();
            size_t last_valid = 0;
            for (size_t c = 0; c < tmpl.size(); ++c) {
                if (!std::isnan(tmpl[c])) {
                    if (c < first_valid) first_valid = c;
                    last_valid = c;
                }
            }

            if (first_valid > last_valid) {
                ppg_templates[i] = {};
                continue;
            }

            // Use foot as start if it's within valid range
            size_t start = (foot >= first_valid && foot <= last_valid) ? foot : first_valid;
            ppg_templates[i].reserve(last_valid - start + 1);
            for (size_t c = start; c <= last_valid; ++c) {
                ppg_templates[i].push_back(std::isnan(tmpl[c]) ? 0.0 : tmpl[c]);
            }
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
        info.index = i;
        info.ppg_bin_indexs = wave_data[i].ppg_bin_indexs;
        info.ecg_bin_indexs = wave_data[i].ecg_bin_indexs;
        info.bad_segment = wave_data[i].bad_segment;

        if (has_ppg && template_good[i]) {
            info.ppgTemplate = ppg_templates[i];
        }
        else {
            info.ppgTemplate = {};
        }

        if (!wave_data[i].bad_segment) {
            fill_channel(info.ch1, ecg_res.ch1, i);
            fill_channel(info.ch2, ecg_res.ch2, i);
            fill_channel(info.ch3, ecg_res.ch3, i);
        }
        else {
            clear_channel(info.ch1);
            clear_channel(info.ch2);
            clear_channel(info.ch3);
        }
    }
    return result;
}