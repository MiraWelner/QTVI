/**
 * @file   GenerateTemplates.hpp
 * @brief  Orchestrate the full template generation pipeline.
 *         Port of GenerateTemplates.m
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
#include "WaveDiff.hpp"
#include "CombineTemplatesGraph.hpp"
#include <iostream>

inline vector<TemplateInfo> GenerateTemplates(const vector<output_binfile_data>& wave_data) {
    constexpr double std_multiplier = 2.5;
    constexpr double threshold_percent = 15.0;

    size_t n = wave_data.size();

    // Check if any bin has PPG data
    bool has_ppg = false;
    for (size_t i = 0; i < n; ++i) {
        if (!wave_data[i].ppgSignal.empty() && !wave_data[i].ppgMinAmps.empty()) {
            has_ppg = true;
            break;
        }
    }

    vector<bool> template_good(n, false);
    CombineResult combined;

    if (has_ppg) {
        // 1. PPG templates (n x max_cols, NaN-padded)
        std::cerr << "    PPG templates (" << n << " bins)..." << std::endl;
        vector<vector<double>> template_matrix = CreatePPGTemplates(wave_data, std_multiplier);

        // Which templates have actual data
        for (size_t i = 0; i < n; ++i) {
            for (double v : template_matrix[i]) {
                if (!std::isnan(v)) { template_good[i] = true; break; }
            }
        }

        // 2. Find template feet
        std::cerr << "    Finding template feet..." << std::endl;
        FootResult feet = find_foot_pulseox(template_matrix);

        // 3. Align by foot
        std::cerr << "    Aligning templates..." << std::endl;
        vector<size_t> foot_idx(n);
        for (size_t i = 0; i < n; ++i) foot_idx[i] = feet.idx[i];
        AlignWavesResult aligned = AlignWaves(template_matrix, foot_idx);

        // 4. Diff matrix
        std::cerr << "    Computing diff matrix..." << std::endl;
        auto diff_matrix = WaveDiff(aligned.alignedWaves);

        // 5. Combine
        std::cerr << "    Combining templates..." << std::endl;
        combined = CombineTemplatesGraph(
            aligned.alignedWaves, diff_matrix, threshold_percent);
    }
    else {
        std::cerr << "    No PPG data found, skipping PPG pipeline." << std::endl;
    }

    // 6. ECG templates (3 methods per bin)
    std::cerr << "    ECG templates (3 methods x " << n << " bins)..." << std::endl;
    EcgTemplateResult ecg_res = CreateEcgTemplates(wave_data, std_multiplier);

    // 7. Assemble TemplateInfo
    std::cerr << "    Assembling output..." << std::endl;

    auto fill_channel = [](ChannelTemplates& dst, const EcgChannelResult& src, size_t i) {
        dst.ecgTemplate_raw = src.ecgTemplates_raw[i];
        dst.ecgTemplate_squared = src.ecgTemplates_squared[i];
        dst.ecgTemplate_absval = src.ecgTemplates_absval[i];
        dst.alignment_point_raw = std::isnan(src.ppg_alignment_point_raw[i]) ? 0.0 : src.ppg_alignment_point_raw[i];
        dst.alignment_point_squared = std::isnan(src.ppg_alignment_point_squared[i]) ? 0.0 : src.ppg_alignment_point_squared[i];
        dst.alignment_point_absval = std::isnan(src.ppg_alignment_point_absval[i]) ? 0.0 : src.ppg_alignment_point_absval[i];
        dst.avg_r_expand_raw = src.avg_r_expand_raw[i];
        dst.avg_r_expand_squared = src.avg_r_expand_squared[i];
        dst.avg_r_expand_absval = src.avg_r_expand_absval[i];
        };

    auto clear_channel = [](ChannelTemplates& dst) {
        dst.ecgTemplate_raw = {};
        dst.ecgTemplate_squared = {};
        dst.ecgTemplate_absval = {};
        dst.alignment_point_raw = NaN;
        dst.alignment_point_squared = NaN;
        dst.alignment_point_absval = NaN;
        dst.avg_r_expand_raw = 0.0;
        dst.avg_r_expand_squared = 0.0;
        dst.avg_r_expand_absval = 0.0;
        };

    vector<TemplateInfo> result(n);
    for (size_t i = 0; i < n; ++i) {
        auto& info = result[i];
        info.index = i;
        info.ppg_bin_indexs = wave_data[i].ppg_bin_indexs;
        info.ecg_bin_indexs = wave_data[i].ecg_bin_indexs;
        info.bad_segment = wave_data[i].bad_segment;

        if (has_ppg && template_good[i]) {
            size_t bin_id = combined.bin_numbers[i];
            info.ppgTemplate = (bin_id < combined.bin_templates.size())
                ? combined.bin_templates[bin_id]
                : vector<double>{};
        }
        else {
            info.ppgTemplate = {};
        }

        // ECG templates are always filled regardless of PPG
        bool has_ecg = !wave_data[i].bad_segment;
        if (has_ecg) {
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