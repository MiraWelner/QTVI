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
 *         std vectors: per-sample std for the ECG raw method and the
 *         PPG template are produced inside EnsembleTemplate. For PPG,
 *         the std vector has to ride through the same AlignWaves shift
 *         and NaN-strip the template gets -- otherwise the band wouldn't
 *         line up with the line.
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

 // FAST: PPG templates + ECG raw/unfiltered templates. The squared/absval
 // fields of each returned TemplateInfo are left empty for
 // AugmentTemplatesSlow() to fill. This is everything the viewer displays.
inline vector<TemplateInfo> GenerateTemplatesFast(const vector<output_binfile_data>& wave_data) {
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

    // PPG templates (+ per-sample std, parallel shape)
    vector<vector<double>> ppg_templates;
    vector<vector<double>> ppg_template_stds;
    vector<bool> template_good(n, false);

    if (has_ppg) {
        PPGTemplatesResult ppg_res = CreatePPGTemplates(wave_data, std_multiplier);
        auto& template_matrix = ppg_res.templates;
        auto& template_std_matrix = ppg_res.stds;

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

        // The std vector for each bin has to undergo the same shift as the
        // template, so std[k] stays paired with template[k] after the
        // align. We replicate AlignWaves' index math here directly --
        // there's no peak-shift (that's a vertical adjustment to the
        // template values, not the std), so this is just a horizontal
        // shift with NaN padding.
        const size_t aligned_cols = aligned.alignedWaves.empty()
            ? 0 : aligned.alignedWaves[0].size();
        vector<vector<double>> aligned_stds(n, vector<double>(aligned_cols, NaN));
        for (size_t i = 0; i < n; ++i) {
            if (!template_good[i]) continue;
            const int mv = (i < aligned.move_dist.size()) ? aligned.move_dist[i] : 0;
            const size_t dst_start = static_cast<size_t>(std::max(0, mv));
            const auto& src = template_std_matrix[i];
            for (size_t j = 0; j < src.size() && dst_start + j < aligned_cols; ++j) {
                aligned_stds[i][dst_start + j] = src[j];
            }
        }

        ppg_templates.resize(n);
        ppg_template_stds.resize(n);

#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < static_cast<int>(n); ++i) {
            if (!template_good[i]) {
                ppg_templates[i] = {};
                ppg_template_stds[i] = {};
                continue;
            }

            const auto& row = aligned.alignedWaves[i];
            const auto& sd_row = aligned_stds[i];

            // Strip NaN positions from the template -- and strip the SAME
            // positions from the std vector so they stay aligned.
            vector<double> stripped;
            vector<double> stripped_sd;
            stripped.reserve(row.size());
            stripped_sd.reserve(row.size());
            for (size_t k = 0; k < row.size(); ++k) {
                if (!std::isnan(row[k])) {
                    stripped.push_back(row[k]);
                    // sd_row may be NaN in padding columns; 0 is a safe
                    // replacement (band collapses to the line there).
                    stripped_sd.push_back(std::isnan(sd_row[k]) ? 0.0 : sd_row[k]);
                }
            }
            ppg_templates[i] = std::move(stripped);
            ppg_template_stds[i] = std::move(stripped_sd);
        }
    }

    // ECG templates -- FAST methods only (raw + unfiltered). The
    // squared/absval columns stay empty here; fill_channel copies those
    // empty vectors through harmlessly, and AugmentTemplatesSlow fills
    // them later.
    EcgTemplateResult ecg_res = CreateEcgTemplatesFast(wave_data, std_multiplier);

    // Assemble TemplateInfo
    auto fill_channel = [](ChannelTemplates& dst, const EcgChannelResult& src, size_t i) {
        dst.ecgTemplate_raw = src.ecgTemplates_raw[i];
        // Per-sample std for the raw method only -- the other three
        // methods are never displayed, so they don't have std computed.
        if (i < src.ecgTemplates_raw_std.size())
            dst.ecgTemplate_raw_std = src.ecgTemplates_raw_std[i];

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
        dst.ecgTemplate_raw_std = {};
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
        // ECG quality is independent of PPG availability. Previously this
        // function gated ECG fill on ppg_template_good, which cleared the
        // ECG templates for every bin of datasets without PPG (Bittium).
        // The top-of-file comment block documents the correct gate as
        // bad_segment, so we use that here.
        const bool ecg_good = (i < wave_data.size()) && !wave_data[i].bad_segment;

        if (ecg_good) {
            fill_channel(info.ch1, ecg_res.ch1, i);
            fill_channel(info.ch2, ecg_res.ch2, i);
            fill_channel(info.ch3, ecg_res.ch3, i);

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

        if (ppg_template_good) {
            info.ppgTemplate = ppg_templates[i];
            info.ppgTemplate_std = ppg_template_stds[i];
        }
        // else: info.ppgTemplate / ppgTemplate_std stay default-empty,
        // which the viewer already interprets as "no PPG for this bin".
    }
    return result;
}

// SLOW: fill the squared/absval ECG templates onto an existing
// vector<TemplateInfo> produced by GenerateTemplatesFast. Applies the same
// bad_segment gate as the fast pass, so squared/absval stay empty on bins
// the fast pass cleared.
inline void AugmentTemplatesSlow(const vector<output_binfile_data>& wave_data,
    vector<TemplateInfo>& templates)
{
    constexpr double std_multiplier = 2.5;
    size_t n = wave_data.size();

    EcgTemplateResult ecg_res;
    init_channel_result(ecg_res.ch1, n);
    init_channel_result(ecg_res.ch2, n);
    init_channel_result(ecg_res.ch3, n);
    CreateEcgTemplatesSlow(wave_data, std_multiplier, ecg_res);

    auto fill_slow = [](ChannelTemplates& dst, const EcgChannelResult& src, size_t i) {
        dst.ecgTemplate_squared = src.ecgTemplates_squared[i];
        dst.ecgTemplate_absval = src.ecgTemplates_absval[i];
        dst.alignment_point_squared = std::isnan(src.ppg_alignment_point_squared[i]) ? 0.0 : src.ppg_alignment_point_squared[i];
        dst.alignment_point_absval = std::isnan(src.ppg_alignment_point_absval[i]) ? 0.0 : src.ppg_alignment_point_absval[i];
        dst.avg_r_expand_squared = src.avg_r_expand_squared[i];
        dst.avg_r_expand_absval = src.avg_r_expand_absval[i];
        };

    for (size_t i = 0; i < n && i < templates.size(); ++i) {
        const bool ecg_good = (i < wave_data.size()) && !wave_data[i].bad_segment;
        if (!ecg_good) continue;
        fill_slow(templates[i].ch1, ecg_res.ch1, i);
        fill_slow(templates[i].ch2, ecg_res.ch2, i);
        fill_slow(templates[i].ch3, ecg_res.ch3, i);
    }
}

// Original all-methods entry point, preserved by composition.
inline vector<TemplateInfo> GenerateTemplates(const vector<output_binfile_data>& wave_data) {
    vector<TemplateInfo> templates = GenerateTemplatesFast(wave_data);
    AugmentTemplatesSlow(wave_data, templates);
    return templates;
}