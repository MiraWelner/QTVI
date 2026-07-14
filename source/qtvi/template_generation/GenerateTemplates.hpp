/**
 * @file   GenerateTemplates.hpp
 * @brief  Orchestrate the full template generation pipeline.
 *         Port of GenerateTemplates.m
 *
 *         Under Patch B, PPG (and arterial) templates are built by the
 *         same [R_i - pad, R_{i+1} + pad] slicer as the ECG templates,
 *         driven by ch1.raw R-peaks. They come out R-anchored by
 *         construction, so the old find_foot -> AlignWaves -> NaN-strip
 *         PPG alignment pipeline is gone; PPG per-sample std rides through
 *         directly.
 *
 *         Rates for every channel arrive via a SignalRates struct
 *         (defined in TemplateTypes.hpp). A rate of 0 means the channel is
 *         absent from this dataset -- the slicer silently produces empty
 *         templates for those.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-26
 */
#pragma once

#include "TemplateTypes.hpp"
#include "CreatePPGTemplates.hpp"
#include "CreateEcgTemplates.hpp"
#include <iostream>

 // FAST: PPG templates + ECG raw/unfiltered templates. The squared/absval
 // fields of each returned TemplateInfo are left empty for
 // AugmentTemplatesSlow() to fill. This is everything the viewer displays.
inline vector<TemplateInfo> GenerateTemplatesFast(const vector<output_binfile_data>& wave_data,
    const SignalRates& rates) {
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

    // PPG templates (+ per-sample std, parallel shape). Under Patch B they
    // come out already R-anchored by construction (slice = [R_i-pad, R_i+1+pad]
    // at ppgRate, R_first at column pad*ppgRate), so the old find_foot ->
    // AlignWaves -> NaN-strip pipeline is unnecessary. We just hand the
    // templates through.
    vector<vector<double>> ppg_templates;
    vector<vector<double>> ppg_template_stds;
    vector<vector<vector<double>>> ppg_kept(n);
    vector<bool> template_good(n, false);

    if (has_ppg && rates.ppg > 0.0) {
        PPGTemplatesResult ppg_res = CreatePPGTemplates(wave_data, rates.ecg, rates.ppg);
        ppg_templates = std::move(ppg_res.templates);
        ppg_template_stds = std::move(ppg_res.stds);
        ppg_kept = std::move(ppg_res.kept);

        for (size_t i = 0; i < n; ++i)
            for (double v : ppg_templates[i])
                if (!std::isnan(v)) { template_good[i] = true; break; }
    }

    // ECG templates -- FAST methods only (raw + unfiltered). The
    // squared/absval columns stay empty here; fill_channel copies those
    // empty vectors through harmlessly, and AugmentTemplatesSlow fills
    // them later.
    EcgTemplateResult ecg_res = CreateEcgTemplatesFast(wave_data, std_multiplier, rates.ecg);

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

        dst.n_beats_raw = (i < src.n_beats_raw.size()) ? src.n_beats_raw[i] : 0;
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

            if (i < ecg_res.ch1.kept_beats_raw.size())
                info.kept_beats_by_channel["CH1"] = std::move(ecg_res.ch1.kept_beats_raw[i]);
            if (i < ecg_res.ch2.kept_beats_raw.size())
                info.kept_beats_by_channel["CH2"] = std::move(ecg_res.ch2.kept_beats_raw[i]);
            if (i < ecg_res.ch3.kept_beats_raw.size())
                info.kept_beats_by_channel["CH3"] = std::move(ecg_res.ch3.kept_beats_raw[i]);
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
            if (i < ppg_kept.size()) {
                info.ppg_n_beats = ppg_kept[i].size();
                info.kept_beats_by_channel["PPG"] = std::move(ppg_kept[i]);
            }
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
    vector<TemplateInfo>& templates,
    const SignalRates& rates)
{
    constexpr double std_multiplier = 2.5;
    size_t n = wave_data.size();

    EcgTemplateResult ecg_res;
    init_channel_result(ecg_res.ch1, n);
    init_channel_result(ecg_res.ch2, n);
    init_channel_result(ecg_res.ch3, n);
    CreateEcgTemplatesSlow(wave_data, std_multiplier, rates.ecg, ecg_res);

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
inline vector<TemplateInfo> GenerateTemplates(const vector<output_binfile_data>& wave_data,
    const SignalRates& rates) {
    vector<TemplateInfo> templates = GenerateTemplatesFast(wave_data, rates);
    AugmentTemplatesSlow(wave_data, templates, rates);
    return templates;
}