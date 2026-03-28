/**
 * @file   CreateEcgTemplates.hpp
 * @brief  Create ECG templates for each bin using EnsembleTemplate.
 *         Builds templates from 3 channels x 3 preprocessing methods.
 *         Port of CreateEcgTemplates.m (extended for multi-channel + 3 methods).
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-27
 */
#pragma once

#include "TemplateTypes.hpp"
#include "EnsembleTemplate.hpp"
#include "StatsUtils.h"

struct SingleMethodResult {
    vector<double> ecgTemplate;
    double ppg_alignment_point;
    double avg_r_expand;
};

static inline SingleMethodResult build_ecg_template_for_method(
    const vector<double>& ecgSignal,
    const vector<size_t>& rpeaks,
    const vector<vector<double>>& pairs,
    double std_multiplier)
{
    SingleMethodResult res;
    res.ecgTemplate = {};
    res.ppg_alignment_point = NaN;
    res.avg_r_expand = 0.0;

    if (rpeaks.empty() || ecgSignal.empty()) return res;

    // Use first ~10% of the signal for template building
    size_t reduced_size = ecgSignal.size() / 10;
    if (reduced_size < 2) reduced_size = ecgSignal.size();

    vector<double> ecg(ecgSignal.begin(),
        ecgSignal.begin() +
        std::min(reduced_size, ecgSignal.size()));

    // R peaks within the reduced region
    vector<size_t> r;
    for (auto idx : rpeaks) {
        if (idx < reduced_size - 1) r.push_back(idx);
    }
    if (r.size() < 2) return res;

    // expand = floor(diff(r) / 5)
    vector<size_t> lens;
    for (size_t j = 0; j + 1 < r.size(); ++j) {
        lens.push_back(static_cast<size_t>((r[j + 1] - r[j]) / 5));
    }

    // avg_r_expand = median(lens)
    {
        vector<double> ld(lens.begin(), lens.end());
        res.avg_r_expand = median(ld);
    }

    try {
        res.ecgTemplate = EnsembleTemplate(ecg, r, std_multiplier, "ecg", lens);
    }
    catch (...) {
        res.ecgTemplate = {};
        return res;
    }

    // ppg_alignment_point = median(pairs(:,1) - pairs(:,2))
    if (!pairs.empty()) {
        vector<double> diffs;
        for (const auto& p : pairs) {
            if (p.size() >= 2 && p[0] >= 0 && p[1] >= 0) {
                diffs.push_back(p[0] - p[1]);
            }
        }
        if (!diffs.empty()) res.ppg_alignment_point = median(diffs);
    }

    return res;
}

static inline void init_channel_result(EcgChannelResult& cr, size_t n) {
    cr.ecgTemplates_raw.resize(n);
    cr.ecgTemplates_squared.resize(n);
    cr.ecgTemplates_absval.resize(n);
    cr.ppg_alignment_point_raw.resize(n, NaN);
    cr.ppg_alignment_point_squared.resize(n, NaN);
    cr.ppg_alignment_point_absval.resize(n, NaN);
    cr.avg_r_expand_raw.resize(n, 0.0);
    cr.avg_r_expand_squared.resize(n, 0.0);
    cr.avg_r_expand_absval.resize(n, 0.0);
}

static inline void process_channel(
    EcgChannelResult& cr,
    const vector<output_binfile_data>& bins,
    size_t i,
    const vector<double>& ecgSignal,
    const ChannelRPeaks& ch,
    double std_multiplier)
{
    const auto& bin = bins[i];

    auto raw_res = build_ecg_template_for_method(ecgSignal, ch.raw, bin.pairs, std_multiplier);
    cr.ecgTemplates_raw[i] = raw_res.ecgTemplate;
    cr.ppg_alignment_point_raw[i] = raw_res.ppg_alignment_point;
    cr.avg_r_expand_raw[i] = raw_res.avg_r_expand;

    auto sq_res = build_ecg_template_for_method(ecgSignal, ch.squared, bin.pairs, std_multiplier);
    cr.ecgTemplates_squared[i] = sq_res.ecgTemplate;
    cr.ppg_alignment_point_squared[i] = sq_res.ppg_alignment_point;
    cr.avg_r_expand_squared[i] = sq_res.avg_r_expand;

    auto abs_res = build_ecg_template_for_method(ecgSignal, ch.absval, bin.pairs, std_multiplier);
    cr.ecgTemplates_absval[i] = abs_res.ecgTemplate;
    cr.ppg_alignment_point_absval[i] = abs_res.ppg_alignment_point;
    cr.avg_r_expand_absval[i] = abs_res.avg_r_expand;
}

inline EcgTemplateResult CreateEcgTemplates(
    const vector<output_binfile_data>& bins,
    double std_multiplier)
{
    size_t n = bins.size();
    EcgTemplateResult res;
    init_channel_result(res.ch1, n);
    init_channel_result(res.ch2, n);
    init_channel_result(res.ch3, n);

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < static_cast<int>(n); ++i) {
        const auto& bin = bins[i];

        // Channel 1 (always present)
        process_channel(res.ch1, bins, i, bin.ecgSignal, bin.ch1, std_multiplier);

        // Channel 2 (Bittium has this)
        if (!bin.ecgSignal2.empty()) {
            process_channel(res.ch2, bins, i, bin.ecgSignal2, bin.ch2, std_multiplier);
        }

        // Channel 3 (Bittium has this)
        if (!bin.ecgSignal3.empty()) {
            process_channel(res.ch3, bins, i, bin.ecgSignal3, bin.ch3, std_multiplier);
        }
    }

    return res;
}