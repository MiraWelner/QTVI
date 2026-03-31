/**
 * @file   CreateEcgTemplates.hpp
 * @brief  Create ECG templates for each bin using EnsembleTemplate.
 *         Builds templates from 3 channels x 3 preprocessing methods.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-27
 */
#pragma once

#include "TemplateTypes.hpp"
#include "EnsembleTemplate.hpp"
#include "StatsUtils.h"
#include <atomic>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

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

    if (rpeaks.size() < 2 || ecgSignal.empty()) return res;

    // Peaks are already indices into the bin's signal — no extraction needed.
    // Just filter out any that are out of bounds.
    vector<size_t> r;
    r.reserve(rpeaks.size());
    for (auto idx : rpeaks)
        if (idx < ecgSignal.size()) r.push_back(idx);
    if (r.size() < 2) return res;

    vector<size_t> lens;
    lens.reserve(r.size() - 1);
    for (size_t j = 0; j + 1 < r.size(); ++j)
        lens.push_back(static_cast<size_t>((r[j + 1] - r[j]) / 5));

    {
        vector<double> ld(lens.begin(), lens.end());
        res.avg_r_expand = median(ld);
    }


    try {
        res.ecgTemplate = EnsembleTemplate(ecgSignal, r, std_multiplier, "ecg", lens);
    }
    catch (...) {
        res.ecgTemplate = {};
        return res;
    }

    if (!res.ecgTemplate.empty()) {
        double minVal = *std::min_element(res.ecgTemplate.begin(), res.ecgTemplate.end());
        for (auto& v : res.ecgTemplate) v -= minVal;
    }

    if (!pairs.empty()) {
        vector<double> diffs;
        for (const auto& p : pairs) {
            if (p.size() >= 2 && p[0] >= 0 && p[1] >= 0)
                diffs.push_back(p[0] - p[1]);
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

    const auto& sq_sig = ch.squared_signal.empty() ? ecgSignal : ch.squared_signal;
    auto sq_res = build_ecg_template_for_method(sq_sig, ch.squared, bin.pairs, std_multiplier);
    cr.ecgTemplates_squared[i] = sq_res.ecgTemplate;
    cr.ppg_alignment_point_squared[i] = sq_res.ppg_alignment_point;
    cr.avg_r_expand_squared[i] = sq_res.avg_r_expand;

    const auto& abs_sig = ch.absval_signal.empty() ? ecgSignal : ch.absval_signal;
    auto abs_res = build_ecg_template_for_method(abs_sig, ch.absval, bin.pairs, std_multiplier);
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

    std::atomic<int> done{ 0 };

    int max_threads = std::min(8, (int)n);
    #pragma omp parallel for schedule(dynamic) num_threads(max_threads)
    for (int i = 0; i < static_cast<int>(n); ++i) {

        auto bin_start = std::chrono::steady_clock::now();

        const auto& bin = bins[i];

        process_channel(res.ch1, bins, i, bin.ecgSignal, bin.ch1, std_multiplier);

        if (!bin.ecgSignal2.empty())
            process_channel(res.ch2, bins, i, bin.ecgSignal2, bin.ch2, std_multiplier);

        if (!bin.ecgSignal3.empty())
            process_channel(res.ch3, bins, i, bin.ecgSignal3, bin.ch3, std_multiplier);

        auto bin_end = std::chrono::steady_clock::now();
        double bin_sec = std::chrono::duration<double>(bin_end - bin_start).count();

        int d = ++done;
    }
    return res;
}