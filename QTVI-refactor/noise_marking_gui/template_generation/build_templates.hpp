/**
* @file   build_templates.hpp
* @brief  Run the template-generation pipeline on peak-finding output
*         and produce a TemplateFile ready to write to disk.
*/

#pragma once

#include <cmath>
#include <vector>

#include "template_io.hpp"
#include "GenerateTemplates.hpp"
#include "TemplateTypes.hpp"
#include "peak_finding/peakfinding_io.hpp"

namespace template_generation_detail {

    inline void copyMethod(template_io::ChannelMethodTemplate& dst,
        const std::vector<double>& tmpl,
        double alignment, double rExpand)
    {
        dst.ecgTemplate = tmpl;
        dst.alignment_point = std::isnan(alignment) ? 0.0 : alignment;
        dst.avg_r_expand = rExpand;
    }

    inline void packBin(template_io::BinTemplates& bt,
        const TemplateInfo& info, bool bad_segment)
    {
        bt.bad_segment = bad_segment;
        if (bad_segment) return;

        copyMethod(bt.ch1_raw, info.ch1.ecgTemplate_raw,
            info.ch1.alignment_point_raw, info.ch1.avg_r_expand_raw);
        copyMethod(bt.ch1_squared, info.ch1.ecgTemplate_squared,
            info.ch1.alignment_point_squared, info.ch1.avg_r_expand_squared);
        copyMethod(bt.ch1_absval, info.ch1.ecgTemplate_absval,
            info.ch1.alignment_point_absval, info.ch1.avg_r_expand_absval);
        copyMethod(bt.ch1_unfiltered, info.ch1.ecgTemplate_unfiltered,
            info.ch1.alignment_point_unfiltered, info.ch1.avg_r_expand_unfiltered);

        copyMethod(bt.ch2_raw, info.ch2.ecgTemplate_raw,
            info.ch2.alignment_point_raw, info.ch2.avg_r_expand_raw);
        copyMethod(bt.ch2_squared, info.ch2.ecgTemplate_squared,
            info.ch2.alignment_point_squared, info.ch2.avg_r_expand_squared);
        copyMethod(bt.ch2_absval, info.ch2.ecgTemplate_absval,
            info.ch2.alignment_point_absval, info.ch2.avg_r_expand_absval);
        copyMethod(bt.ch2_unfiltered, info.ch2.ecgTemplate_unfiltered,
            info.ch2.alignment_point_unfiltered, info.ch2.avg_r_expand_unfiltered);

        copyMethod(bt.ch3_raw, info.ch3.ecgTemplate_raw,
            info.ch3.alignment_point_raw, info.ch3.avg_r_expand_raw);
        copyMethod(bt.ch3_squared, info.ch3.ecgTemplate_squared,
            info.ch3.alignment_point_squared, info.ch3.avg_r_expand_squared);
        copyMethod(bt.ch3_absval, info.ch3.ecgTemplate_absval,
            info.ch3.alignment_point_absval, info.ch3.avg_r_expand_absval);
        copyMethod(bt.ch3_unfiltered, info.ch3.ecgTemplate_unfiltered,
            info.ch3.alignment_point_unfiltered, info.ch3.avg_r_expand_unfiltered);

        bt.ppgTemplate = info.ppgTemplate;
    }

    inline template_io::AveragedTemplate
        averageOne(const std::vector<template_io::BinTemplates>& bins,
            std::vector<double>(*pick)(const template_io::BinTemplates&))
    {
        std::vector<std::vector<double>> bucket;
        for (const auto& b : bins) {
            if (b.bad_segment) continue;
            std::vector<double> t = pick(b);
            if (!t.empty()) bucket.push_back(std::move(t));
        }

        template_io::AveragedTemplate avg;
        if (bucket.empty()) return avg;

        size_t maxLen = 0;
        for (const auto& t : bucket) maxLen = std::max(maxLen, t.size());
        avg.waveform.assign(maxLen, 0.0);
        std::vector<uint64_t> counts(maxLen, 0);
        for (const auto& t : bucket) {
            for (size_t j = 0; j < t.size(); ++j) {
                if (!std::isnan(t[j])) { avg.waveform[j] += t[j]; ++counts[j]; }
            }
        }
        for (size_t j = 0; j < maxLen; ++j) {
            avg.waveform[j] = counts[j] > 0
                ? avg.waveform[j] / counts[j]
                : std::numeric_limits<double>::quiet_NaN();
        }
        avg.n_contributing = static_cast<uint64_t>(bucket.size());
        return avg;
    }

}  // namespace template_generation_detail

inline std::pair<template_io::TemplateFile, template_io::BeatsFile>
buildTemplatesAndBeatsFromPeakResults(const std::vector<output_binfile_data>& peakResults)
{
    using namespace template_generation_detail;

    auto templates = GenerateTemplates(peakResults);

    template_io::TemplateFile tmpl;
    tmpl.bins.resize(peakResults.size());
    for (size_t i = 0; i < peakResults.size(); ++i) {
        bool bad = peakResults[i].bad_segment;
        const TemplateInfo& info = (i < templates.size()) ? templates[i] : TemplateInfo{};
        packBin(tmpl.bins[i], info, bad);
    }

    // SAECG averages -- copy the existing logic from your current
    // buildTemplatesFromPeakResults (the 13 averageOne calls).
    tmpl.saecg.ch1_raw = averageOne(tmpl.bins, [](const template_io::BinTemplates& b) { return b.ch1_raw.ecgTemplate; });
    tmpl.saecg.ch1_squared = averageOne(tmpl.bins, [](const template_io::BinTemplates& b) { return b.ch1_squared.ecgTemplate; });
    tmpl.saecg.ch1_absval = averageOne(tmpl.bins, [](const template_io::BinTemplates& b) { return b.ch1_absval.ecgTemplate; });
    tmpl.saecg.ch1_unfiltered = averageOne(tmpl.bins, [](const template_io::BinTemplates& b) { return b.ch1_unfiltered.ecgTemplate; });
    tmpl.saecg.ch2_raw = averageOne(tmpl.bins, [](const template_io::BinTemplates& b) { return b.ch2_raw.ecgTemplate; });
    tmpl.saecg.ch2_squared = averageOne(tmpl.bins, [](const template_io::BinTemplates& b) { return b.ch2_squared.ecgTemplate; });
    tmpl.saecg.ch2_absval = averageOne(tmpl.bins, [](const template_io::BinTemplates& b) { return b.ch2_absval.ecgTemplate; });
    tmpl.saecg.ch2_unfiltered = averageOne(tmpl.bins, [](const template_io::BinTemplates& b) { return b.ch2_unfiltered.ecgTemplate; });
    tmpl.saecg.ch3_raw = averageOne(tmpl.bins, [](const template_io::BinTemplates& b) { return b.ch3_raw.ecgTemplate; });
    tmpl.saecg.ch3_squared = averageOne(tmpl.bins, [](const template_io::BinTemplates& b) { return b.ch3_squared.ecgTemplate; });
    tmpl.saecg.ch3_absval = averageOne(tmpl.bins, [](const template_io::BinTemplates& b) { return b.ch3_absval.ecgTemplate; });
    tmpl.saecg.ch3_unfiltered = averageOne(tmpl.bins, [](const template_io::BinTemplates& b) { return b.ch3_unfiltered.ecgTemplate; });
    tmpl.saecg.ppg = averageOne(tmpl.bins, [](const template_io::BinTemplates& b) { return b.ppgTemplate; });

    // Beats: ch1 raw kept beats from each TemplateInfo.
    template_io::BeatsFile beats;
    beats.per_bin_beats.resize(templates.size());
    beats.bad_segment.resize(templates.size(), false);
    for (size_t i = 0; i < templates.size(); ++i) {
        beats.bad_segment[i] = (i < peakResults.size())
            ? peakResults[i].bad_segment : false;
        if (!beats.bad_segment[i]) {
            beats.per_bin_beats[i] = std::move(templates[i].kept_beats_ch1_raw);
        }
    }

    return { tmpl, beats };
}