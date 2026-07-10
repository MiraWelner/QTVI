/**
* @file   build_templates.hpp
* @brief  Run the template-generation pipeline on peak-finding output
*         and produce a TemplateFile ready to write to disk.
*
*         Carries per-sample std for the ECG raw method (the one the
*         viewer displays) and for PPG, from the in-memory TemplateInfo
*         into the on-disk BinTemplates layout. The other three ECG
*         methods don't have std computed and write empty std vectors.
*/

#pragma once

#include <cmath>
#include <vector>

#include "template_io.hpp"
#include "GenerateTemplates.hpp"
#include "TemplateTypes.hpp"
#include "template_generation/CreateArterialTemplates.hpp"
#include "peak_finding/peakfinding_io.hpp"

namespace template_generation_detail {

    // Copy one template (+ optional std) into the on-disk method block.
    // The non-raw methods pass an empty tmpl_std and the on-disk std
    // field stays sz=0 (no payload). Both writer and reader handle that
    // uniformly, so there's only ever one code path.
    inline void copyMethod(template_io::ChannelMethodTemplate& dst,
        const std::vector<double>& tmpl,
        const std::vector<double>& tmpl_std,
        double alignment, double rExpand)
    {
        dst.ecgTemplate = tmpl;
        dst.ecgTemplate_std = tmpl_std;
        dst.alignment_point = std::isnan(alignment) ? 0.0 : alignment;
        dst.avg_r_expand = rExpand;
    }

    inline void packBin(template_io::BinTemplates& bt,
        const TemplateInfo& info, bool bad_segment)
    {
        bt.bad_segment = bad_segment;
        if (bad_segment) return;

        // Only the raw methods carry std; the other three pass an empty
        // vector via the default-constructed std::vector<double>{}.
        const std::vector<double> noStd;

        copyMethod(bt.ch1_raw, info.ch1.ecgTemplate_raw,
            info.ch1.ecgTemplate_raw_std,
            info.ch1.alignment_point_raw, info.ch1.avg_r_expand_raw);
        copyMethod(bt.ch1_squared, info.ch1.ecgTemplate_squared, noStd,
            info.ch1.alignment_point_squared, info.ch1.avg_r_expand_squared);
        copyMethod(bt.ch1_absval, info.ch1.ecgTemplate_absval, noStd,
            info.ch1.alignment_point_absval, info.ch1.avg_r_expand_absval);
        copyMethod(bt.ch1_unfiltered, info.ch1.ecgTemplate_unfiltered, noStd,
            info.ch1.alignment_point_unfiltered, info.ch1.avg_r_expand_unfiltered);

        copyMethod(bt.ch2_raw, info.ch2.ecgTemplate_raw,
            info.ch2.ecgTemplate_raw_std,
            info.ch2.alignment_point_raw, info.ch2.avg_r_expand_raw);
        copyMethod(bt.ch2_squared, info.ch2.ecgTemplate_squared, noStd,
            info.ch2.alignment_point_squared, info.ch2.avg_r_expand_squared);
        copyMethod(bt.ch2_absval, info.ch2.ecgTemplate_absval, noStd,
            info.ch2.alignment_point_absval, info.ch2.avg_r_expand_absval);
        copyMethod(bt.ch2_unfiltered, info.ch2.ecgTemplate_unfiltered, noStd,
            info.ch2.alignment_point_unfiltered, info.ch2.avg_r_expand_unfiltered);

        copyMethod(bt.ch3_raw, info.ch3.ecgTemplate_raw,
            info.ch3.ecgTemplate_raw_std,
            info.ch3.alignment_point_raw, info.ch3.avg_r_expand_raw);
        copyMethod(bt.ch3_squared, info.ch3.ecgTemplate_squared, noStd,
            info.ch3.alignment_point_squared, info.ch3.avg_r_expand_squared);
        copyMethod(bt.ch3_absval, info.ch3.ecgTemplate_absval, noStd,
            info.ch3.alignment_point_absval, info.ch3.avg_r_expand_absval);
        copyMethod(bt.ch3_unfiltered, info.ch3.ecgTemplate_unfiltered, noStd,
            info.ch3.alignment_point_unfiltered, info.ch3.avg_r_expand_unfiltered);

        bt.ppgTemplate = info.ppgTemplate;
        bt.ppgTemplate_std = info.ppgTemplate_std;
    }

    // FAST pack: raw + unfiltered ECG blocks + PPG. Leaves the squared and
    // absval blocks default-empty for packBinSlow.
    inline void packBinFast(template_io::BinTemplates& bt,
        const TemplateInfo& info, bool bad_segment)
    {
        bt.bad_segment = bad_segment;
        if (bad_segment) return;

        copyMethod(bt.ch1_raw, info.ch1.ecgTemplate_raw,
            info.ch1.ecgTemplate_raw_std,
            info.ch1.alignment_point_raw, info.ch1.avg_r_expand_raw);
        copyMethod(bt.ch1_unfiltered, info.ch1.ecgTemplate_unfiltered, {},
            info.ch1.alignment_point_unfiltered, info.ch1.avg_r_expand_unfiltered);

        copyMethod(bt.ch2_raw, info.ch2.ecgTemplate_raw,
            info.ch2.ecgTemplate_raw_std,
            info.ch2.alignment_point_raw, info.ch2.avg_r_expand_raw);
        copyMethod(bt.ch2_unfiltered, info.ch2.ecgTemplate_unfiltered, {},
            info.ch2.alignment_point_unfiltered, info.ch2.avg_r_expand_unfiltered);

        copyMethod(bt.ch3_raw, info.ch3.ecgTemplate_raw,
            info.ch3.ecgTemplate_raw_std,
            info.ch3.alignment_point_raw, info.ch3.avg_r_expand_raw);
        copyMethod(bt.ch3_unfiltered, info.ch3.ecgTemplate_unfiltered, {},
            info.ch3.alignment_point_unfiltered, info.ch3.avg_r_expand_unfiltered);

        bt.ppgTemplate = info.ppgTemplate;
        bt.ppgTemplate_std = info.ppgTemplate_std;
    }

    // SLOW pack: squared + absval blocks onto an already fast-packed bin.
    inline void packBinSlow(template_io::BinTemplates& bt,
        const TemplateInfo& info)
    {
        if (bt.bad_segment) return;
        const std::vector<double> noStd;
        copyMethod(bt.ch1_squared, info.ch1.ecgTemplate_squared, noStd,
            info.ch1.alignment_point_squared, info.ch1.avg_r_expand_squared);
        copyMethod(bt.ch1_absval, info.ch1.ecgTemplate_absval, noStd,
            info.ch1.alignment_point_absval, info.ch1.avg_r_expand_absval);
        copyMethod(bt.ch2_squared, info.ch2.ecgTemplate_squared, noStd,
            info.ch2.alignment_point_squared, info.ch2.avg_r_expand_squared);
        copyMethod(bt.ch2_absval, info.ch2.ecgTemplate_absval, noStd,
            info.ch2.alignment_point_absval, info.ch2.avg_r_expand_absval);
        copyMethod(bt.ch3_squared, info.ch3.ecgTemplate_squared, noStd,
            info.ch3.alignment_point_squared, info.ch3.avg_r_expand_squared);
        copyMethod(bt.ch3_absval, info.ch3.ecgTemplate_absval, noStd,
            info.ch3.alignment_point_absval, info.ch3.avg_r_expand_absval);
    }

}  // namespace template_generation_detail

// Carries the in-memory state from the fast build to the slow merge.
struct FastTemplateBuild {
    template_io::TemplateFile tmpl;   // raw/unfiltered/ppg blocks + their SAECG
    template_io::BeatsFile beats;     // ch1 raw kept beats (final)
    std::vector<TemplateInfo> info;   // retained so mergeTemplatesSlow can pack squared/absval
};

// FAST build: produces a TemplateFile the viewer can open immediately
// (raw/unfiltered/ppg + their SAECG averages) and the final beats file.
// The squared/absval per-bin blocks and their SAECG entries stay empty
// until mergeTemplatesSlow runs.
inline FastTemplateBuild
buildTemplatesAndBeatsFast(const std::vector<output_binfile_data>& peakResults)
{
    using namespace template_generation_detail;

    FastTemplateBuild out;
    out.info = GenerateTemplatesFast(peakResults);

    out.tmpl.bins.resize(peakResults.size());
    for (size_t i = 0; i < peakResults.size(); ++i) {
        bool bad = peakResults[i].bad_segment;
        const TemplateInfo& info = (i < out.info.size()) ? out.info[i] : TemplateInfo{};
        packBinFast(out.tmpl.bins[i], info, bad);
    }

    // Arterial background-context templates (ABP / ART / ART_PULM),
    // foot-anchored on the PPG feet, present-only. Built directly from the
    // peak results (which carry the raw arterial signals + ppgMinAmps) and
    // stored into each good bin. They survive mergeTemplatesSlow untouched.
    {
        constexpr double kArtStdMult = 2.5;
        auto abp = CreateArterialTemplates(
            peakResults, &output_binfile_data::abpSignal, kArtStdMult);
        auto art = CreateArterialTemplates(
            peakResults, &output_binfile_data::artSignal, kArtStdMult);
        auto artp = CreateArterialTemplates(
            peakResults, &output_binfile_data::artPulmSignal, kArtStdMult);
        for (size_t i = 0; i < out.tmpl.bins.size(); ++i) {
            if (out.tmpl.bins[i].bad_segment) continue;
            if (i < abp.templates.size()) {
                out.tmpl.bins[i].abpTemplate = std::move(abp.templates[i]);
                out.tmpl.bins[i].abpTemplate_std = std::move(abp.stds[i]);
            }
            if (i < art.templates.size()) {
                out.tmpl.bins[i].artTemplate = std::move(art.templates[i]);
                out.tmpl.bins[i].artTemplate_std = std::move(art.stds[i]);
            }
            if (i < artp.templates.size()) {
                out.tmpl.bins[i].artPulmTemplate = std::move(artp.templates[i]);
                out.tmpl.bins[i].artPulmTemplate_std = std::move(artp.stds[i]);
            }
        }
    }

    // Beats: ch1 raw kept beats (moves them out of out.info; the slow merge
    // doesn't read kept beats).
    out.beats.per_bin_beats.resize(out.info.size());
    out.beats.bad_segment.resize(out.info.size(), false);
    for (size_t i = 0; i < out.info.size(); ++i) {
        out.beats.bad_segment[i] = (i < peakResults.size())
            ? peakResults[i].bad_segment : false;
        if (!out.beats.bad_segment[i])
            out.beats.per_bin_beats[i] = std::move(out.info[i].kept_beats_ch1_raw);
    }

    return out;
}

// SLOW merge: fills the squared/absval per-bin blocks and their SAECG
// averages into a TemplateFile produced by buildTemplatesAndBeatsFast.
// peakResults must carry the squared/absval R-peaks + preprocessed signals
// (run augment_ecg_ppg_pairs_sqabs first, or load them from wave_markings).
inline void mergeTemplatesSlow(const std::vector<output_binfile_data>& peakResults,
    template_io::TemplateFile& tmpl,
    std::vector<TemplateInfo>& info)
{
    using namespace template_generation_detail;

    AugmentTemplatesSlow(peakResults, info);

    for (size_t i = 0; i < tmpl.bins.size(); ++i) {
        const TemplateInfo& bi = (i < info.size()) ? info[i] : TemplateInfo{};
        packBinSlow(tmpl.bins[i], bi);
    }
}

inline std::pair<template_io::TemplateFile, template_io::BeatsFile>
buildTemplatesAndBeatsFromPeakResults(const std::vector<output_binfile_data>& peakResults)
{
    FastTemplateBuild fast = buildTemplatesAndBeatsFast(peakResults);
    mergeTemplatesSlow(peakResults, fast.tmpl, fast.info);
    return { std::move(fast.tmpl), std::move(fast.beats) };
}