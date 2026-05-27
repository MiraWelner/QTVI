/**
 * @file   CreatePPGTemplates.hpp
 * @brief  Create PPG templates for each bin using EnsembleTemplate.
 *         Port of CreatePPGTemplates.m
 *
 *         Each bin stores the full-night PPG signal but ppgMinAmps are
 *         absolute indices.  We extract just the segment spanning the
 *         peaks (with margin) and adjust indices to be local, matching
 *         MATLAB's bins{i}.ppgSeg behaviour.
 *
 *         Returns a struct carrying both the templates and a parallel
 *         vector of per-sample std vectors (same shape). The std is
 *         computed inside EnsembleTemplate over the surviving aligned
 *         beats, and downstream code in GenerateTemplates carries it
 *         through the same AlignWaves shift as the templates.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-26
 */
#pragma once

#include "TemplateTypes.hpp"
#include "EnsembleTemplate.hpp"

struct PPGTemplatesResult {
    vector<vector<double>> templates;   // [bin][sample]
    vector<vector<double>> stds;        // [bin][sample], same shape as templates
};

/**
 * @brief  Extract the minimal segment of the full signal that contains
 *         all peaks, with a margin of the median peak-to-peak interval
 *         on each side. Returns the segment and adjusts peaks to local indices.
 */
static inline vector<double> extract_ppg_segment(
    const vector<double>& fullSignal,
    const vector<size_t>& peaks,
    vector<size_t>& localPeaks)
{
    localPeaks.clear();
    if (peaks.empty() || fullSignal.empty()) return {};

    size_t minIdx = *std::min_element(peaks.begin(), peaks.end());
    size_t maxIdx = *std::max_element(peaks.begin(), peaks.end());

    // Margin: median peak-to-peak interval (or 2000 samples ~ 1 sec fallback)
    size_t margin = 2000;
    if (peaks.size() >= 2) {
        vector<double> pp;
        for (size_t j = 0; j + 1 < peaks.size(); ++j)
            pp.push_back(static_cast<double>(peaks[j + 1] - peaks[j]));
        std::sort(pp.begin(), pp.end());
        margin = static_cast<size_t>(pp[pp.size() / 2]);
    }

    size_t start = (minIdx > margin) ? minIdx - margin : 0;
    size_t end = std::min(maxIdx + margin, fullSignal.size());

    vector<double> segment(fullSignal.begin() + start, fullSignal.begin() + end);

    for (auto idx : peaks) {
        if (idx >= start && idx < end) {
            localPeaks.push_back(idx - start);
        }
    }

    return segment;
}

inline PPGTemplatesResult CreatePPGTemplates(
    const vector<output_binfile_data>& bins,
    double std_multiplier)
{
    size_t n = bins.size();
    PPGTemplatesResult out;
    out.templates.assign(n, {});
    out.stds.assign(n, {});

    int ppg_threads = std::min(8, static_cast<int>(n));
#pragma omp parallel for schedule(dynamic) num_threads(ppg_threads)
    for (int i = 0; i < static_cast<int>(n); ++i) {
        if (bins[i].bad_segment || bins[i].ppgMinAmps.empty()) {
            continue;  // templates[i] and stds[i] stay empty
        }
        try {
            vector<size_t> localPeaks;
            vector<double> segment = extract_ppg_segment(
                bins[i].ppgSignal, bins[i].ppgMinAmps, localPeaks);

            if (segment.empty() || localPeaks.size() < 2) {
                continue;
            }

            // Capture per-sample std alongside the template.
            vector<double> sd;
            out.templates[i] = EnsembleTemplate(
                segment,
                localPeaks,
                std_multiplier,
                "ppg",
                {},        // expand
                nullptr,   // out_kept_beats (not used for PPG)
                &sd);      // out_std
            out.stds[i] = std::move(sd);
        }
        catch (...) {
            out.templates[i] = {};
            out.stds[i] = {};
        }
    }

    // Pad templates to a common length with NaN, as before. Pad stds the
    // same way but with 0 so the band collapses to the line in padded
    // regions (no spurious gray bulge at the tail).
    size_t max_len = 0;
    for (const auto& t : out.templates)
        if (t.size() > max_len) max_len = t.size();

    for (auto& t : out.templates) t.resize(max_len, NaN);
    for (auto& s : out.stds)      s.resize(max_len, 0.0);

    return out;
}