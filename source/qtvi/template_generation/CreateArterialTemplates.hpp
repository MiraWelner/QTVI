/**
 * @file   CreateArterialTemplates.hpp
 * @brief  R-anchored per-beat templates for arterial channels
 *         (ABP / ART / ART_PULM), used as faint background-context traces
 *         in the template-marking viewer.
 *
 *         Under Patch B, arterial channels are sliced the SAME way as PPG
 *         and ECG: [t_R_i - pad, t_R_{i+1} + pad] driven by ch1.raw
 *         (ECG-frame R-peaks). The R sample index maps into this channel's
 *         own sample space via `channelRate / ecgRate`. This gives every
 *         channel identical real-time slice windows -- ECG, PPG, ABP, ART,
 *         ART_PULM all overlay by construction, with R at column
 *         `pad * channelRate` in each.
 *
 *         The old foot-anchored / ppgMinAmps-based logic is gone (and with
 *         it the awkward "rescale PPG feet into arterial index space" hack).
 *         Delegates to build_pulse_template_pair_windowed() in
 *         CreatePPGTemplates.hpp.
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "TemplateTypes.hpp"
#include "CreatePPGTemplates.hpp"   // build_pulse_template_pair_windowed
#include "peak_finding/peakfinding_io.hpp"

struct ArterialTemplatesResult {
    std::vector<std::vector<double>> templates;
    std::vector<std::vector<double>> stds;
    std::vector<std::vector<std::vector<double>>> kept; // [bin][beat][sample]
};

/**
 * @brief  One R-anchored averaged template per bin for a single arterial
 *         channel, selected by a member pointer into output_binfile_data
 *         (e.g. &output_binfile_data::abpSignal).
 *
 * @param bins          Input bins.
 * @param sigMember     Which signal vector to average (e.g. abpSignal).
 * @param ecgRate       ECG sample rate (for R-peak time base).
 * @param channelRate   This arterial channel's sample rate.
 * @param padSeconds    Slice pad on each side of the R-pair (default 0.25).
 */
inline ArterialTemplatesResult CreateArterialTemplates(
    const std::vector<output_binfile_data>& bins,
    std::vector<double> output_binfile_data::* sigMember,
    double ecgRate,
    double channelRate,
    double padSeconds = 0.3)
{
    const size_t n = bins.size();
    ArterialTemplatesResult out;
    out.templates.assign(n, {});
    out.stds.assign(n, {});
    out.kept.assign(n, {});

    if (channelRate <= 0.0) return out;   // channel absent from this dataset

    int threads = std::min(8, static_cast<int>(n > 0 ? n : 1));
#pragma omp parallel for schedule(dynamic) num_threads(threads)
    for (int i = 0; i < static_cast<int>(n); ++i) {
        const output_binfile_data& b = bins[i];
        const std::vector<double>& sig = b.*sigMember;

        if (b.bad_segment || sig.empty() || b.ch1.raw.size() < 2)
            continue;   // out.templates[i]/stds[i]/kept[i] stay empty

        try {
            build_pulse_template_pair_windowed(
                sig, channelRate,
                b.ch1.raw, ecgRate,
                padSeconds,
                out.templates[i], out.stds[i], out.kept[i]);
        }
        catch (...) {
            out.templates[i].clear();
            out.stds[i].clear();
            out.kept[i].clear();
        }
    }

    return out;
}