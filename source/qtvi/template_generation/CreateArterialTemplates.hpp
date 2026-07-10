/**
 * @file   CreateArterialTemplates.hpp
 * @brief  Foot-anchored averaged templates for arterial channels
 *         (ABP / ART / ART_PULM), for display as faint background-context
 *         traces in the template-marking viewer.
 *
 *         These are built exactly like the PPG templates (CreatePPGTemplates),
 *         keyed off the PPG feet (ppgMinAmps), but on a different raw signal
 *         (the arterial channel). Because channels are NOT guaranteed to share
 *         a sample rate, the PPG-frame foot indices are rescaled into the
 *         arterial channel's own index space by the length ratio (both signals
 *         cover the same bin time window, so length ratio == rate ratio).
 *
 *         No per-sample std is produced -- these are background context only.
 *         A bin with no arterial signal, no feet, or a bad segment yields an
 *         empty template (=> the viewer simply doesn't draw it).
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "TemplateTypes.hpp"
#include "EnsembleTemplate.hpp"
#include "CreatePPGTemplates.hpp"   // reuse extract_ppg_segment
#include "peak_finding/peakfinding_io.hpp"

 /**
  * @brief  One foot-anchored averaged template per bin for a single arterial
  *         channel, selected by a member pointer into output_binfile_data
  *         (e.g. &output_binfile_data::abpSignal).
  *
  *         Returned vector is parallel to `bins`. Entry i is empty when the
  *         channel is absent for that bin, the bin is bad, or too few feet
  *         land inside the signal.
  */
  // Result of CreateArterialTemplates: per-bin averaged template and its
  // per-sample std, both parallel to the input bins.
struct ArterialTemplatesResult {
    std::vector<std::vector<double>> templates;
    std::vector<std::vector<double>> stds;
};

inline ArterialTemplatesResult CreateArterialTemplates(
    const std::vector<output_binfile_data>& bins,
    std::vector<double> output_binfile_data::* sigMember,
    double std_multiplier)
{
    const size_t n = bins.size();
    ArterialTemplatesResult out;
    out.templates.assign(n, {});
    out.stds.assign(n, {});

    int threads = std::min(8, static_cast<int>(n > 0 ? n : 1));
#pragma omp parallel for schedule(dynamic) num_threads(threads)
    for (int i = 0; i < static_cast<int>(n); ++i) {
        const output_binfile_data& b = bins[i];
        const std::vector<double>& sig = b.*sigMember;

        if (b.bad_segment || sig.empty()
            || b.ppgMinAmps.empty() || b.ppgSignal.empty())
            continue;   // out.templates[i]/stds[i] stay empty

        // Rescale PPG feet into this channel's index space. Both the PPG
        // signal and the arterial signal span the same bin time window, so
        // their length ratio is the sample-rate ratio.
        const double scale = static_cast<double>(sig.size())
            / static_cast<double>(b.ppgSignal.size());

        std::vector<size_t> feet;
        feet.reserve(b.ppgMinAmps.size());
        for (size_t f : b.ppgMinAmps) {
            const long long r = std::llround(static_cast<double>(f) * scale);
            if (r >= 0 && static_cast<size_t>(r) < sig.size())
                feet.push_back(static_cast<size_t>(r));
        }
        if (feet.size() < 2) continue;

        try {
            std::vector<size_t> localPeaks;
            std::vector<double> segment = extract_ppg_segment(sig, feet, localPeaks);
            if (segment.empty() || localPeaks.size() < 2) continue;

            std::vector<double> sd;
            out.templates[i] = EnsembleTemplate(
                segment,
                localPeaks,
                std_multiplier,
                "ppg",     // foot-anchored averaging, same as PPG
                {},        // expand
                nullptr,   // no kept beats
                &sd);      // capture per-sample std for the band
            out.stds[i] = std::move(sd);
        }
        catch (...) {
            out.templates[i].clear();
            out.stds[i].clear();
        }
    }

    return out;
}