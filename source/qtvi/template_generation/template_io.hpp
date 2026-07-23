/**
 * @file   template_io.hpp
 * @brief  Data structures and I/O for the template-generation output
 *         file. Independent of peak_finding's binfile_handling.hpp:
 *         template generation has its own deliverable, written to its
 *         own .bin alongside the wave-markings file.
 *
 *   On-disk layout:
 *
 *     [uint64 n_bins]
 *
 *     For each bin:
 *       12 ECG method blocks (3 channels x 4 methods, fixed order):
 *         ch1_raw, ch1_squared, ch1_absval, ch1_unfiltered,
 *         ch2_raw, ch2_squared, ch2_absval, ch2_unfiltered,
 *         ch3_raw, ch3_squared, ch3_absval, ch3_unfiltered
 *
 *       Each block:
 *         [uint64 sz][sz x double ecgTemplate]
 *         [uint64 sz][sz x double ecg_template_iqr]   (sz=0 for non-raw methods)
 *         [double alignment_point]
 *         [int32  r_col]
 *
 *       Then PPG template + its std:
 *         [uint64 sz][sz x double ppgTemplate]
 *         [uint64 sz][sz x double ppg_template_iqr]
 *
 *       Then arterial background templates, each with its std:
 *         [uint64 sz][sz x double abpTemplate]
 *         [uint64 sz][sz x double abpTemplate_iqr]
 *         [uint64 sz][sz x double artTemplate]
 *         [uint64 sz][sz x double artTemplate_iqr]
 *         [uint64 sz][sz x double artPulmTemplate]
 *         [uint64 sz][sz x double artPulmTemplate_iqr]
 *       (each sz=0 when that channel was absent in the dataset)
 *
 *       Then per-channel beat counts (slices that survived drop rules and
 *       fed each raw-method median). All driven by ch1.raw R-pairs under
 *       Patch B, so they normally read equal; per-channel storage lets a
 *       future filter drop them per-channel:
 *         [uint64 ch1_n_beats_raw]
 *         [uint64 ch2_n_beats_raw]
 *         [uint64 ch3_n_beats_raw]
 *         [uint64 ppg_n_beats]
 *         [int32  ppg_peak_col][int32 ppg_onset_col]
 *
 *       Then bad-segment flag:
 *         [uint8 bad_segment]
 *
 *   Only the ECG "raw" method writes a populated std vector; the other
 *   three methods (squared, absval, unfiltered) write an empty vector
 *   (sz=0, no payload) since the viewer doesn't display them.
 */

#pragma once

#include <cstdint>
#include <string>
#include <map>
#include <vector>
namespace template_io {

    struct ChannelMethodTemplate {
        std::vector<double> ecgTemplate;
        // Per-sample standard deviation across the beats that contributed
        // to this template. Same length as ecgTemplate, OR empty when not
        // computed (e.g. the squared/absval/unfiltered methods, which the
        // viewer doesn't display).
        std::vector<double> ecg_template_iqr;
        double alignment_point = 0.0;
        int    r_col = -1;   // true R column in the template (was avg_r_expand)
    };

    struct BinTemplates {
        ChannelMethodTemplate ch1_raw, ch1_squared, ch1_absval, ch1_unfiltered;
        ChannelMethodTemplate ch2_raw, ch2_squared, ch2_absval, ch2_unfiltered;
        ChannelMethodTemplate ch3_raw, ch3_squared, ch3_absval, ch3_unfiltered;
        std::vector<double>   ppgTemplate;
        // Per-sample std for the PPG template, same length as ppgTemplate
        // (or empty if no PPG / not computed).
        std::vector<double>   ppg_template_iqr;
        // Foot-anchored averaged arterial templates (ABP / ART / ART_PULM),
        // shown as faint background-context traces in the viewer. Empty when
        // the channel wasn't present in the dataset. No std (background only).
        std::vector<double>   abpTemplate;
        std::vector<double>   artTemplate;
        std::vector<double>   artPulmTemplate;
        // Per-sample std for each arterial template (same length when
        // present, or empty). Written right after each template vector.
        std::vector<double>   abpTemplate_iqr;
        std::vector<double>   artTemplate_iqr;
        std::vector<double>   artPulmTemplate_iqr;
        // Per-channel slice counts (post drop-rules) fed to each raw-method
        // median. Under Patch B they're driven by ch1.raw R-pairs, so they
        // normally read equal, but any per-channel drop (short slice, bad
        // signal) would diverge. 0 = unknown (channel absent or bad).
        uint64_t              ch1_n_beats_raw = 0;
        uint64_t              ch2_n_beats_raw = 0;
        uint64_t              ch3_n_beats_raw = 0;
        uint64_t              ppg_n_beats = 0;
        int                   ppg_peak_col = -1;   // construction-time fiducials
        int                   ppg_onset_col = -1;
        bool                  bad_segment = false;

    };

    struct AveragedTemplate {
        std::vector<double> waveform;
        uint64_t            n_contributing = 0;
    };

    struct TemplateFile {
        std::vector<BinTemplates> bins;
    };

    struct BeatsFile {
        std::vector<std::vector<std::vector<double>>> per_bin_beats;
        std::vector<bool> bad_segment;
        std::map<std::string, std::vector<std::vector<std::vector<double>>>> per_channel_beats;
    };

    void write_template_binfile(const std::string& path, const TemplateFile& data);
    void write_snips_csv(const std::string& path, const BeatsFile& beats);

    TemplateFile read_template_binfile(const std::string& path);
}  // namespace template_io