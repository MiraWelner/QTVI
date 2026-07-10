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
 *         [uint64 sz][sz x double ecgTemplate_std]   (sz=0 for non-raw methods)
 *         [double alignment_point]
 *         [double avg_r_expand]
 *
 *       Then PPG template + its std:
 *         [uint64 sz][sz x double ppgTemplate]
 *         [uint64 sz][sz x double ppgTemplate_std]
 *
 *       Then arterial background templates (foot-anchored, no std):
 *         [uint64 sz][sz x double abpTemplate]
 *         [uint64 sz][sz x double artTemplate]
 *         [uint64 sz][sz x double artPulmTemplate]
 *       (each sz=0 when that channel was absent in the dataset)
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
#include <vector>
namespace template_io {

    struct ChannelMethodTemplate {
        std::vector<double> ecgTemplate;
        // Per-sample standard deviation across the beats that contributed
        // to this template. Same length as ecgTemplate, OR empty when not
        // computed (e.g. the squared/absval/unfiltered methods, which the
        // viewer doesn't display).
        std::vector<double> ecgTemplate_std;
        double alignment_point = 0.0;
        double avg_r_expand = 0.0;
    };

    struct BinTemplates {
        ChannelMethodTemplate ch1_raw, ch1_squared, ch1_absval, ch1_unfiltered;
        ChannelMethodTemplate ch2_raw, ch2_squared, ch2_absval, ch2_unfiltered;
        ChannelMethodTemplate ch3_raw, ch3_squared, ch3_absval, ch3_unfiltered;
        std::vector<double>   ppgTemplate;
        // Per-sample std for the PPG template, same length as ppgTemplate
        // (or empty if no PPG / not computed).
        std::vector<double>   ppgTemplate_std;
        // Foot-anchored averaged arterial templates (ABP / ART / ART_PULM),
        // shown as faint background-context traces in the viewer. Empty when
        // the channel wasn't present in the dataset. No std (background only).
        std::vector<double>   abpTemplate;
        std::vector<double>   artTemplate;
        std::vector<double>   artPulmTemplate;
        // Per-sample std for each arterial template (same length when
        // present, or empty). Written right after each template vector.
        std::vector<double>   abpTemplate_std;
        std::vector<double>   artTemplate_std;
        std::vector<double>   artPulmTemplate_std;
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
    };

    void write_template_binfile(const std::string& path, const TemplateFile& data);

    TemplateFile read_template_binfile(const std::string& path);

    void write_beats_binfile(const std::string& path, const BeatsFile& data);
    void write_template_csvfile(const std::string& path, const TemplateFile& data, const std::string& fileID, double sampleRateHz = 0.0);

}  // namespace template_io