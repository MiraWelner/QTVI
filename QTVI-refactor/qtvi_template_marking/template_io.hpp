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
 *         [double alignment_point]
 *         [double avg_r_expand]
 *
 *       Then PPG template:
 *         [uint64 sz][sz x double ppgTemplate]
 *
 *       Then bad-segment flag:
 *         [uint8 bad_segment]
 *
 *     SAECG tail (recording-wide averages, written once):
 *       13 averaged waveforms in the same order as per-bin blocks
 *       (12 ECG methods + PPG). Each:
 *         [uint64 sz][sz x double avg]
 *         [uint64 n_contributing_bins]
 */

#pragma once

#include <cstdint>
#include <string>
#include <vector>
namespace template_io {

    struct ChannelMethodTemplate {
        std::vector<double> ecgTemplate;
        double alignment_point = 0.0;
        double avg_r_expand = 0.0;
    };

    struct BinTemplates {
        ChannelMethodTemplate ch1_raw, ch1_squared, ch1_absval, ch1_unfiltered;
        ChannelMethodTemplate ch2_raw, ch2_squared, ch2_absval, ch2_unfiltered;
        ChannelMethodTemplate ch3_raw, ch3_squared, ch3_absval, ch3_unfiltered;
        std::vector<double>   ppgTemplate;
        bool                  bad_segment = false;
    };

    struct AveragedTemplate {
        std::vector<double> waveform;
        uint64_t            n_contributing = 0;
    };

    struct SAECG {
        AveragedTemplate ch1_raw, ch1_squared, ch1_absval, ch1_unfiltered;
        AveragedTemplate ch2_raw, ch2_squared, ch2_absval, ch2_unfiltered;
        AveragedTemplate ch3_raw, ch3_squared, ch3_absval, ch3_unfiltered;
        AveragedTemplate ppg;
    };

    struct TemplateFile {
        std::vector<BinTemplates> bins;
        SAECG                     saecg;
    };

    struct BeatsFile {
        std::vector<std::vector<std::vector<double>>> per_bin_beats;
        std::vector<bool> bad_segment;
    };

    void write_template_binfile(const std::string& path, const TemplateFile& data);
    TemplateFile read_template_binfile(const std::string& path);

    void write_beats_binfile(const std::string& path, const BeatsFile& data);

}  // namespace template_io