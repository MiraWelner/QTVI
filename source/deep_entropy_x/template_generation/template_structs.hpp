/**
 * @file   template_structs.hpp
 * @brief  The in-memory template-generation structures: ChannelMethodTemplate,
 *         BinTemplates, TemplateFile and BeatsFile.
 */

#pragma once

#include <cstdint>
#include <string>
#include <map>
#include <array>
#include <vector>
#include "template_bank.hpp"
namespace template_structs {

    struct ChannelMethodTemplate {
        std::vector<double> ecgTemplate;
        // Per-sample standard deviation across the beats that contributed
        // to this template. Same length as ecgTemplate, OR empty when not
        // computed (e.g. the squared/absval/unfiltered methods, which the
        // viewer doesn't display).
        std::vector<double> ecg_template_std;
        int    r_col = -1;
        int    median_rr_samples = -1;
    };

    struct BinTemplates {
        ChannelMethodTemplate ch1_raw, ch1_squared, ch1_absval, ch1_unfiltered;
        ChannelMethodTemplate ch2_raw, ch2_squared, ch2_absval, ch2_unfiltered;
        ChannelMethodTemplate ch3_raw, ch3_squared, ch3_absval, ch3_unfiltered;
        std::vector<double>   ppgTemplate;
        // Per-sample std for the PPG template, same length as ppgTemplate
        // (or empty if no PPG / not computed).
        std::vector<double>   ppg_template_std;
        // Foot-anchored averaged arterial templates (ABP / ART / ART_PULM),
        // shown as faint background-context traces in the viewer. Empty when
        // the channel wasn't present in the dataset. No std (background only).
        std::vector<double>   abpTemplate;
        std::vector<double>   artTemplate;
        std::vector<double>   artPulmTemplate;
        // Per-sample std for each arterial template (same length when
        // present, or empty). Written right after each template vector.
        std::vector<double>   abp_template_std;
        std::vector<double>   art_template_std;
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
        int                   ppg_r_col = -1;
        int                   abp_r_col = -1;
        int                   art_r_col = -1;
        int                   art_pulm_r_col = -1;
        bool                  bad_segment = false;


        std::array<tbank::TemplateBank, 3> ecg_bank;
        tbank::TemplateBank ppg_bank;
    };

    struct AveragedTemplate {
        std::vector<double> waveform;
        uint64_t            n_contributing = 0;
    };

    struct TemplateFile {
        std::vector<BinTemplates> bins;
        std::map<int, std::vector<std::array<ChannelMethodTemplate, 3>>> raw_anchors;
        // Indexing: bank_anchors[anchorTag][bin][channel][slot].
        struct BankSlotTemplate {
            std::vector<double> tmpl;
            std::vector<double> tmpl_std;
            uint32_t n_members = 0;
        };
        std::map<int, std::vector<std::array<std::vector<BankSlotTemplate>, 3>>> bank_anchors;
    };

    struct BeatsFile {
        std::vector<bool> bad_segment;
        std::map<std::string, std::vector<std::vector<std::vector<double>>>> per_channel_beats;
        std::map<std::string, std::vector<int>> per_channel_ref_index;
        // Rhythm verdict per kept beat, [bin][beat], keyed like
        // per_channel_beats: 0 = NORMAL, 1 = PVC, 2 = VOTED_PVC. Assigned in
        // alignment.hpp after the slice and before the pruning, and carried
        // here because it cannot be recomputed once the R-peak vector and the
        // kept-beat matrix stop corresponding.
        std::map<std::string, std::vector<std::vector<uint8_t>>> per_channel_rhythm;
    };
}  // namespace template_structs