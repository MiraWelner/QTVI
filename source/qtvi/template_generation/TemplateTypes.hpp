/**
 * @file   TemplateTypes.hpp
 * @brief  Common types for the template generation pipeline.
 *
 *         std vectors: we only carry per-sample std for the ECG "raw"
 *         method (the one the viewer displays) and for PPG. The other
 *         three ECG methods (squared, absval, unfiltered) don't get
 *         std fields -- the viewer doesn't display them, so computing
 *         and storing std for them would just inflate disk usage.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-26
 */
#pragma once
#include <map>
#include "template_io.hpp"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <utility>
#include <cstdint>
#ifdef _OPENMP
#include <omp.h>
#endif

 // Per-channel sample rates (Hz), threaded through the template-generation
 // pipeline. Each channel is at its own rate; the slicer converts between
 // them via the ratio channelRate / ecgRate. A rate of 0 means the channel
 // is absent from this dataset (skip it).
struct SignalRates {
    double ecg = 0.0;
    double ppg = 0.0;
    double abp = 0.0;
    double art = 0.0;
    double artPulm = 0.0;
};

// Per-channel, per-method ECG template results
struct ChannelTemplates {
    vector<double> ecgTemplate_raw;
    vector<double> ecgTemplate_raw_std;       // per-sample std of the beats
    // contributing to ecgTemplate_raw.
    // Same length as ecgTemplate_raw,
    // or empty if not computed.
    vector<double> ecgTemplate_squared;
    vector<double> ecgTemplate_absval;
    vector<double> ecgTemplate_unfiltered;
    double alignment_point_raw = 0.0;
    double alignment_point_squared = 0.0;
    double alignment_point_absval = 0.0;
    double alignment_point_unfiltered = 0.0;
    double avg_r_expand_raw = 0.0;
    double avg_r_expand_squared = 0.0;
    double avg_r_expand_absval = 0.0;
    double avg_r_expand_unfiltered = 0.0;
    // Slice count fed to the raw-method median for this bin/channel.
    // Only tracked for the raw method since that's what the viewer shows.
    size_t n_beats_raw = 0;
};

struct TemplateInfo {
    ChannelTemplates ch1;
    ChannelTemplates ch2;
    ChannelTemplates ch3;
    std::vector<double> ppgTemplate;
    std::vector<double> ppgTemplate_std;      // per-sample std of the beats
    // contributing to ppgTemplate
    // (post AlignWaves shift).
    // Slice count that fed the PPG median (post drop rules).
    size_t ppg_n_beats = 0;
    // Surviving beats from ch1 raw method (each entry is one beat's
    // samples, all of equal length, possibly with NaN tails). Only
    // populated when capture_beats was requested for ch1 in
    // CreateEcgTemplates.
    std::vector<std::vector<double>> kept_beats_ch1_raw;
    // Retained per-channel beats for the snips CSV. Key is the channel label
    // ("CH1"/"CH2"/"CH3"/"PPG"); value is [beat][sample] for this bin.
    std::map<std::string, std::vector<std::vector<double>>> kept_beats_by_channel;
};

struct AlignWavesResult {
    vector<vector<double>> alignedWaves;
    vector<int> move_dist;
};

struct CombineResult {
    vector<size_t> bin_numbers;
    vector<vector<double>> bin_templates;
    vector<size_t> foot_locations;
};

struct FootResult {
    vector<double> val;
    vector<size_t> idx;
};

struct EcgChannelResult {
    vector<vector<double>> ecgTemplates_raw;
    vector<vector<double>> ecgTemplates_raw_std;   // parallel to ecgTemplates_raw
    vector<vector<double>> ecgTemplates_squared;
    vector<vector<double>> ecgTemplates_absval;
    vector<vector<double>> ecgTemplates_unfiltered;

    vector<double> ppg_alignment_point_raw;
    vector<double> ppg_alignment_point_squared;
    vector<double> ppg_alignment_point_absval;
    vector<double> ppg_alignment_point_unfiltered;

    vector<double> avg_r_expand_raw;
    vector<double> avg_r_expand_squared;
    vector<double> avg_r_expand_absval;
    vector<double> avg_r_expand_unfiltered;

    // Per-bin slice count that fed the raw-method median (post drop rules).
    // Only tracked for the raw method since that's what the viewer displays.
    vector<size_t> n_beats_raw;

    // Per-bin captured surviving beats from the "raw" method.
    // Only populated when the caller asks for it (currently: ch1 only).
    // [bin][beat][sample]
    vector<vector<vector<double>>> kept_beats_raw;
};

struct EcgTemplateResult {
    EcgChannelResult ch1;
    EcgChannelResult ch2;
    EcgChannelResult ch3;
};