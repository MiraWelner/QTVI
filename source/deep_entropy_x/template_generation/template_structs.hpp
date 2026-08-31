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
#include "template_io.hpp"
#include "template_morphology_grouping/bin_pipeline.hpp"

#include <map>
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
    vector<double> ecgTemplate_raw_iqr;       // per-sample std of the beats
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
    // True R column in the template (from alignment's r_aligned_col). This is
    // the detected-R fiducial the template was built around -- used directly
    // as the R marker, replacing the old avg_r_expand positioning constant.
    int r_col_raw = -1;
    int r_col_squared = -1;
    int r_col_absval = -1;
    int r_col_unfiltered = -1;
    // Slice count fed to the raw-method median for this bin/channel.
    // Only tracked for the raw method since that's what the viewer shows.
    size_t n_beats_raw = 0;
};

struct TemplateInfo {
    ChannelTemplates ch1;
    ChannelTemplates ch2;
    ChannelTemplates ch3;
    std::vector<double> ppgTemplate;
    std::vector<double> ppg_template_iqr;      // per-sample std of the beats
    // contributing to ppgTemplate
    // (post AlignWaves shift).
    // Slice count that fed the PPG median (post drop rules).
    size_t ppg_n_beats = 0;
    // Deterministic PPG fiducials computed at construction from the real
    // R-pair interval: peak = max in [R1,R2], foot = min in [R1,peak].
    // -1 when no PPG for this bin.
    int ppg_peak_col = -1;
    int ppg_onset_col = -1;
    // Surviving beats from ch1 raw method (each entry is one beat's
    // samples, all of equal length, possibly with NaN tails). Only
    // populated when capture_beats was requested for ch1 in
    // CreateEcgTemplates.
    std::vector<std::vector<double>> kept_beats_ch1_raw;
    std::map<std::string, int> ref_index_by_channel;   // channel -> ref beat idx
    // Retained per-channel beats for the snips CSV. Key is the channel label
    // ("CH1"/"CH2"/"CH3"/"PPG"); value is [beat][sample] for this bin.
    std::map<std::string, std::vector<std::vector<double>>> kept_beats_by_channel;
    // Rhythm verdict per kept beat, parallel to kept_beats_by_channel[ch]:
    // 0 = NORMAL, 1 = PVC (premature), 2 = VOTED_PVC (5-of-8 vote).
    // Assigned in alignment.hpp after the slice and before any pruning, and
    // carried here because it cannot be recomputed downstream: the R-peak
    // vector and the kept-beat matrix stop corresponding the moment the Tukey
    // passes run.
    std::map<std::string, std::vector<uint8_t>> kept_rhythm_by_channel;
    std::map<std::string, bin_pipeline::ChannelOutput> bank_by_channel;
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
    vector<vector<double>> ecgTemplates_raw_iqr;   // parallel to ecgTemplates_raw
    vector<vector<double>> ecgTemplates_squared;
    vector<vector<double>> ecgTemplates_absval;
    vector<vector<double>> ecgTemplates_unfiltered;
    std::vector<int> ref_index_raw;

    vector<double> ppg_alignment_point_raw;
    vector<double> ppg_alignment_point_squared;
    vector<double> ppg_alignment_point_absval;
    vector<double> ppg_alignment_point_unfiltered;

    vector<int> r_col_raw;
    vector<int> r_col_squared;
    vector<int> r_col_absval;
    vector<int> r_col_unfiltered;

    vector<size_t> n_beats_raw;//the viewer displays the number of beats contributing to template for each channel

    vector<vector<vector<double>>> kept_beats_raw;
    vector<bin_pipeline::ChannelOutput> bank_out_raw;

    // Rhythm verdict per kept beat, [bin][beat]: 0 NORMAL, 1 PVC, 2 VOTED_PVC.
    vector<vector<uint8_t>> kept_rhythm_raw;

    // Per-bin, per-beat vertical DC leveling shifts (two-stage TP/PQ),
    // indexed [bin][beat]. Written to the beat-move log post-loop.
    vector<vector<double>> tp_shift_raw;
    vector<vector<double>> pq_shift_raw;
};

struct EcgTemplateResult {
    EcgChannelResult ch1;
    EcgChannelResult ch2;
    EcgChannelResult ch3;
};