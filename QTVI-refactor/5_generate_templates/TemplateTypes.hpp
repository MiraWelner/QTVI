/**
 * @file   TemplateTypes.hpp
 * @brief  Common types for the template generation pipeline.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-26
 */
#pragma once

#include "SignalProcessingTypes.h"
#include "binfile_handling.h"

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <utility>
#include <cstdint>

#ifdef _OPENMP
#include <omp.h>
#endif

 // Per-channel, per-method ECG template results
struct ChannelTemplates {
    vector<double> ecgTemplate_raw;
    vector<double> ecgTemplate_squared;
    vector<double> ecgTemplate_absval;
    vector<double> ecgTemplate_unfiltered;  ///< Template from original signal + raw R-peaks

    double alignment_point_raw = 0.0;
    double alignment_point_squared = 0.0;
    double alignment_point_absval = 0.0;
    double alignment_point_unfiltered = 0.0;

    double avg_r_expand_raw = 0.0;
    double avg_r_expand_squared = 0.0;
    double avg_r_expand_absval = 0.0;
    double avg_r_expand_unfiltered = 0.0;
};

struct TemplateInfo {
    size_t index;
    vector<std::pair<uint64_t, uint64_t>> ppg_bin_indexs;
    vector<std::pair<uint64_t, uint64_t>> ecg_bin_indexs;
    bool bad_segment;

    // ECG templates: 3 channels x 4 methods (raw, squared, absval, unfiltered)
    ChannelTemplates ch1;
    ChannelTemplates ch2;
    ChannelTemplates ch3;

    vector<double> ppgTemplate;
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
};

struct EcgTemplateResult {
    EcgChannelResult ch1;
    EcgChannelResult ch2;
    EcgChannelResult ch3;
};