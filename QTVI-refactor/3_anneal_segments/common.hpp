/**
 * @file   common.hpp
 * @brief  Shared data structures and utilities for the annealing pipeline.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-22
 */
#pragma once

#include <vector>
#include <cstdint>
#include <cmath>
#include <string>
#include <algorithm>
#include <set>

 // ============================================================================
 // Data structures
 // ============================================================================

struct RawData {
    std::vector<double> ppg, ecg1, ecg2, ecg3, sleepStages;
    double ppgSR = 0, ecgSR = 0, scoringEpochSec = 0;
};

struct FinalSegment {
    std::vector<double> ppg;
    std::vector<double> ecg1, ecg2, ecg3;
    std::vector<double> sleep_stages;
    std::vector<std::pair<uint64_t, uint64_t>> ppg_bin_indexs;
    std::vector<std::pair<uint64_t, uint64_t>> ecg_bin_indexs;
    double ppgSampleRate = 0;
    double ecgSampleRate = 0;
    double scoring_epoch_size_sec = 0;
};

/// Per-channel noise markings as time intervals (seconds).
/// ECG noise is only excluded when all 3 channels overlap.
/// PPG noise is excluded independently.
struct NoiseMarkings {
    std::vector<std::pair<double, double>> ecg1;
    std::vector<std::pair<double, double>> ecg2;
    std::vector<std::pair<double, double>> ecg3;
    std::vector<std::pair<double, double>> ppg;
};

struct Exclusion {
    uint64_t idx_start, idx_end;
    int bin_start, bin_end;
};

/// A contiguous good section of signal.
///   dir:  0 = stationary, 1 = move left, 2 = move right
///   flag: 1 = too small, needs merging with a neighbour
struct Section {
    uint64_t begin, end;
    int dir = 0;
    int flag = 0;
};

struct ProjectConfig {
    std::string dataType, binPath, noisePath, annealedPath;
};

struct BinBreaksResult {
    std::vector<uint64_t> bin_breaks;
    int bin_count = 0;
};

// ============================================================================
// MATLAB-compatible rounding
// ============================================================================

/// Banker's rounding (round half to even) to match MATLAB's round().
inline double matlab_round(double x) {
    double r = std::round(x);
    double frac = x - std::floor(x);
    if (std::abs(frac - 0.5) < 1e-12) {
        double f = std::floor(x);
        r = (std::fmod(std::abs(f), 2.0) < 0.5) ? f : f + 1.0;
    }
    return r;
}

/// Convert a time in seconds to a 1-based sample index.
inline uint64_t closest_idx(double time_sec, double sr) {
    double rounded = matlab_round(time_sec * sr);
    if (rounded < 0.0) rounded = 0.0;
    return static_cast<uint64_t>(rounded) + 1;
}