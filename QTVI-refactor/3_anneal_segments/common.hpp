#pragma once

#include <vector>
#include <cstdint>
#include <cmath>
#include <string>
#include <algorithm>
#include <set>

// ============================================================================
// DATA STRUCTURES
// ============================================================================

struct RawData {
    std::vector<double> ppg, ecg1, ecg2, ecg3, sleepStages;
    double ppgSR, ecgSR, scoringEpochSec;
};

struct FinalSegment {
    std::vector<double> ppg;
    std::vector<double> ecg1, ecg2, ecg3;
    std::vector<double> sleep_stages;
    std::vector<std::pair<uint64_t, uint64_t>> ppg_bin_indexs;
    std::vector<std::pair<uint64_t, uint64_t>> ecg_bin_indexs;
    double ppgSampleRate;
    double ecgSampleRate;
    double scoring_epoch_size_sec;
};

/**
 * @brief Per-channel noise markings as time intervals (seconds).
 *
 * ECG noise is only excluded when all 3 ECG channels have overlapping
 * noise marks for a given time region. PPG noise is excluded independently.
 */
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

struct Section {
    uint64_t begin, end;
    int dir;   // 1=left, 2=right, 0=none
    int flag;
};

struct ProjectConfig {
    std::string dataType, binPath, noisePath, annealedPath;
};

struct BinBreaksResult {
    std::vector<uint64_t> bin_breaks;
    int bin_count;
};

// ============================================================================
// MATLAB-compatible round and closest_idx
// ============================================================================

inline double matlab_round(double x) {
    double r = std::round(x);
    double diff = x - std::floor(x);
    if (std::abs(diff - 0.5) < 1e-12) {
        double f = std::floor(x);
        if (std::fmod(std::abs(f), 2.0) < 0.5)
            r = f;
        else
            r = f + 1.0;
    }
    return r;
}

inline uint64_t closest_idx(double target_time, double sr) {
    double raw = target_time * sr;
    double rounded = matlab_round(raw);
    if (rounded < 0.0) rounded = 0.0;
    return static_cast<uint64_t>(rounded) + 1;
}