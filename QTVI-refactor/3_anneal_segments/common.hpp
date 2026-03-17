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

double matlab_round(double x);
uint64_t closest_idx(double target_time, double sr);
