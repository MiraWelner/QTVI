// ============================================================================
// File: pan_tompkin.h
// Pan-Tompkins QRS detection algorithm
// ============================================================================
#pragma once

#include <vector>
#include <string>
#include <cstddef>

using namespace std;

struct PanTompkinResult {
    vector<size_t> qrs_i_raw;   // R-wave sample indices
    vector<double> qrs_amp_raw; // R-wave amplitudes
    int delay;                  // processing delay in samples
};

PanTompkinResult pan_tompkin(const std::vector<double>& ecg_input, double fs, int gr, const std::string& fileID);