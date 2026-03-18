// ============================================================================
// File: filterECG256Hz.hpp
// Bandpass FIR filter for ECG at 256 Hz sampling rate
// ============================================================================
#pragma once

#include <vector>

std::vector<double> filterECG256Hz(const std::vector<double>& x);