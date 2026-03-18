// ============================================================================
// File: filterECG128Hz.hpp
// Bandpass FIR filter for ECG at 128 Hz sampling rate
// ============================================================================
#pragma once

#include <vector>

std::vector<double> filterECG128Hz(const std::vector<double>& x);