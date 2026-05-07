// ============================================================================
// File: SignalProcessingTypes.h
// Common types and utilities for signal processing
// ============================================================================
#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <stdexcept>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <functional>
#include <tuple>

using std::vector;
using std::string;
using std::map;
using std::pair;
using std::tuple;

constexpr double ECG_SAMPLE_RATE = 1000.0;
constexpr double PPG_SAMPLE_RATE = 1000.0;

inline constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
inline constexpr double Inf = std::numeric_limits<double>::infinity();