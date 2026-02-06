// ============================================================================
// File: SignalProcessingTypes.h
// Common types and utilities for signal processing
// ============================================================================
#ifndef SIGNAL_PROCESSING_TYPES_H
#define SIGNAL_PROCESSING_TYPES_H

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

// Type aliases for clarity
using std::vector;
using std::string;
using std::map;
using std::pair;
using std::tuple;

// Constants
const double NaN = std::numeric_limits<double>::quiet_NaN();
const double Inf = std::numeric_limits<double>::infinity();

// Helper function to check if value is NaN
inline bool isNaN(double value) {
    return std::isnan(value);
}

// Helper function to get size of vector
template<typename T>
inline size_t numel(const vector<T>& v) {
    return v.size();
}

// Helper function to get length of vector
template<typename T>
inline size_t length(const vector<T>& v) {
    return v.size();
}

#endif // SIGNAL_PROCESSING_TYPES_H