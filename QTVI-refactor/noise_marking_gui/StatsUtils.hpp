// ============================================================================
// File: StatsUtils.hpp
// Statistical utility functions (header-only)
// ============================================================================
#pragma once
#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>
#include <limits> 

using std::vector; 
using std::pair;
inline constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
inline constexpr double Inf = std::numeric_limits<double>::infinity();

inline void detrend(std::vector<double>& x) {
    size_t n = x.size();
    if (n < 2) return;

    double sum_x = 0, sum_y = 0, sum_xy = 0, sum_xx = 0;
    for (size_t i = 0; i < n; ++i) {
        sum_x += (double)i;
        sum_y += x[i];
        sum_xy += (double)i * x[i];
        sum_xx += (double)i * (double)i;
    }

    double denom = (n * sum_xx - sum_x * sum_x);
    if (denom == 0) return;

    double slope = (n * sum_xy - sum_x * sum_y) / denom;
    double intercept = (sum_y - slope * sum_x) / (double)n;

    for (size_t i = 0; i < n; ++i)
        x[i] -= (slope * (double)i + intercept);
}

inline double mean(const vector<double>& x) {
    if (x.empty()) return 0.0;
    double sum = 0.0;
    size_t count = 0;
    for (const auto& val : x) {
        if (!std::isnan(val)) { sum += val; count++; }
    }
    return count > 0 ? sum / count : 0.0;
}

inline double std_dev(const vector<double>& x) {
    if (x.size() <= 1) return 0.0;
    double m = mean(x);
    double sum_sq = 0.0;
    size_t count = 0;
    for (const auto& val : x) {
        if (!std::isnan(val)) {
            double diff = val - m;
            sum_sq += diff * diff;
            count++;
        }
    }
    return count > 1 ? std::sqrt(sum_sq / (count - 1)) : 0.0;
}

inline double median(const vector<double>& x) {
    if (x.empty()) return NaN;
    vector<double> sorted;
    for (const auto& val : x)
        if (!std::isnan(val)) sorted.push_back(val);
    if (sorted.empty()) return NaN;
    std::sort(sorted.begin(), sorted.end());
    size_t n = sorted.size();
    return (n % 2 == 0) ? (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0 : sorted[n / 2];
}

inline pair<double, size_t> max_element_index(const vector<double>& x, size_t start, size_t end) {
    if (start >= x.size() || end > x.size() || start >= end)
        return { NaN, 0 };
    double maxVal = -Inf;
    size_t maxIdx = start;
    for (size_t i = start; i < end; ++i) {
        if (!std::isnan(x[i]) && x[i] > maxVal) { maxVal = x[i]; maxIdx = i; }
    }
    return { maxVal, maxIdx - start };
}

inline pair<double, size_t> min_element_index(const vector<double>& x, size_t start, size_t end) {
    if (start >= x.size() || end > x.size() || start >= end)
        return { NaN, 0 };
    double minVal = Inf;
    size_t minIdx = start;
    for (size_t i = start; i < end; ++i) {
        if (!std::isnan(x[i]) && x[i] < minVal) { minVal = x[i]; minIdx = i; }
    }
    return { minVal, minIdx - start };
}

inline vector<double> movmean(const vector<double>& data, size_t window) {
    vector<double> result(data.size());
    size_t back = (window - 1) / 2;
    size_t front = window / 2;
    for (size_t i = 0; i < data.size(); ++i) {
        size_t start = (i >= back) ? i - back : 0;
        size_t end = std::min(i + front + 1, data.size());
        double sum = 0.0;
        size_t count = 0;
        for (size_t j = start; j < end; ++j) {
            if (!std::isnan(data[j])) { sum += data[j]; count++; }
        }
        result[i] = count > 0 ? sum / count : NaN;
    }
    return result;
}

inline vector<double> diff(const vector<double>& x) {
    if (x.size() <= 1) return {};
    vector<double> result(x.size() - 1);
    for (size_t i = 0; i < result.size(); ++i)
        result[i] = x[i + 1] - x[i];
    return result;
}

inline double sum(const vector<double>& x) {
    double s = 0.0;
    for (const auto& val : x)
        if (!std::isnan(val)) s += val;
    return s;
}

inline vector<double> sort(const vector<double>& x) {
    vector<double> sorted = x;
    std::sort(sorted.begin(), sorted.end());
    return sorted;
}

inline vector<size_t> find(const vector<bool>& condition) {
    vector<size_t> indices;
    for (size_t i = 0; i < condition.size(); ++i)
        if (condition[i]) indices.push_back(i);
    return indices;
}