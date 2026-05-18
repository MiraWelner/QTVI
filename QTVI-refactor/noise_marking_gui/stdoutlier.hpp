// ============================================================================
// File: stdoutlier.hpp
// Detect outliers using moving mean and standard deviation of differences
// ============================================================================
#pragma once

#include "StatsUtils.hpp"
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

inline vector<bool> stdoutlier(const vector<double>& data,
    double multiplier,
    size_t mean_window,
    const string& direction,
    bool debug_plot = false)
{
    vector<double> d = diff(data);
    vector<double> x = movmean(d, mean_window);
    double s = std_dev(d);

    vector<double> upper_bound(d.size());
    vector<double> lower_bound(d.size());
    for (size_t i = 0; i < d.size(); ++i) {
        upper_bound[i] = x[i] + s * multiplier;
        lower_bound[i] = x[i] - s * multiplier;
    }

    vector<bool> outlierDiffs(d.size(), false);
    if (direction == "lower") {
        for (size_t i = 0; i < d.size(); ++i)
            outlierDiffs[i] = (d[i] < lower_bound[i]);
    }
    else if (direction == "upper") {
        for (size_t i = 0; i < d.size(); ++i)
            outlierDiffs[i] = (d[i] > upper_bound[i]);
    }
    else {
        for (size_t i = 0; i < d.size(); ++i)
            outlierDiffs[i] = (d[i] > upper_bound[i] || d[i] < lower_bound[i]);
    }

    // Mark both endpoints of each outlier diff
    vector<bool> outliers(data.size(), false);
    for (size_t i = 0; i + 1 < data.size(); ++i) {
        if (i < outlierDiffs.size() && outlierDiffs[i]) {
            outliers[i] = true;
            outliers[i + 1] = true;
        }
    }

    return outliers;
}