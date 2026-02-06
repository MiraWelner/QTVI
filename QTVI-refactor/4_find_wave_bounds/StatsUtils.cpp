
// ============================================================================
// File: StatsUtils.cpp
// ============================================================================
#include "StatsUtils.h"

double mean(const vector<double>& x) {
    if (x.empty()) return 0.0;
    double sum = 0.0;
    size_t count = 0;
    for (const auto& val : x) {
        if (!std::isnan(val)) {
            sum += val;
            count++;
        }
    }
    return count > 0 ? sum / count : 0.0;
}

double std(const vector<double>& x) {
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

double median(const vector<double>& x) {
    if (x.empty()) return NaN;

    vector<double> sorted;
    for (const auto& val : x) {
        if (!std::isnan(val)) {
            sorted.push_back(val);
        }
    }

    if (sorted.empty()) return NaN;

    std::sort(sorted.begin(), sorted.end());
    size_t n = sorted.size();

    if (n % 2 == 0) {
        return (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0;
    }
    else {
        return sorted[n / 2];
    }
}

pair<double, size_t> max_element_index(const vector<double>& x, size_t start, size_t end) {
    if (start >= x.size() || end > x.size() || start >= end) {
        return { NaN, 0 };
    }

    double maxVal = -Inf;
    size_t maxIdx = start;

    for (size_t i = start; i < end; ++i) {
        if (!std::isnan(x[i]) && x[i] > maxVal) {
            maxVal = x[i];
            maxIdx = i;
        }
    }

    return { maxVal, maxIdx - start };  // Return relative index
}

pair<double, size_t> min_element_index(const vector<double>& x, size_t start, size_t end) {
    if (start >= x.size() || end > x.size() || start >= end) {
        return { NaN, 0 };
    }

    double minVal = Inf;
    size_t minIdx = start;

    for (size_t i = start; i < end; ++i) {
        if (!std::isnan(x[i]) && x[i] < minVal) {
            minVal = x[i];
            minIdx = i;
        }
    }

    return { minVal, minIdx - start };  // Return relative index
}

vector<double> movmean(const vector<double>& data, size_t window) {
    vector<double> result(data.size());

    for (size_t i = 0; i < data.size(); ++i) {
        size_t start = i >= window / 2 ? i - window / 2 : 0;
        size_t end = std::min(i + window / 2 + 1, data.size());

        double sum = 0.0;
        size_t count = 0;
        for (size_t j = start; j < end; ++j) {
            if (!std::isnan(data[j])) {
                sum += data[j];
                count++;
            }
        }

        result[i] = count > 0 ? sum / count : NaN;
    }

    return result;
}

vector<double> diff(const vector<double>& x) {
    if (x.size() <= 1) return vector<double>();

    vector<double> result(x.size() - 1);
    for (size_t i = 0; i < result.size(); ++i) {
        result[i] = x[i + 1] - x[i];
    }

    return result;
}

double sum(const vector<double>& x) {
    double s = 0.0;
    for (const auto& val : x) {
        if (!std::isnan(val)) {
            s += val;
        }
    }
    return s;
}

vector<double> sort(const vector<double>& x) {
    vector<double> sorted = x;
    std::sort(sorted.begin(), sorted.end());
    return sorted;
}

vector<size_t> find(const vector<bool>& condition) {
    vector<size_t> indices;
    for (size_t i = 0; i < condition.size(); ++i) {
        if (condition[i]) {
            indices.push_back(i);
        }
    }
    return indices;
}