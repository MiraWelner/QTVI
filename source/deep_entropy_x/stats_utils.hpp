/*
 * StatsUtils.hpp
 * @brief Statistical utility functions including bandpass
 * 
 * @author: Mira Welner
 * @date: 2026-10-01
 * @email: MEW386@pitt.edu
 */

#pragma once
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>
#include <limits> 
#include <stdexcept>

using std::vector;
using std::pair;
inline constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
inline constexpr double Inf = std::numeric_limits<double>::infinity();

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
    const size_t n = data.size();
    vector<double> result(n, NaN);
    if (n == 0 || window == 0) return result;
    const size_t back = (window - 1) / 2;
    const size_t front = window / 2;
    double sum = 0.0;
    size_t count = 0, lo = 0, hi = 0;
    for (size_t i = 0; i < n; ++i) {
        const size_t start = (i >= back) ? i - back : 0;
        const size_t end = std::min(i + front + 1, n);
        while (hi < end) { if (!std::isnan(data[hi])) { sum += data[hi]; ++count; } ++hi; }
        while (lo < start) { if (!std::isnan(data[lo])) { sum -= data[lo]; --count; } ++lo; }
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


struct TukeyFences {
    double fence_lo = std::numeric_limits<double>::quiet_NaN();
    double fence_hi = std::numeric_limits<double>::quiet_NaN();
};

inline vector<bool> keep_within_tukey(const vector<double>& values, double k, TukeyFences* fences = nullptr)
{
    vector<bool> keep(values.size(), true);
    if (fences) *fences = TukeyFences{};
    if (values.size() < 4) return keep;

    vector<double> sorted;
    sorted.reserve(values.size());
    for (double v : values) if (!std::isnan(v)) sorted.push_back(v);
    if (sorted.size() < 4) return keep;
    std::sort(sorted.begin(), sorted.end());

    auto quantile = [&](double q) {
        const double h = q * (sorted.size() - 1);
        const size_t lo = static_cast<size_t>(std::floor(h));
        const size_t hi = static_cast<size_t>(std::ceil(h));
        const double frac = h - lo;
        return sorted[lo] * (1.0 - frac) + sorted[hi] * frac;
        };
    const double q1 = quantile(0.25);
    const double q3 = quantile(0.75);
    const double iqr = q3 - q1;
    if (iqr <= 0.0) return keep;

    const double lo_b = q1 - k * iqr;
    const double hi_b = q3 + k * iqr;
    for (size_t i = 0; i < values.size(); ++i) {
        if (std::isnan(values[i])) continue;
        if (values[i] < lo_b || values[i] > hi_b) keep[i] = false;
    }
    if (fences) { fences->fence_lo = lo_b;  fences->fence_hi = hi_b; }
    return keep;
}

struct PearsonResult {
    //NaN aware Pearson Correlation
    double r = std::numeric_limits<double>::quiet_NaN();
    int    n_overlap = 0;      // finite pairs actually used
    bool   defined() const { return !std::isnan(r); }
};

inline PearsonResult pearson(const std::vector<double>& a, const std::vector<double>& b, int lo = 0, int hi = -1)
{
    PearsonResult out;
    const int n = static_cast<int>(std::min(a.size(), b.size()));
    if (hi < 0 || hi > n) hi = n;
    lo = std::max(0, lo);

    double sa = 0, sb = 0, saa = 0, sbb = 0, sab = 0;
    int cnt = 0;
    for (int k = lo; k < hi; ++k) {
        const double av = a[k], bv = b[k];
        if (std::isnan(av) || std::isnan(bv)) continue;
        sa += av; sb += bv; saa += av * av; sbb += bv * bv; sab += av * bv;
        ++cnt;
    }
    out.n_overlap = cnt;
    if (cnt < 2) return out;                  // r stays NaN

    const double ma = sa / cnt, mb = sb / cnt;
    const double cov = sab / cnt - ma * mb;
    const double va = saa / cnt - ma * ma;
    const double vb = sbb / cnt - mb * mb;
    if (va <= 0.0 || vb <= 0.0) return out;   // flat on one side: r undefined
    out.r = cov / std::sqrt(va * vb);
    return out;
}

struct Biquad {
    //Second order biquad coefficients 
    double b0, b1, b2;
    double a1, a2;
};

namespace bandpass_detail {

    inline std::vector<double> butterworth_poles(int order) {
        std::vector<double> angles;
        for (int k = 0; k < order; k++)
            angles.push_back(M_PI * (2.0 * k + order + 1.0) / (2.0 * order));
        return angles;
    }

} // namespace bandpass_detail

inline std::vector<Biquad> butterworth_lowpass(int order, double cutoff_hz, double sample_rate) {
    //Lowpass design via bilinear transform - cascaded second-order sections
    if (order < 1) throw std::invalid_argument("Order must be >= 1");
    if (cutoff_hz <= 0 || cutoff_hz >= sample_rate / 2.0)
        throw std::invalid_argument("Cutoff must be in (0, Nyquist)");

    double wc = std::tan(M_PI * cutoff_hz / sample_rate);
    auto angles = bandpass_detail::butterworth_poles(order);
    std::vector<Biquad> sections;

    int i = 0;
    while (i < order) {
        Biquad bq;
        if (i + 1 < order) {
            double re = std::cos(angles[i]);
            double A = 1.0, B = -2.0 * re * wc, C = wc * wc;
            double a0 = A + B + C;
            double a1_coeff = -2.0 * A + 2.0 * C;
            double a2_coeff = A - B + C;
            double b0 = C, b1_val = 2.0 * C, b2 = C;

            bq.b0 = b0 / a0; bq.b1 = b1_val / a0; bq.b2 = b2 / a0;
            bq.a1 = a1_coeff / a0; bq.a2 = a2_coeff / a0;
            sections.push_back(bq);
            i += 2;
        }
        else {
            double a0 = 1.0 + wc;
            bq.b0 = wc / a0; bq.b1 = wc / a0; bq.b2 = 0.0;
            bq.a1 = (-1.0 + wc) / a0; bq.a2 = 0.0;
            sections.push_back(bq);
            i += 1;
        }
    }
    return sections;
}

inline std::vector<Biquad> butterworth_highpass(int order, double cutoff_hz, double sample_rate) {
    //Highpass butterworth filter via bilinear transform
    if (order < 1) throw std::invalid_argument("Order must be >= 1");
    if (cutoff_hz <= 0 || cutoff_hz >= sample_rate / 2.0)
        throw std::invalid_argument("Cutoff must be in (0, Nyquist)");

    double wc = std::tan(M_PI * cutoff_hz / sample_rate);
    auto angles = bandpass_detail::butterworth_poles(order);
    std::vector<Biquad> sections;

    int i = 0;
    while (i < order) {
        Biquad bq;
        if (i + 1 < order) {
            double re = std::cos(angles[i]);
            double A = 1.0, B = -2.0 * re * wc, C = wc * wc;
            double a0 = A + B + C;
            double a1_coeff = -2.0 * A + 2.0 * C;
            double a2_coeff = A - B + C;
            double b0 = A, b1_val = -2.0 * A, b2 = A;

            bq.b0 = b0 / a0; bq.b1 = b1_val / a0; bq.b2 = b2 / a0;
            bq.a1 = a1_coeff / a0; bq.a2 = a2_coeff / a0;
            sections.push_back(bq);
            i += 2;
        }
        else {
            double a0 = 1.0 + wc;
            bq.b0 = 1.0 / a0; bq.b1 = -1.0 / a0; bq.b2 = 0.0;
            bq.a1 = (-1.0 + wc) / a0; bq.a2 = 0.0;
            sections.push_back(bq);
            i += 1;
        }
    }
    return sections;
}

inline std::vector<double> apply_biquad(const Biquad& bq, const std::vector<double>& x) {
    //Apply single biquad - Direct Form II Transposed
    size_t n = x.size();
    std::vector<double> y(n);
    double z1 = 0.0, z2 = 0.0;
    for (size_t i = 0; i < n; i++) {
        double in = x[i];
        double out = bq.b0 * in + z1;
        z1 = bq.b1 * in - bq.a1 * out + z2;
        z2 = bq.b2 * in - bq.a2 * out;
        y[i] = out;
    }
    return y;
}

inline std::vector<double> apply_sos(const std::vector<Biquad>& sos, const std::vector<double>& x) {
    //Apply cascaded second-order sections forward
    std::vector<double> y = x;
    for (const auto& bq : sos)
        y = apply_biquad(bq, y);
    return y;
}

inline std::vector<double> filtfilt(const std::vector<Biquad>& sos, const std::vector<double>& x) {
    // filtfilt - zero-phase filtering (forward + reverse) for biquad cascade
    if (x.size() < 4) return x;

    size_t pad_len = 3 * sos.size();
    if (pad_len >= x.size()) pad_len = x.size() - 1;

    std::vector<double> padded(pad_len + x.size() + pad_len);

    for (size_t i = 0; i < pad_len; i++)
        padded[i] = 2.0 * x[0] - x[pad_len - i];
    for (size_t i = 0; i < x.size(); i++)
        padded[pad_len + i] = x[i];
    for (size_t i = 0; i < pad_len; i++)
        padded[pad_len + x.size() + i] = 2.0 * x.back() - x[x.size() - 2 - i];

    std::vector<double> y = apply_sos(sos, padded);
    std::reverse(y.begin(), y.end());
    y = apply_sos(sos, y);
    std::reverse(y.begin(), y.end());

    std::vector<double> result(x.size());
    for (size_t i = 0; i < x.size(); i++)
        result[i] = y[pad_len + i];
    return result;
}

inline std::vector<double> bandpass_filtfilt(int order, double low_hz, double high_hz, double sample_rate, const std::vector<double>& x)
{
    // Convenience: bandpass via cascaded HP + LP with filtfilt
    auto hp = butterworth_highpass(order, low_hz, sample_rate);
    auto lp = butterworth_lowpass(order, high_hz, sample_rate);
    auto y = filtfilt(hp, x);
    return filtfilt(lp, y);
}