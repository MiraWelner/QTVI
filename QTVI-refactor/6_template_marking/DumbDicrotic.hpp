#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------
namespace dicrotic_detail {

    inline std::vector<double> simpleSmooth(const std::vector<double>& Y, int w) {
        int L = static_cast<int>(Y.size());
        int halfw = w / 2;
        std::vector<double> out(L);
        for (int i = 0; i < L; ++i) {
            double sum = 0; int cnt = 0;
            for (int j = std::max(0, i - halfw); j <= std::min(L - 1, i + halfw); ++j) {
                if (!std::isnan(Y[j])) { sum += Y[j]; ++cnt; }
            }
            out[i] = (cnt > 0) ? sum / cnt : std::nan("");
        }
        return out;
    }

    inline std::vector<double> shearTransform(const std::vector<double>& time,
        const std::vector<double>& slice) {
        if (slice.size() < 2) return slice;
        double m = (slice.front() - slice.back()) / (time.back() - time.front());
        double c = -m * time.front();
        std::vector<double> result(slice.size());
        for (size_t i = 0; i < slice.size(); ++i) {
            result[i] = slice[i] + m * (i + 1) + c;
        }
        return result;
    }

    inline double pointToLineDist(double px, double py,
        double x1, double y1, double x2, double y2) {
        double dx = x2 - x1, dy = y2 - y1;
        double len = std::sqrt(dx * dx + dy * dy);
        if (len < 1e-15) return std::sqrt((px - x1) * (px - x1) + (py - y1) * (py - y1));
        return std::abs(dy * px - dx * py + x2 * y1 - y2 * x1) / len;
    }

    inline bool orthogonalDistThresh(const std::vector<double>& time,
        const std::vector<double>& line,
        const std::vector<double>& curve,
        double thresh) {
        for (size_t i = 1; i < time.size(); ++i) {
            double d = pointToLineDist(time[i], curve[i],
                time[i - 1], line[i - 1],
                time[i], line[i]);
            if (d > thresh) return true;
        }
        return false;
    }

    inline void normalize(std::vector<double>& v) {
        double mn = *std::min_element(v.begin(), v.end());
        double mx = *std::max_element(v.begin(), v.end());
        double r = mx - mn;
        if (r < 1e-15) { std::fill(v.begin(), v.end(), 0.0); return; }
        for (auto& x : v) x = (x - mn) / r;
    }

} // namespace dicrotic_detail

// ---------------------------------------------------------------------------
// dumbDicrotic — returns sample index of dicrotic notch, or -1 if not found
// ---------------------------------------------------------------------------
inline int dumbDicrotic(const std::vector<double>& beat_in) {
    using namespace dicrotic_detail;

    if (beat_in.size() < 4) return -1;

    std::vector<double> beat = simpleSmooth(beat_in, 15);
    if (beat.empty()) return -1;

    // Find peak
    auto it = std::max_element(beat.begin(), beat.end());
    int pmax = static_cast<int>(it - beat.begin());
    int pmin = static_cast<int>(beat.size()) - 1;

    if (pmax == pmin || pmin - pmax < 2) return -1;

    // Regions
    int p_min_dpdt_region = pmax + static_cast<int>(std::round((pmin - pmax) / 3.0));
    int init_EP = pmax + static_cast<int>(std::round((pmin - pmax) * 3.0 / 4.0));

    // Find steepest descent point
    int p_min_dpdt = pmax;
    {
        double minDiff = 1e30;
        for (int i = pmax; i < p_min_dpdt_region && i + 1 < (int)beat.size(); ++i) {
            double d = beat[i + 1] - beat[i];
            if (d < minDiff) { minDiff = d; p_min_dpdt = i; }
        }
    }

    // Half-pressure level
    double p_half = beat[pmax] - (beat[pmax] - beat[p_min_dpdt]) / 2.0;

    // Find SP candidates: points near p_half in the descent
    std::vector<std::pair<double, int>> candidates;
    for (int i = pmax; i <= p_min_dpdt_region && i < (int)beat.size(); ++i) {
        candidates.push_back({ std::abs(beat[i] - p_half), i });
    }
    std::sort(candidates.begin(), candidates.end());

    // Find SP using shear transform test
    int SP = pmax;
    for (const auto& [_, idx] : candidates) {
        std::vector<double> time_v, slice;
        for (int i = idx; i <= init_EP && i < (int)beat.size(); ++i) {
            time_v.push_back(i);
            slice.push_back(beat[i]);
        }
        if (slice.size() < 2) continue;
        auto transform = shearTransform(time_v, slice);
        int below = 0;
        for (size_t i = 0; i < transform.size(); ++i) {
            if (transform[i] < slice[i]) ++below;
        }
        if (below < (int)slice.size() / 2) {
            SP = idx;
            break;
        }
    }

    // Refine EP
    int EP = init_EP;
    while (EP > SP + 1) {
        int len = EP - SP + 1;
        std::vector<double> pressure(len), time_v(len), shearline(len);

        for (int i = 0; i < len; ++i) {
            pressure[i] = beat[SP + i];
            time_v[i] = SP + i;
        }

        double m = (pressure.back() - pressure.front()) / (len > 1 ? (len - 1) : 1);
        for (int i = 0; i < len; ++i) shearline[i] = m * i;

        auto nl = shearline; normalize(nl);
        auto np_ = pressure; normalize(np_);
        auto nt = time_v; normalize(nt);

        if (orthogonalDistThresh(nt, nl, np_, 0.3)) {
            --EP;
        }
        else {
            break;
        }
    }

    // Final shear transform to find notch
    {
        std::vector<double> time_v, slice;
        for (int i = SP; i <= EP && i < (int)beat.size(); ++i) {
            time_v.push_back(i);
            slice.push_back(beat[i]);
        }
        if (slice.size() < 2) return -1;

        auto transform = shearTransform(time_v, slice);

        auto minIt = std::min_element(transform.begin(), transform.end());
        int min_shear = SP + static_cast<int>(minIt - transform.begin());

        int start_relax = min_shear;
        double mx = -1e30;
        for (int i = min_shear; i <= pmin && i < (int)beat.size(); ++i) {
            if (beat[i] > mx) { mx = beat[i]; start_relax = i; }
        }

        int notch = min_shear;
        double mn = 1e30;
        for (int i = min_shear; i <= start_relax && i < (int)beat.size(); ++i) {
            if (beat[i] < mn) { mn = beat[i]; notch = i; }
        }

        return notch;
    }
}