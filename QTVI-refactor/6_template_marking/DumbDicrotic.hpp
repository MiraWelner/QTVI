#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------
namespace dicrotic_detail {

    // Running-sum smooth - O(n) instead of O(n*w)
    inline void smoothInPlace(const std::vector<double>& Y, int w,
        std::vector<double>& out) {
        int L = static_cast<int>(Y.size());
        out.resize(L);
        if (L == 0) return;
        int halfw = w / 2;

        double sum = 0.0;
        int cnt = 0;

        // Seed: sum of Y[0..halfw]
        for (int j = 0; j <= std::min(halfw, L - 1); ++j) {
            if (!std::isnan(Y[j])) { sum += Y[j]; ++cnt; }
        }
        out[0] = cnt > 0 ? sum / cnt : std::nan("");

        for (int i = 1; i < L; ++i) {
            // Remove element leaving the window (i - halfw - 1)
            int drop = i - halfw - 1;
            if (drop >= 0 && drop < L && !std::isnan(Y[drop])) {
                sum -= Y[drop]; --cnt;
            }
            // Add element entering the window (i + halfw)
            int add = i + halfw;
            if (add < L && !std::isnan(Y[add])) {
                sum += Y[add]; ++cnt;
            }
            out[i] = cnt > 0 ? sum / cnt : std::nan("");
        }
    }

    inline void shearTransformInPlace(int startIdx, int endIdx,
        const std::vector<double>& beat,
        std::vector<double>& result) {
        int len = endIdx - startIdx + 1;
        result.resize(len);
        if (len < 2) { if (len == 1) result[0] = beat[startIdx]; return; }

        double m = (beat[startIdx] - beat[endIdx])
            / static_cast<double>(endIdx - startIdx);
        double c = -m * startIdx;
        for (int i = 0; i < len; ++i) {
            result[i] = beat[startIdx + i] + m * (i + 1) + c;
        }
    }

    inline double pointToLineDist(double px, double py,
        double x1, double y1,
        double x2, double y2) {
        double dx = x2 - x1, dy = y2 - y1;
        double len2 = dx * dx + dy * dy;
        if (len2 < 1e-30)
            return std::sqrt((px - x1) * (px - x1) + (py - y1) * (py - y1));
        return std::abs(dy * px - dx * py + x2 * y1 - y2 * x1) / std::sqrt(len2);
    }

    // Works on raw pointers/offsets to avoid allocation
    inline bool orthogonalDistThresh(const double* time, const double* line,
        const double* curve, int len, double thresh) {
        for (int i = 1; i < len; ++i) {
            double d = pointToLineDist(time[i], curve[i],
                time[i - 1], line[i - 1],
                time[i], line[i]);
            if (d > thresh) return true;
        }
        return false;
    }

    inline void normalizeVec(double* v, int len) {
        double mn = *std::min_element(v, v + len);
        double mx = *std::max_element(v, v + len);
        double r = mx - mn;
        if (r < 1e-15) { std::fill(v, v + len, 0.0); return; }
        double inv = 1.0 / r;
        for (int i = 0; i < len; ++i) v[i] = (v[i] - mn) * inv;
    }

} // namespace dicrotic_detail

inline int dumbDicrotic(const std::vector<double>& beat_in) {
    /**
    * @brief returns sample index of dicrotic notch, or -1 if not found
    */
    using namespace dicrotic_detail;

    int N = static_cast<int>(beat_in.size());
    if (N < 4) return -1;

    // Smooth
    std::vector<double> beat;
    smoothInPlace(beat_in, 15, beat);
    if (beat.empty()) return -1;

    // Find peak
    int pmax = static_cast<int>(
        std::max_element(beat.begin(), beat.end()) - beat.begin());
    int pmin = N - 1;

    if (pmax == pmin || pmin - pmax < 2) return -1;

    // Regions
    int p_min_dpdt_region = pmax + static_cast<int>(
        std::round((pmin - pmax) / 3.0));
    int init_EP = pmax + static_cast<int>(
        std::round((pmin - pmax) * 3.0 / 4.0));

    // Steepest descent (fix 1: +1 to match MATLAB indexing)
    int p_min_dpdt = pmax;
    {
        double minDiff = 1e30;
        int lim = std::min(p_min_dpdt_region, N - 2);
        for (int i = pmax; i <= lim; ++i) {
            double d = beat[i + 1] - beat[i];
            if (d < minDiff) { minDiff = d; p_min_dpdt = i + 1; }
        }
    }

    // Half-pressure level
    double p_half = beat[pmax] - (beat[pmax] - beat[p_min_dpdt]) / 2.0;

    // SP candidates sorted by distance to p_half
    int sliceLen = std::min(p_min_dpdt_region, N - 1) - pmax + 1;
    std::vector<std::pair<double, int>> candidates(sliceLen);
    for (int i = 0; i < sliceLen; ++i) {
        candidates[i] = { std::abs(beat[pmax + i] - p_half), pmax + i };
    }
    std::sort(candidates.begin(), candidates.end());

    // Reusable buffer for shear transform
    std::vector<double> transformBuf;

    // Find SP (fix 2: float division)
    int SP = pmax;
    for (const auto& [dist, idx] : candidates) {
        int len = std::min(init_EP, N - 1) - idx + 1;
        if (len < 2) continue;

        shearTransformInPlace(idx, std::min(init_EP, N - 1), beat, transformBuf);

        int below = 0;
        for (int i = 0; i < len; ++i) {
            if (transformBuf[i] < beat[idx + i]) ++below;
        }
        if (static_cast<double>(below) / len < 0.5) {
            SP = idx;
            break;
        }
    }

    // EP refinement (fix 3: replicate MATLAB polyfit logic)
    // Pre-allocate all buffers at max size to avoid per-iteration allocation
    int maxLen = init_EP - SP + 1;
    std::vector<double> pressure(maxLen), timeVec(maxLen);
    std::vector<double> shearline(maxLen), nl(maxLen), np_(maxLen), nt(maxLen);

    int EP = init_EP;
    while (EP > SP + 1) {
        int len = EP - SP + 1;

        for (int i = 0; i < len; ++i) {
            pressure[i] = beat[SP + i];
            timeVec[i] = static_cast<double>(SP + i);
        }

        // MATLAB's crossed polyfit: x=[EP, SP], y=[pressure[0], pressure[len-1]]
        double pfit_x1 = static_cast<double>(EP);
        double pfit_x2 = static_cast<double>(SP);
        double slope = (pressure[len - 1] - pressure[0]) / (pfit_x2 - pfit_x1);
        double intercept = pressure[0] - slope * pfit_x1;

        double m_used = (std::abs(slope) > 1e-15) ? intercept / slope : 0.0;

        for (int i = 0; i < len; ++i) shearline[i] = m_used * (i + 1);

        // Copy into normalize buffers (avoid repeated allocation)
        std::copy_n(shearline.data(), len, nl.data());
        std::copy_n(pressure.data(), len, np_.data());
        std::copy_n(timeVec.data(), len, nt.data());

        normalizeVec(nl.data(), len);
        normalizeVec(np_.data(), len);
        normalizeVec(nt.data(), len);

        if (orthogonalDistThresh(nt.data(), nl.data(), np_.data(), len, 0.3)) {
            --EP;
        }
        else {
            break;
        }
    }

    // Final: shear transform to find notch
    {
        int len = std::min(EP, N - 1) - SP + 1;
        if (len < 2) return -1;

        shearTransformInPlace(SP, std::min(EP, N - 1), beat, transformBuf);

        int minShearLocal = 0;
        double minVal = transformBuf[0];
        for (int i = 1; i < len; ++i) {
            if (transformBuf[i] < minVal) { minVal = transformBuf[i]; minShearLocal = i; }
        }
        int min_shear = SP + minShearLocal;

        int start_relax = min_shear;
        double mx = -1e30;
        for (int i = min_shear; i <= pmin && i < N; ++i) {
            if (beat[i] > mx) { mx = beat[i]; start_relax = i; }
        }

        int notch = min_shear;
        double mn = 1e30;
        for (int i = min_shear; i <= start_relax && i < N; ++i) {
            if (beat[i] < mn) { mn = beat[i]; notch = i; }
        }

        return notch;
    }
}