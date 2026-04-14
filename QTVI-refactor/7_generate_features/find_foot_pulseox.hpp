#pragma once
// ═══════════════════════════════════════════════════════════════════════════════
// find_foot_pulseox.hpp — Exact match of find_foot_pulseox.m
// ═══════════════════════════════════════════════════════════════════════════════
#include "ppg_features.hpp"

namespace ppg {

namespace detail {
    inline void findpeaks(const std::vector<double>& data,
                          std::vector<double>& peak_vals,
                          std::vector<int>& peak_locs) {
        peak_vals.clear(); peak_locs.clear();
        for (int i = 1; i < (int)data.size() - 1; i++) {
            if (data[i] > data[i-1] && data[i] > data[i+1]) {
                peak_vals.push_back(data[i]);
                peak_locs.push_back(i);
            }
        }
    }
} // namespace detail

// Returns {value_at_foot, foot_index} (0-based)
// Matches MATLAB: [val,idx] = find_foot_pulseox(data', 0)
// MATLAB passes data as a row vector; we accept column vector.
inline std::pair<double, int> find_foot_pulseox(const std::vector<double>& data) {
    if (data.empty() || data.size() < 4) return {0.0, 0};

    int n = (int)data.size();
    double mx_val = *std::max_element(data.begin(), data.end());

    // MATLAB: data = data - max(data)
    std::vector<double> shifted(n);
    for (int i = 0; i < n; i++) shifted[i] = data[i] - mx_val;

    // MATLAB: [~,m] = max(data)
    int m_idx = 0;
    for (int i = 1; i < n; i++)
        if (shifted[i] > shifted[m_idx]) m_idx = i;

    // MATLAB: diff(data,1,2) — for row vector this is diff along columns
    std::vector<double> ddata(n - 1);
    for (int i = 0; i < n - 1; i++) ddata[i] = shifted[i+1] - shifted[i];

    // MATLAB: [a,b] = findpeaks(diff(data))
    std::vector<double> pvals; std::vector<int> plocs;
    detail::findpeaks(ddata, pvals, plocs);

    // MATLAB: points_of_intrest: peaks where b <= m
    std::vector<std::pair<double, int>> poi;
    for (size_t i = 0; i < plocs.size(); i++)
        if (plocs[i] <= m_idx) poi.push_back({pvals[i], plocs[i]});

    if (poi.empty()) {
        // MATLAB: tmp(:,1) = m; tmp(:,2) = 0;
        // But then the overshoot loop below would use this.
        // Actually MATLAB does: if isempty(tmp), set tmp = [m(i), 0]
        // which means the diffpeak location defaults to m_idx
        auto it = std::max_element(ddata.begin(), ddata.end());
        poi.push_back({*it, (int)(it - ddata.begin())});
    }

    // MATLAB: sorted = sortrows(x, 'descend') — sort by amplitude descending
    std::sort(poi.begin(), poi.end(),
              [](auto& a, auto& b) { return a.first > b.first; });

    // MATLAB: overshoot logic — walk through sorted candidates
    double best_val = 0;
    int best_loc = poi[0].second;
    int overshoot = 0;
    for (auto& [cury, curx] : poi) {
        if (cury > best_val) {
            best_val = cury;
            best_loc = curx;
            overshoot = 0;
        } else {
            overshoot++;
        }
        if (overshoot > 2) break;
    }

    // MATLAB: if diffpeaks(i,2) == 0 → use max(diff(data))
    if (best_loc == 0) {
        auto it = std::max_element(ddata.begin(), ddata.end());
        best_loc = (int)(it - ddata.begin());
    }

    int diffpeak_loc = std::max(1, best_loc);

    // MATLAB: p1_prime = [1 0] - beginPoints → translation so point starts at (1,0)
    // beginPoints = [1, data(1)] → p1_prime = [0, -data(1)]
    double p1_y_offset = -shifted[0];

    // MATLAB: moved = data + p1_prime(:,2)
    std::vector<double> moved(diffpeak_loc + 1);
    for (int i = 0; i <= diffpeak_loc; i++) moved[i] = shifted[i] + p1_y_offset;

    // MATLAB: p2_prime = endPoints + p1_prime
    // endPoints = [diffpeaks(i,2), data(i, diffpeaks(i,2))]
    double p2_x = (double)diffpeak_loc;  // 0-based index used as coordinate
    double p2_y = shifted[diffpeak_loc] + p1_y_offset;

    // MATLAB: theta = atand(p2_prime(i,1)/p2_prime(i,2))
    // then: if theta < 0, theta = theta * -1
    double theta_deg = std::atan2(p2_x, p2_y) * 180.0 / M_PI;
    // MATLAB uses atand(p2_x/p2_y) which is atan(p2_x/p2_y) in radians→degrees
    theta_deg = std::atan(p2_x / p2_y) * 180.0 / M_PI;
    if (theta_deg < 0) theta_deg = -theta_deg;

    double theta_rad = theta_deg * M_PI / 180.0;
    double cos_t = std::cos(theta_rad), sin_t = std::sin(theta_rad);

    // MATLAB: R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]
    // tmp = R * [1:diffpeaks(i,2); moved(i,1:diffpeaks(i,2))]
    // rotated = tmp(1,:)
    // [val,idx] = max(rotated)
    // Note: MATLAB uses 1:diffpeaks(i,2) as x-coords (1-based, integers)
    int best_idx = 0; double best_rotated = -1e18;
    for (int i = 0; i <= diffpeak_loc; i++) {
        // MATLAB x-coord is (i+1) since it uses 1-based indexing
        double rx = cos_t * (double)(i + 1) - sin_t * moved[i];
        if (rx > best_rotated) { best_rotated = rx; best_idx = i; }
    }

    return {data[best_idx], best_idx};
}

} // namespace ppg
