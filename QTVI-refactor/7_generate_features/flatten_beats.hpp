#pragma once

#include "data_types.hpp"
#include "common.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

/**
 * @file flatten_beats.hpp
 * @brief Merge per-segment beat features into a single flattened output,
 *        correct indices for segment boundaries, compute inter-beat intervals,
 *        and reclassify sleep stages.
 *
 * Port of flatten_beat_idx.m.
 */

namespace ppg {

    /**
     * @brief Correct beat indices to account for concatenated segment offsets.
     * @param flat     Flattened beats (modified in place).
     * @param segments Annealed segments for bin-index lookup.
     */
    inline void correct_times(FlattenedBeats& flat, const std::vector<AnnealedSegment>& segments) {
        const int n = static_cast<int>(flat.beats.size());
        flat.corrected_time_sec.resize(n, 0.0);
        flat.edge_beat_mask.resize(n, 0);
        std::vector<int> correct_idx_begin(n, 0);

        int beg_idx = 0;
        for (size_t seg = 0; seg < segments.size(); ++seg) {
            for (size_t row = 0; row < segments[seg].ppg_bin_indices.size(); ++row) {
                int bin_start = segments[seg].ppg_bin_indices[row].first;
                int bin_end = segments[seg].ppg_bin_indices[row].second;
                int len = bin_end - bin_start;
                int end_idx = beg_idx + len;

                int first_hit = -1, last_hit = -1;
                for (int i = 0; i < n; ++i) {
                    if (flat.beats[i].idx_begin >= beg_idx && flat.beats[i].idx_begin <= end_idx) {
                        correct_idx_begin[i] = flat.beats[i].idx_begin - beg_idx + bin_start;
                        if (first_hit < 0) first_hit = i;
                        last_hit = i;
                    }
                }
                if (first_hit >= 0) flat.edge_beat_mask[first_hit] = 1;
                if (last_hit >= 0)  flat.edge_beat_mask[last_hit] = 1;

                beg_idx = end_idx;
            }
        }

        // Remove outermost edge marks
        {
            int first = -1, last = -1;
            for (int i = 0; i < n; ++i) {
                if (flat.edge_beat_mask[i] == 1) { if (first < 0) first = i; last = i; }
            }
            if (first >= 0) flat.edge_beat_mask[first] = 0;
            if (last >= 0)  flat.edge_beat_mask[last] = 0;
        }

        double sr = flat.ppg_sample_rate;
        for (int i = 0; i < n; ++i)
            flat.corrected_time_sec[i] = (correct_idx_begin[i] - 1) / sr;
    }

    /**
     * @brief Flatten per-segment beat vectors into a single FlattenedBeats.
     * @param bins      Per-segment beat feature vectors.
     * @param segments  Annealed segments (for offset computation and sample rate).
     * @return FlattenedBeats with corrected indices, inter-beat intervals, and sleep-state classification.
     */
    inline FlattenedBeats flatten_beat_idx(
        const std::vector<std::vector<BeatFeatures>>& bins,
        const std::vector<AnnealedSegment>& segments)
    {
        FlattenedBeats flat;
        if (segments.empty()) return flat;
        flat.ppg_sample_rate = segments[0].ppg_sample_rate;

        // Count total beats and concatenate
        int offset = 0;
        for (size_t seg = 0; seg < bins.size(); ++seg) {
            if (seg > 0) offset += static_cast<int>(segments[seg - 1].po.size()) - 1;
            for (const auto& b : bins[seg]) {
                BeatFeatures bc = b;
                // Shift global indices by segment offset
                auto shift = [&](int& idx) { if (idx >= 0) idx += offset; };
                shift(bc.idx_begin); shift(bc.idx_end); shift(bc.idx_foot);
                shift(bc.idx_pos_slope); shift(bc.idx_systolic);
                shift(bc.idx_neg_slope_b4); shift(bc.idx_dnotch);
                shift(bc.idx_diastolic); shift(bc.idx_neg_slope_after);
                flat.beats.push_back(std::move(bc));
            }
        }

        const int n = static_cast<int>(flat.beats.size());
        const double sr = flat.ppg_sample_rate;

        // Inter-beat intervals
        flat.sec_valley_to_valley.resize(n, 0.0);
        flat.sec_foot_to_foot.resize(n, 0.0);
        for (int i = 1; i < n; ++i) {
            flat.sec_valley_to_valley[i] =
                (flat.beats[i].idx_begin - flat.beats[i - 1].idx_begin) / sr;
            flat.sec_foot_to_foot[i] =
                (flat.beats[i].idx_foot - flat.beats[i - 1].idx_foot) / sr;
        }

        // ── Sleep-state reclassification ──
        // Original: 0=Awake, -1=REM, -2..-5=NREM1-4, NaN=unknown
        // Target:   0=awake-before-sleep, 1=NREM, 2=REM, 3=awake-during-sleep, 4=awake-after-sleep
        std::vector<double> raw_ss(n);
        for (int i = 0; i < n; ++i) raw_ss[i] = flat.beats[i].sleep_stage;

        auto rle = run_length_encode(raw_ss);

        // Find first run of value 0 → before-sleep wake count
        int before_sleep_count = 0;
        for (auto& e : rle) {
            if (e.value == 0.0) { before_sleep_count = e.count; break; }
        }

        // Find last non-zero, non-NaN run → after-sleep count
        int after_sleep_count = 0;
        for (int i = static_cast<int>(rle.size()) - 1; i >= 0; --i) {
            if (rle[i].value != 0.0 && !std::isnan(rle[i].value)) {
                if (i + 1 < static_cast<int>(rle.size()))
                    after_sleep_count = rle[i + 1].count;
                break;
            }
        }

        std::vector<double> ss = raw_ss;
        // Set after-sleep wake
        if (after_sleep_count > 0) {
            for (int i = n - after_sleep_count; i < n; ++i) ss[i] = 4.0;
        }
        // REM → 2
        for (auto& v : ss) if (v == -1.0) v = 2.0;
        // All remaining NREM → 1
        for (auto& v : ss) if (v != 0.0 && v != 2.0 && v != 4.0) v = 1.0;
        // Mid-sleep wake → 3
        for (auto& v : ss) if (v == 0.0) v = 3.0;
        // Correct before-sleep back to 0
        for (int i = 0; i < before_sleep_count && i < n; ++i) ss[i] = 0.0;

        flat.adjusted_sleep_state = ss;

        // Correct times
        correct_times(flat, segments);

        // Time-to-sleep-onset
        flat.sec_to_first_onset.resize(n, 0.0);
        flat.sec_from_last_onset.resize(n, 0.0);
        double onset_time = (before_sleep_count < n) ? flat.corrected_time_sec[before_sleep_count] : 0.0;
        double last_time = (n - after_sleep_count - 1 >= 0) ? flat.corrected_time_sec[n - after_sleep_count - 1] : 0.0;
        for (int i = 0; i < n; ++i) {
            flat.sec_to_first_onset[i] = -(flat.corrected_time_sec[i] - onset_time);
            flat.sec_from_last_onset[i] = flat.corrected_time_sec[i] - last_time;
        }

        // Concatenated PPG
        for (auto& seg : segments) {
            flat.ppg_wout_noise.insert(flat.ppg_wout_noise.end(),
                seg.po.begin(), seg.po.end());
        }

        return flat;
    }

} // namespace ppg