/**
 * @file   beat_log.hpp
 * @brief  Per-beat log across the markable channels (ECG1, ECG2, ECG3,
 *         PPG, ABP). One row per beat; each channel contributes an (x, y)
 *         pair (time, amplitude) plus the blanking/threshold used.
 *
 *         Column order written to CSV:
 *           beat,
 *           ecg1_x, ecg1_y, ecg2_x, ecg2_y, ecg3_x, ecg3_y,
 *           ppg_x,  ppg_y,  abp_x,  abp_y,
 *           blanking_ecg1, threshold_ecg1, ... blanking_abp, threshold_abp
 *
 *         Detected peaks are buffered per channel in a time-keyed map
 *         (m_pending) as the GUI finds them, keyed by global time so order
 *         is independent of the order windows are viewed and duplicates of
 *         the same beat collapse. flushPending() merges the buffer into the
 *         fixed beat table (m_beats) -- a beat with the same time overwrites,
 *         new ones insert in time order -- then empties the buffer. writeCsv
 *         dumps whatever the table currently holds; rows past the last beat
 *         keep x = 0 as the unfilled-tail marker.
 *
 *         Intentionally Qt-free so it can be unit-tested / reused outside
 *         the GUI.
 */
#pragma once

#include <array>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

class beat_log {
public:
    // Channels in column order. Keep in sync with the CSV header.
    enum ChannelIdx { ECG1 = 0, ECG2, ECG3, PPG, ABP, NUM_CHANNELS };

    // One channel's reading for a beat: time (x), amplitude (y), and the
    // blanking / threshold in effect when it was detected.
    struct sample {
        double x = 0.0;
        double y = 0.0;
        double blanking = 0.0;
        double threshold = 0.0;
    };

    // A beat = one sample per channel.
    struct beat {
        std::array<sample, NUM_CHANNELS> chan;
    };

    // Seed every beat's blanking/threshold columns (all channels) with the
    // config defaults. x/y stay zero until detection fills them.
    void setDefaultParams(double blanking, double threshold) {
        for (beat& b : m_beats)
            for (sample& s : b.chan) {
                s.blanking = blanking;
                s.threshold = threshold;
            }
    }

    // Record a peak for a channel in the pending buffer, keyed by global
    // time. A peak within ~10 ms of one already buffered (the same beat seen
    // through an overlapping window) is ignored.
    void logPeak(ChannelIdx ch, double x, double y) {
        auto& m = m_pending[ch];
        auto it = m.lower_bound(x - 0.01);
        if (it != m.end() && it->first <= x + 0.01) return;
        m[x] = y;
    }

    // Merge the pending buffer into the persistent table: a beat with the
    // same time overwrites; new ones insert in time order. Then empty the
    // buffer. Blanking/threshold columns are left as seeded.
    void flushPending() {
        for (int ch = 0; ch < NUM_CHANNELS; ++ch) {
            auto& pend = m_pending[ch];
            if (pend.empty()) continue;

            // Gather what's already in the table (the file content so far),
            // fold the buffer in (same time overwrites), and re-sort.
            std::map<double, double> merged;
            for (const beat& b : m_beats)
                if (b.chan[ch].x != 0.0) merged[b.chan[ch].x] = b.chan[ch].y;
            for (const auto& [t, y] : pend) merged[t] = y;

            for (beat& b : m_beats) { b.chan[ch].x = 0.0; b.chan[ch].y = 0.0; }
            std::size_t i = 0;
            for (const auto& [t, y] : merged) {
                if (i >= m_beats.size()) break;
                m_beats[i].chan[ch].x = t;
                m_beats[i].chan[ch].y = y;
                ++i;
            }
            pend.clear();   // empty the map for the next interval
        }
    }

    /**
     * @brief Write the log to a CSV file.
     * @param path  Output file path.
     * @return true on success, false if the file could not be opened.
     */
    bool writeCsv(const std::string& path) const;

private:
    // 57600 = 8 * 120 * 60, an 8 h recording at 120 bpm -- likely an overshoot.
    std::vector<beat> m_beats = std::vector<beat>(57600);

    // Rolling buffer of detected peaks per channel (global time -> amplitude),
    // flushed into m_beats and cleared every 30 s by the GUI.
    std::array<std::map<double, double>, NUM_CHANNELS> m_pending;
};