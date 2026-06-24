/**
 * @file   beat_log.hpp
 * @brief  Logs the x and y locations, as well as the peak finding constants like blanking and threshold, and certain
 *         attributes pertaining to arrhythmia, of each beat as it is observed. If a beat is skipped over and not ever
 *         shown in the window, it will not be included in the log. One row per beat.
 *
 *         The output log is stored in the output folder in a subfolder log_initials. The header is:
 *
 *         "beat,ecg1_x,ecg1_y,ecg2_x,ecg2_y,ecg3_x,ecg3_y,ppg_x,ppg_y,abp_x,abp_y,"
 *         "blanking_ecg1,threshold_ecg1,blanking_ecg2,threshold_ecg2,"
 *         "blanking_ecg1,threshold_ecg1,blanking_ecg2,threshold_ecg2,"
 *         "blanking_ecg3,threshold_ecg3,blanking_ppg,threshold_ppg,blanking_abp,threshold_abp,"
 *         "marked_ecg1..marked_abp (annotation type containing the beat),"
 *         "post_ecg1..post_abp   (post_pvc/post_af/... if the beat follows an eligible arrhythmia)
 *
 *         The data is collected in a map which is keyed by the global x location of the beat, and is added to every time the user scrolls.
 *         Peaks are not duplicated because they are keyed by global x location.
 *
 *         The map is flushed into the log every 30 seconds. This empties the map.
 *
 *
 */
#pragma once

#include <array>
#include <cstddef>
#include <fstream>
#include <iomanip>
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
        int markType = 0;   // annotation type containing this beat, or "0"
        int postType = 0;   // "post_pvc"/"post_af"/... if this beat
        // immediately follows an eligible arrhythmia
    };

    // A beat = one sample per channel.
    struct beat {
        std::array<sample, NUM_CHANNELS> chan;
    };

    // Seed every beat's blanking/threshold columns (all channels) with the
    // config defaults. x/y stay zero until detection fills them.
    void setDefaultParams(double blanking, double threshold) {
        for (beat& b : all_beats_in_log)
            for (sample& s : b.chan) {
                s.blanking = blanking;
                s.threshold = threshold;
            }
    }

    struct pendingPeak { double y, blanking, threshold; int markType, postType; };
    void logPeak(ChannelIdx ch, double x, double y, double blanking, double threshold, int markType = 0, int postType = 0) {
        //Record a peak for a channel in the pending buffer, keyed by global time 
        auto& m = m_pending[ch];
        m[x] = { y, blanking, threshold, markType, postType };
    }


    // Merge the pending buffer into the persistent table: a beat with the
    // same time overwrites; new ones insert in time order. Then empty the
    // buffer. Blanking/threshold columns are left as seeded.
    void flushPending() {
        for (int ch = 0; ch < NUM_CHANNELS; ++ch) {
            auto& pend = m_pending[ch];
            if (pend.empty()) continue;

            std::map<double, pendingPeak> merged;
            for (const beat& b : all_beats_in_log)
                if (b.chan[ch].x != 0.0)
                    merged[b.chan[ch].x] = { b.chan[ch].y, b.chan[ch].blanking, b.chan[ch].threshold, b.chan[ch].markType, b.chan[ch].postType };
            for (const auto& [t, pv] : pend) merged[t] = pv;

            for (beat& b : all_beats_in_log) { b.chan[ch].x = 0.0; b.chan[ch].y = 0.0; b.chan[ch].markType = 0; b.chan[ch].postType = 0; }
            std::size_t i = 0;
            for (const auto& [t, pv] : merged) {
                if (i >= all_beats_in_log.size()) break;
                all_beats_in_log[i].chan[ch].x = t;
                all_beats_in_log[i].chan[ch].y = pv.y;
                all_beats_in_log[i].chan[ch].blanking = pv.blanking;
                all_beats_in_log[i].chan[ch].threshold = pv.threshold;
                all_beats_in_log[i].chan[ch].markType = pv.markType;
                all_beats_in_log[i].chan[ch].postType = pv.postType;
                ++i;
            }
            pend.clear();
        }
    }

    // Remove any pending or committed peaks for a channel whose time falls in
    // [t0, t1]. Used when a threshold/blanking override changes detection in a
    // region: dropped beats must not linger in the log -- survivors are re-logged
    // by the next redraw. Committed slots are cleared by zeroing x (the
    // empty-slot marker); blanking/threshold columns are left as seeded.
    void removeInRange(ChannelIdx ch, double t0, double t1) {
        if (t0 > t1) { double tmp = t0; t0 = t1; t1 = tmp; }
        auto& pend = m_pending[ch];
        for (auto it = pend.begin(); it != pend.end(); ) {
            if (it->first >= t0 && it->first <= t1) it = pend.erase(it);
            else ++it;
        }
        for (beat& b : all_beats_in_log) {
            sample& s = b.chan[ch];
            if (s.x != 0.0 && s.x >= t0 && s.x <= t1) { s.x = 0.0; s.y = 0.0; s.markType = 0; s.postType = 0; }
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
    std::vector<beat> all_beats_in_log = std::vector<beat>(57600);

    // Rolling buffer of detected peaks per channel (global time -> amplitude),
    // flushed into all_beats_in_log and cleared every 30 s by the GUI.
    std::array<std::map<double, pendingPeak>, NUM_CHANNELS> m_pending;
    std::map<double, std::array<double, 3>> m_accel;
    void logAccel(double x, double ax, double ay, double az) { m_accel[x] = { ax, ay, az }; }

};

inline bool beat_log::writeCsv(const std::string& path) const {
    /*
        This empties the map and fills the CSV file
    */
    std::ofstream f(path);
    if (!f.is_open()) return false;

    f << "beat,ecg1_x,ecg1_y,ecg2_x,ecg2_y,ecg3_x,ecg3_y,ppg_x,ppg_y,abp_x,abp_y,"
        "blanking_ecg1,threshold_ecg1,blanking_ecg2,threshold_ecg2,"
        "blanking_ecg3,threshold_ecg3,blanking_ppg,threshold_ppg,blanking_abp,threshold_abp,"
        "marked_ecg1,marked_ecg2,marked_ecg3,marked_ppg,marked_abp,marked_accel,accelx_val, accely_val, accelz_val,"
        "post_ecg1,post_ecg2,post_ecg3,post_ppg,post_abp\n";


    f << std::fixed << std::setprecision(6);

    for (size_t i = 0; i < all_beats_in_log.size(); ++i) {
        f << i;
        for (const sample& s : all_beats_in_log[i].chan)      // x,y
            f << ',' << s.x << ',' << s.y;
        for (const sample& s : all_beats_in_log[i].chan)      // blanking,threshold
            f << ',' << s.blanking << ',' << s.threshold;
        for (const sample& s : all_beats_in_log[i].chan)      // markType
            f << ',' << s.markType;
        for (const sample& s : all_beats_in_log[i].chan)      // postType  <-- this loop
            f << ',' << s.postType;
        for (const sample& s : all_beats_in_log[i].chan)      // marked_ecg1..marked_abp
            f << ',' << s.markType;

        // marked_accelx/y/z: accel at this row's ECG1 x, zeros if none
        const double ex = all_beats_in_log[i].chan[ECG1].x;
        std::array<double, 3> a{ 0.0, 0.0, 0.0 };
        auto it = (ex != 0.0) ? m_accel.find(ex) : m_accel.end();
        if (it != m_accel.end()) a = it->second;
        f << ',' << a[0] << ',' << a[1] << ',' << a[2];

        for (const sample& s : all_beats_in_log[i].chan)      // post_ecg1..post_abp
            f << ',' << s.postType;

        f << '\n';
    }
    return true;
}