/**
 * @file   beat_log.hpp
 * @brief  Logs the x and y locations, as well as the peak finding constants like blanking and threshold, and certain
 *         attributes pertaining to arrhythmia, of each beat as it is observed. If a beat is skipped over and not ever
 *         shown in the window, it will not be included in the log. One row per beat.
 *
 *         The output log is stored in the output folder in a subfolder log_initials. The columns are, in order:
 *
 *           beat
 *           <chan>_x, <chan>_y                       for ecg1, ecg2, ecg3, ppg, abp, art, art_pulm
 *           blanking_<chan>, threshold_<chan>        for the same channels
 *           marked_<chan>                            annotation type containing the beat (0 if none)
 *           marked_accel, accelx_val, accely_val, accelz_val
 *           post_<chan>                              post_pvc/post_af/... if the beat follows an eligible arrhythmia
 *           inverted_<chan>                          1 if the channel is inverted at this beat (checkbox XOR invert-region), else 0
 *
 *         The data is collected in a map which is keyed by the global x location of the beat, and is added to every time the user scrolls.
 *         Peaks are not duplicated because they are keyed by global x location.
 *
 *         The map is flushed into the log every 30 seconds. This empties the map.
 *
 *         The table starts small (kInitialRows) to save memory and grows a chunk
 *         at a time (kGrowChunk) only when a flush needs more room, so growth is
 *         rare rather than per-beat. writeCsv skips fully-empty rows so the write
 *         cost tracks the real beat count, not the allocated capacity.
 *
 *         NOTE: marked_accel is currently written as a 0 placeholder. Accel is drawn on the shared PPG_ACCEL chart and
 *               is not yet its own markable channel, so no accel annotation type is recorded. Once PPG and accel are
 *               split into separate charts, accel becomes a normal ChannelIdx and this column fills itself.
 */
#pragma once

#include <algorithm>
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
    enum ChannelIdx { ECG1 = 0, ECG2, ECG3, PPG, ABP, ART, ART_PULM, NUM_CHANNELS };

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
        int inverted = 0;   // 1 if this channel was inverted at this beat, else 0
    };

    // A beat = one sample per channel.
    struct beat {
        std::array<sample, NUM_CHANNELS> chan;
    };

    // Seed every beat's blanking/threshold columns (all channels) with the
    // config defaults, and remember them so rows added by later chunked growth
    // inherit the same defaults. x/y stay zero until detection fills them.
    void setDefaultParams(double blanking, double threshold) {
        m_defaultBlanking = blanking;
        m_defaultThreshold = threshold;
        for (beat& b : all_beats_in_log)
            for (sample& s : b.chan) {
                s.blanking = blanking;
                s.threshold = threshold;
            }
    }

    struct pendingPeak { double y, blanking, threshold; int markType, postType, inverted; };
    void logPeak(ChannelIdx ch, double x, double y, double blanking, double threshold,
        int markType = 0, int postType = 0, int inverted = 0) {
        //Record a peak for a channel in the pending buffer, keyed by global time 
        auto& m = m_pending[ch];
        m[x] = { y, blanking, threshold, markType, postType, inverted };
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
                    merged[b.chan[ch].x] = { b.chan[ch].y, b.chan[ch].blanking, b.chan[ch].threshold, b.chan[ch].markType, b.chan[ch].postType, b.chan[ch].inverted };
            for (const auto& [t, pv] : pend) merged[t] = pv;

            // Grow (a chunk at a time) so every merged beat has a row -- no
            // silent drop of beats past a fixed capacity.
            ensureCapacity(merged.size());

            for (beat& b : all_beats_in_log) { b.chan[ch].x = 0.0; b.chan[ch].y = 0.0; b.chan[ch].markType = 0; b.chan[ch].postType = 0; b.chan[ch].inverted = 0; }
            std::size_t i = 0;
            for (const auto& [t, pv] : merged) {
                all_beats_in_log[i].chan[ch].x = t;
                all_beats_in_log[i].chan[ch].y = pv.y;
                all_beats_in_log[i].chan[ch].blanking = pv.blanking;
                all_beats_in_log[i].chan[ch].threshold = pv.threshold;
                all_beats_in_log[i].chan[ch].markType = pv.markType;
                all_beats_in_log[i].chan[ch].postType = pv.postType;
                all_beats_in_log[i].chan[ch].inverted = pv.inverted;
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
            if (s.x != 0.0 && s.x >= t0 && s.x <= t1) { s.x = 0.0; s.y = 0.0; s.markType = 0; s.postType = 0; s.inverted = 0; }
        }
    }

    /**
     * @brief Write the log to a CSV file.
     * @param path  Output file path.
     * @return true on success, false if the file could not be opened.
     */
    bool writeCsv(const std::string& path) const;

private:
    // Start small to save memory; grow a chunk at a time on demand. 57600 rows
    // (8 h at 120 bpm) was the old fixed overshoot -- most recordings use far
    // fewer, so we allocate a fraction up front and only grow if a busy file
    // actually needs it.
    static constexpr std::size_t kInitialRows = 5700;
    static constexpr std::size_t kGrowChunk = 5700;   // grow a chunk at a time, not per beat

    std::vector<beat> all_beats_in_log = std::vector<beat>(kInitialRows);

    // Seeded config defaults, remembered so rows created by ensureCapacity()
    // inherit them instead of writing zeros for blanking/threshold.
    double m_defaultBlanking = 0.0;
    double m_defaultThreshold = 0.0;

    // Grow the table to hold at least `need` rows, a chunk at a time so a flush
    // rarely reallocates. New rows inherit the seeded default params.
    void ensureCapacity(std::size_t need) {
        if (need <= all_beats_in_log.size()) return;
        const std::size_t oldSize = all_beats_in_log.size();
        const std::size_t grown = std::max(need, oldSize + kGrowChunk);
        all_beats_in_log.resize(grown);
        for (std::size_t i = oldSize; i < grown; ++i)
            for (sample& s : all_beats_in_log[i].chan) {
                s.blanking = m_defaultBlanking;
                s.threshold = m_defaultThreshold;
            }
    }

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
        "art_x,art_y,art_pulm_x,art_pulm_y,"
        "blanking_ecg1,threshold_ecg1,blanking_ecg2,threshold_ecg2,"
        "blanking_ecg3,threshold_ecg3,blanking_ppg,threshold_ppg,blanking_abp,threshold_abp,"
        "blanking_art,threshold_art,blanking_art_pulm,threshold_art_pulm,"
        "marked_ecg1,marked_ecg2,marked_ecg3,marked_ppg,marked_abp,marked_art,marked_art_pulm,"
        "marked_accel,accelx_val,accely_val,accelz_val,"
        "post_ecg1,post_ecg2,post_ecg3,post_ppg,post_abp,post_art,post_art_pulm,"
        "inverted_ecg1,inverted_ecg2,inverted_ecg3,inverted_ppg,inverted_abp,inverted_art,inverted_art_pulm\n";


    f << std::fixed << std::setprecision(6);

    for (size_t i = 0; i < all_beats_in_log.size(); ++i) {
        // Skip fully-empty rows: a row with no beat on any channel is just a
        // line of zeros and, at the allocated capacity, dominates the write
        // cost. Writing only populated rows keeps the write proportional to the
        // real beat count.
        bool anyBeat = false;
        for (const sample& s : all_beats_in_log[i].chan)
            if (s.x != 0.0) { anyBeat = true; break; }
        if (!anyBeat) continue;

        f << i;
        for (const sample& s : all_beats_in_log[i].chan)      // <chan>_x, <chan>_y
            f << ',' << s.x << ',' << s.y;
        for (const sample& s : all_beats_in_log[i].chan)      // blanking_<chan>, threshold_<chan>
            f << ',' << s.blanking << ',' << s.threshold;
        for (const sample& s : all_beats_in_log[i].chan)      // marked_<chan> (annotation type)
            f << ',' << s.markType;

        // marked_accel + accelx/y/z: accel value at this row's ECG1 x, zeros if none.
        // marked_accel is a placeholder 0 until accel becomes its own markable channel.
        const double ex = all_beats_in_log[i].chan[ECG1].x;
        std::array<double, 3> a{ 0.0, 0.0, 0.0 };
        auto it = (ex != 0.0) ? m_accel.find(ex) : m_accel.end();
        if (it != m_accel.end()) a = it->second;
        f << ',' << 0 << ',' << a[0] << ',' << a[1] << ',' << a[2];

        for (const sample& s : all_beats_in_log[i].chan)      // post_<chan>
            f << ',' << s.postType;

        for (const sample& s : all_beats_in_log[i].chan)      // inverted_<chan>
            f << ',' << s.inverted;

        f << '\n';
    }
    return true;
}