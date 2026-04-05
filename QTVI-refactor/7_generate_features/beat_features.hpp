#pragma once

#include "common.hpp"
#include "interx.hpp"
#include "find_foot_pulseox.hpp"
#include "dicrotic_notch.hpp"
#include <vector>
#include <cmath>
#include <algorithm>

/**
 * @file beat_features.hpp
 * @brief Extract intra-beat morphological features from a PPG waveform.
 *
 * Port of GetBeatFeaturesFromTemplate.m.
 * Computes timing (tP/tR at 20/50/80%), amplitude, slope, area,
 * R-peak-referenced intervals, and dicrotic-notch features per beat.
 */

namespace ppg {

    // ─── Feature struct ─────────────────────────────────────────────────────────

    /**
     * @brief Complete set of morphological features for a single PPG beat.
     */
    struct BeatFeatures {
        // Indices (global, into full PPG signal)
        int idx_begin = -1;   ///< Valley (onset) index.
        int idx_end = -1;   ///< Next valley index.
        int idx_foot = -1;   ///< Foot-detected onset.
        int idx_pos_slope = -1;   ///< Max positive slope location.
        int idx_systolic = -1;   ///< Systolic peak.
        int idx_neg_slope_b4 = -1; ///< Max negative slope before dicrotic notch.
        int idx_dnotch = -1;   ///< Dicrotic notch.
        int idx_diastolic = -1;   ///< Diastolic peak.
        int idx_neg_slope_after = -1; ///< Max negative slope after dicrotic notch.

        // tP crossing positions (relative to beat start, fractional samples)
        double tP_20_x = kNaN, tP_50_x = kNaN, tP_80_x = kNaN;
        double tP_20_x_inv = kNaN, tP_50_x_inv = kNaN, tP_80_x_inv = kNaN;

        // tR crossing positions
        double tR_20_x = kNaN, tR_50_x = kNaN, tR_80_x = kNaN;
        double tR_20_x_inv = kNaN, tR_50_x_inv = kNaN, tR_80_x_inv = kNaN;

        // tP amplitudes
        double tP_20_y = kNaN, tP_50_y = kNaN, tP_80_y = kNaN;
        double tP_20_y_inv = kNaN, tP_50_y_inv = kNaN, tP_80_y_inv = kNaN;

        // tR amplitudes
        double tR_20_y = kNaN, tR_50_y = kNaN, tR_80_y = kNaN;
        double tR_20_y_inv = kNaN, tR_50_y_inv = kNaN, tR_80_y_inv = kNaN;

        // Timing from tP-50 (msec)
        double msec_tP50_to_valley1 = kNaN, msec_tP50_to_foot = kNaN;
        double msec_tP50_to_tP20 = kNaN, msec_tP50_to_tP80 = kNaN;
        double msec_tP50_to_tP20_inv = kNaN, msec_tP50_to_tP80_inv = kNaN;
        double msec_tP50_to_pos_slope = kNaN, msec_tP50_to_systolic = kNaN;
        double msec_tP50_to_neg_pre = kNaN, msec_tP50_to_dnotch = kNaN;
        double msec_tP50_to_diastolic = kNaN, msec_tP50_to_neg_post = kNaN;
        double msec_tP50_to_tR20 = kNaN, msec_tP50_to_tR50 = kNaN, msec_tP50_to_tR80 = kNaN;
        double msec_tP50_to_tR20_inv = kNaN, msec_tP50_to_tR50_inv = kNaN, msec_tP50_to_tR80_inv = kNaN;

        // Durations (msec)
        double msec_beat_length = kNaN;
        double msec_dur_tP20 = kNaN, msec_dur_tP50 = kNaN, msec_dur_tP80 = kNaN;
        double msec_dur_tR20 = kNaN, msec_dur_tR50 = kNaN, msec_dur_tR80 = kNaN;

        // Raw amplitudes
        double amp_raw_valley = kNaN, amp_raw_foot = kNaN;
        double amp_raw_tP20 = kNaN, amp_raw_tP50 = kNaN, amp_raw_tP80 = kNaN;
        double amp_raw_tP20_inv = kNaN, amp_raw_tP50_inv = kNaN, amp_raw_tP80_inv = kNaN;
        double amp_raw_pos_slope = kNaN, amp_raw_systolic = kNaN;
        double amp_raw_neg_pre = kNaN, amp_raw_dnotch = kNaN;
        double amp_raw_diastolic = kNaN, amp_raw_neg_post = kNaN;
        double amp_raw_tR20 = kNaN, amp_raw_tR50 = kNaN, amp_raw_tR80 = kNaN;
        double amp_raw_tR20_inv = kNaN, amp_raw_tR50_inv = kNaN, amp_raw_tR80_inv = kNaN;

        // Baselined amplitudes (relative to valley)
        double amp_bl_foot = kNaN;
        double amp_bl_tP20 = kNaN, amp_bl_tP50 = kNaN, amp_bl_tP80 = kNaN;
        double amp_bl_tP20_inv = kNaN, amp_bl_tP50_inv = kNaN, amp_bl_tP80_inv = kNaN;
        double amp_bl_pos_slope = kNaN, amp_bl_systolic = kNaN;
        double amp_bl_neg_pre = kNaN, amp_bl_dnotch = kNaN;
        double amp_bl_diastolic = kNaN, amp_bl_neg_post = kNaN;
        double amp_bl_tR20 = kNaN, amp_bl_tR50 = kNaN, amp_bl_tR80 = kNaN;
        double amp_bl_tR20_inv = kNaN, amp_bl_tR50_inv = kNaN, amp_bl_tR80_inv = kNaN;

        // Area & proportional measures
        double area = kNaN, area_baselined = kNaN;
        double amp_delta_systolic = kNaN;
        double proportional_pulse_amp = kNaN;
        double abs_amp_foot = kNaN, abs_amp_peak = kNaN;

        // R-peak referenced timing (msec)
        double msec_R_to_valley1 = kNaN, msec_R_to_foot = kNaN;
        double msec_R_to_tP20 = kNaN, msec_R_to_tP50 = kNaN, msec_R_to_tP80 = kNaN;
        double msec_R_to_tP20_inv = kNaN, msec_R_to_tP50_inv = kNaN, msec_R_to_tP80_inv = kNaN;
        double msec_R_to_pos_slope = kNaN, msec_R_to_systolic = kNaN;
        double msec_R_to_neg_pre = kNaN, msec_R_to_dnotch = kNaN;
        double msec_R_to_diastolic = kNaN, msec_R_to_neg_post = kNaN;
        double msec_R_to_valley2 = kNaN;
        double msec_R_to_tR20 = kNaN, msec_R_to_tR50 = kNaN, msec_R_to_tR80 = kNaN;
        double msec_R_to_tR20_inv = kNaN, msec_R_to_tR50_inv = kNaN, msec_R_to_tR80_inv = kNaN;

        // Sleep stage
        double sleep_stage = kNaN;

        // SQI (filled externally)
        std::vector<double> sqi;
    };

    // ─── Internal helpers ───────────────────────────────────────────────────────

    namespace detail {

        /**
         * @brief Find where a horizontal line crosses a polyline, return the last crossing.
         * @param x  X coordinates (e.g. sample indices 1..N).
         * @param y  Y values of polyline.
         * @param level  Horizontal level to cross.
         * @return {x_crossing, y_level} or {NaN, NaN} if none.
         */
        inline Point2D find_last_crossing(
            const std::vector<double>& x,
            const std::vector<double>& y,
            double level)
        {
            std::vector<double> lx = { x.front(), x.back() };
            std::vector<double> ly = { level, level };
            auto pts = interx(x, y, lx, ly);
            if (pts.empty()) return { kNaN, kNaN };
            return pts.back();
        }

    } // namespace detail

    /**
     * @brief Extract features from a single PPG beat.
     * @param ppg            Full PPG signal.
     * @param valley_start   0-based index of the beat's first valley.
     * @param valley_end     0-based index of the beat's second valley.
     * @param sleep_states   Sleep stage annotation (same length as ppg).
     * @param ppg_r_idx      ECG R-peak sample index paired with this beat (-1 if unavailable).
     * @param ppg_fs         PPG sampling rate (Hz).
     * @param ecg_fs         ECG sampling rate (Hz).
     * @param dnotch_sp_ratio Optional dicrotic-notch ratio hint (NaN for auto).
     * @return Populated BeatFeatures struct.
     */
    inline BeatFeatures extract_beat_features(
        const std::vector<double>& ppg,
        int valley_start,
        int valley_end,
        const std::vector<double>& sleep_states,
        int ppg_r_idx,
        double ppg_fs,
        double ecg_fs,
        double dnotch_sp_ratio = kNaN)
    {
        BeatFeatures f;
        f.idx_begin = valley_start;
        f.idx_end = valley_end;

        int vs = valley_start;
        int ve = std::min(valley_end, static_cast<int>(ppg.size()) - 1);

        // Extract beat
        std::vector<double> beat(ppg.begin() + vs, ppg.begin() + ve + 1);
        const int bn = static_cast<int>(beat.size());
        if (bn < 4) return f;

        // Sleep stage (mode in range)
        {
            std::map<double, int> counts;
            for (int k = vs; k <= ve && k < static_cast<int>(sleep_states.size()); ++k)
                if (!std::isnan(sleep_states[k])) counts[sleep_states[k]]++;
            double mode_val = kNaN; int mode_cnt = 0;
            for (auto& [v, c] : counts) if (c > mode_cnt) { mode_val = v; mode_cnt = c; }
            f.sleep_stage = mode_val;
        }

        // Beat time axis (msec)
        std::vector<double> bt_msec(bn);
        for (int i = 0; i < bn; ++i) bt_msec[i] = i / ppg_fs * 1000.0;

        // Foot detection
        auto foot_res = find_foot_pulseox(beat);
        int beat_foot = foot_res.index;

        // Systolic peak
        auto it_max = std::max_element(beat.begin() + beat_foot, beat.end());
        int systolic_local = static_cast<int>(it_max - beat.begin());
        double max_amp = *it_max;

        // ── tP crossings (ascending limb: beat_foot → systolic) ──
        auto find_crossings = [&](
            const std::vector<double>& seg_x,
            const std::vector<double>& seg_y,
            double level, bool last)
            -> Point2D
            {
                auto pts = interx(seg_x, seg_y,
                    { seg_x.front(), seg_x.back() },
                    { level, level });
                if (pts.empty()) return { kNaN, kNaN };
                return last ? pts.back() : pts.front();
            };

        // tP levels (zeroed from beat start)
        std::vector<double> beat_zeroed(bn);
        for (int i = 0; i < bn; ++i) beat_zeroed[i] = beat[i] - beat[0];
        double max_z = beat_zeroed[systolic_local];
        double dif_tP = max_z - max_amp;
        double tP20_level = max_z * 0.2 - dif_tP;
        double tP50_level = max_z * 0.5 - dif_tP;
        double tP80_level = max_z * 0.8 - dif_tP;

        // Ascending limb crossings (beat start → systolic peak)
        if (systolic_local > 0) {
            std::vector<double> ax, ay;
            for (int i = 0; i <= systolic_local; ++i) {
                ax.push_back(static_cast<double>(i + 1));
                ay.push_back(beat[i]);
            }

            struct CrossTarget { double level; double* px; double* py; };
            CrossTarget asc_targets[] = {
                {tP20_level, &f.tP_20_x, &f.tP_20_y},
                {tP50_level, &f.tP_50_x, &f.tP_50_y},
                {tP80_level, &f.tP_80_x, &f.tP_80_y},
            };
            for (auto& t : asc_targets) {
                auto pt = find_crossings(ax, ay, t.level, true);
                *t.px = pt.x; *t.py = pt.y;
            }
        }

        // Descending limb tP-inverse crossings (systolic peak → end)
        if (systolic_local < bn - 1) {
            int seg_len = bn - systolic_local;
            std::vector<double> dx, dy;
            for (int i = 0; i < seg_len; ++i) {
                dx.push_back(static_cast<double>(i + 1));
                dy.push_back(beat[systolic_local + i]);
            }

            struct CrossTarget { double level; double* px; double* py; };
            CrossTarget desc_targets[] = {
                {tP20_level, &f.tP_20_x_inv, &f.tP_20_y_inv},
                {tP50_level, &f.tP_50_x_inv, &f.tP_50_y_inv},
                {tP80_level, &f.tP_80_x_inv, &f.tP_80_y_inv},
            };
            for (auto& t : desc_targets) {
                auto pt = find_crossings(dx, dy, t.level, true);
                *t.px = std::isnan(pt.x) ? kNaN : systolic_local + pt.x - 1;
                *t.py = pt.y;
            }
        }

        // Positive slope
        {
            auto d = diff(std::vector<double>(beat.begin() + beat_foot,
                beat.begin() + systolic_local + 1));
            if (!d.empty()) {
                int psi = static_cast<int>(std::max_element(d.begin(), d.end()) - d.begin());
                f.idx_pos_slope = vs + beat_foot + psi;
            }
        }
        int pos_slope_local = f.idx_pos_slope - vs;

        // Dicrotic notch (smoothing window scaled for sample rate: ~59ms)
        int dn_smooth = static_cast<int>(std::round(0.059 * ppg_fs));
        if (dn_smooth % 2 == 0) ++dn_smooth; // ensure odd
        int dn_local = dicrotic_notch(beat, dnotch_sp_ratio, dn_smooth);
        if (dn_local >= 0 && dn_local > systolic_local && dn_local < bn) {
            f.idx_dnotch = vs + dn_local;

            // Neg slope before dnotch
            auto seg = diff(std::vector<double>(beat.begin() + systolic_local,
                beat.begin() + dn_local + 1));
            if (!seg.empty()) {
                // max(-diff) = min(diff)
                int mi = static_cast<int>(std::min_element(seg.begin(), seg.end()) - seg.begin());
                f.idx_neg_slope_b4 = vs + systolic_local + mi;
            }

            // Neg slope after dnotch
            auto seg2 = diff(std::vector<double>(beat.begin() + dn_local, beat.end()));
            if (!seg2.empty()) {
                int mi = static_cast<int>(std::min_element(seg2.begin(), seg2.end()) - seg2.begin());
                int ns_after = dn_local + mi;
                f.idx_neg_slope_after = vs + ns_after;

                // Diastolic peak = max between dnotch and neg_slope_after
                auto it_dp = std::max_element(beat.begin() + dn_local,
                    beat.begin() + ns_after + 1);
                int dp_local = static_cast<int>(it_dp - beat.begin());
                if (dp_local != dn_local && dp_local != ns_after)
                    f.idx_diastolic = vs + dp_local;
            }
        }

        // ── tR crossings (zeroed from beat end) ──
        std::vector<double> beat_z_end(bn);
        for (int i = 0; i < bn; ++i) beat_z_end[i] = beat[i] - beat.back();
        double max_z_end = beat_z_end[systolic_local];
        double dif_tR = max_z_end - max_amp;
        double tR20_level = max_z_end * 0.2 - dif_tR;
        double tR50_level = max_z_end * 0.5 - dif_tR;
        double tR80_level = max_z_end * 0.8 - dif_tR;

        // Descending tR crossings (systolic peak → end)
        if (systolic_local < bn - 1) {
            int seg_len = bn - systolic_local;
            std::vector<double> dx, dy;
            for (int i = 0; i < seg_len; ++i) {
                dx.push_back(static_cast<double>(i + 1));
                dy.push_back(beat[systolic_local + i]);
            }

            struct CrossTarget { double level; double* px; double* py; };
            CrossTarget tR_desc[] = {
                {tR80_level, &f.tR_20_x, &f.tR_20_y},
                {tR50_level, &f.tR_50_x, &f.tR_50_y},
                {tR20_level, &f.tR_80_x, &f.tR_80_y},
            };
            for (auto& t : tR_desc) {
                auto pt = find_crossings(dx, dy, t.level, true);
                *t.px = std::isnan(pt.x) ? kNaN : systolic_local + pt.x - 1;
                *t.py = pt.y;
            }
        }

        // Ascending tR-inverse crossings (start → systolic peak)
        if (systolic_local > 0) {
            std::vector<double> ax, ay;
            for (int i = 0; i <= systolic_local; ++i) {
                ax.push_back(static_cast<double>(i + 1));
                ay.push_back(beat[i]);
            }

            struct CrossTarget { double level; double* px; double* py; };
            CrossTarget tR_asc[] = {
                {tR80_level, &f.tR_20_x_inv, &f.tR_20_y_inv},
                {tR50_level, &f.tR_50_x_inv, &f.tR_50_y_inv},
                {tR20_level, &f.tR_80_x_inv, &f.tR_80_y_inv},
            };
            for (auto& t : tR_asc) {
                auto pt = find_crossings(ax, ay, t.level, true);
                *t.px = pt.x; *t.py = pt.y;
            }
        }

        // ── Global indices ──
        f.idx_foot = vs + beat_foot;
        f.idx_systolic = vs + systolic_local;

        // Adjust tP/tR x to be relative to beat start (subtract 1 for 0-based)
        auto adj = [](double v) { return std::isnan(v) ? kNaN : v - 1.0; };
        f.tP_20_x = adj(f.tP_20_x); f.tP_50_x = adj(f.tP_50_x); f.tP_80_x = adj(f.tP_80_x);
        f.tP_20_x_inv = adj(f.tP_20_x_inv); f.tP_50_x_inv = adj(f.tP_50_x_inv); f.tP_80_x_inv = adj(f.tP_80_x_inv);
        f.tR_20_x = adj(f.tR_20_x); f.tR_50_x = adj(f.tR_50_x); f.tR_80_x = adj(f.tR_80_x);
        f.tR_20_x_inv = adj(f.tR_20_x_inv); f.tR_50_x_inv = adj(f.tR_50_x_inv); f.tR_80_x_inv = adj(f.tR_80_x_inv);

        // ── Timings from tP-50 (msec) ──
        double t50 = f.tP_50_x / ppg_fs * 1000.0;
        auto ms = [&](double local_x) { return local_x / ppg_fs * 1000.0 - t50; };
        f.msec_tP50_to_valley1 = bt_msec[0] - t50;
        f.msec_tP50_to_foot = static_cast<double>(beat_foot) / ppg_fs * 1000.0 - t50;
        f.msec_tP50_to_tP20 = ms(f.tP_20_x);
        f.msec_tP50_to_tP80 = ms(f.tP_80_x);
        f.msec_tP50_to_tP20_inv = ms(f.tP_20_x_inv);
        f.msec_tP50_to_tP80_inv = ms(f.tP_80_x_inv);
        f.msec_tP50_to_pos_slope = static_cast<double>(pos_slope_local) / ppg_fs * 1000.0 - t50;
        f.msec_tP50_to_systolic = static_cast<double>(systolic_local) / ppg_fs * 1000.0 - t50;

        int dn_l = (f.idx_dnotch >= 0) ? f.idx_dnotch - vs : -1;
        int ns_b4_l = (f.idx_neg_slope_b4 >= 0) ? f.idx_neg_slope_b4 - vs : -1;
        int ns_af_l = (f.idx_neg_slope_after >= 0) ? f.idx_neg_slope_after - vs : -1;
        int dia_l = (f.idx_diastolic >= 0) ? f.idx_diastolic - vs : -1;

        f.msec_tP50_to_neg_pre = (ns_b4_l >= 0) ? static_cast<double>(ns_b4_l) / ppg_fs * 1000.0 - t50 : kNaN;
        f.msec_tP50_to_dnotch = (dn_l >= 0) ? static_cast<double>(dn_l) / ppg_fs * 1000.0 - t50 : kNaN;
        f.msec_tP50_to_diastolic = (dia_l >= 0) ? static_cast<double>(dia_l) / ppg_fs * 1000.0 - t50 : kNaN;
        f.msec_tP50_to_neg_post = (ns_af_l >= 0) ? static_cast<double>(ns_af_l) / ppg_fs * 1000.0 - t50 : kNaN;
        f.msec_tP50_to_tR20 = ms(f.tR_20_x);
        f.msec_tP50_to_tR50 = ms(f.tR_50_x);
        f.msec_tP50_to_tR80 = ms(f.tR_80_x);
        f.msec_tP50_to_tR20_inv = ms(f.tR_20_x_inv);
        f.msec_tP50_to_tR50_inv = ms(f.tR_50_x_inv);
        f.msec_tP50_to_tR80_inv = ms(f.tR_80_x_inv);

        // Beat length
        f.msec_beat_length = bt_msec.back();

        // Durations
        auto dur = [&](double a, double b) {
            return std::abs(a / ppg_fs - b / ppg_fs) * 1000.0;
            };
        f.msec_dur_tP20 = dur(f.tP_20_x, f.tP_20_x_inv);
        f.msec_dur_tP50 = dur(f.tP_50_x, f.tP_50_x_inv);
        f.msec_dur_tP80 = dur(f.tP_80_x, f.tP_80_x_inv);
        f.msec_dur_tR20 = dur(f.tR_20_x, f.tR_20_x_inv);
        f.msec_dur_tR50 = dur(f.tR_50_x, f.tR_50_x_inv);
        f.msec_dur_tR80 = dur(f.tR_80_x, f.tR_80_x_inv);

        // ── Amplitudes ──
        double valley_amp = beat[0];
        f.amp_raw_valley = valley_amp;
        f.amp_raw_foot = safe_get(beat, beat_foot);
        f.amp_raw_systolic = safe_get(beat, systolic_local);
        f.amp_raw_pos_slope = safe_get(beat, pos_slope_local);
        f.amp_raw_neg_pre = (ns_b4_l >= 0) ? safe_get(beat, ns_b4_l) : kNaN;
        f.amp_raw_dnotch = (dn_l >= 0) ? safe_get(beat, dn_l) : kNaN;
        f.amp_raw_diastolic = (dia_l >= 0) ? safe_get(beat, dia_l) : kNaN;
        f.amp_raw_neg_post = (ns_af_l >= 0) ? safe_get(beat, ns_af_l) : kNaN;
        f.amp_raw_tP20 = f.tP_20_y; f.amp_raw_tP50 = f.tP_50_y; f.amp_raw_tP80 = f.tP_80_y;
        f.amp_raw_tP20_inv = f.tP_20_y_inv; f.amp_raw_tP50_inv = f.tP_50_y_inv; f.amp_raw_tP80_inv = f.tP_80_y_inv;
        f.amp_raw_tR20 = f.tR_20_y; f.amp_raw_tR50 = f.tR_50_y; f.amp_raw_tR80 = f.tR_80_y;
        f.amp_raw_tR20_inv = f.tR_20_y_inv; f.amp_raw_tR50_inv = f.tR_50_y_inv; f.amp_raw_tR80_inv = f.tR_80_y_inv;

        // Baselined
        auto bl = [&](double raw) { return raw - valley_amp; };
        f.amp_bl_foot = bl(f.amp_raw_foot);
        f.amp_bl_systolic = bl(f.amp_raw_systolic);
        f.amp_bl_pos_slope = bl(f.amp_raw_pos_slope);
        f.amp_bl_neg_pre = bl(f.amp_raw_neg_pre);
        f.amp_bl_dnotch = bl(f.amp_raw_dnotch);
        f.amp_bl_diastolic = bl(f.amp_raw_diastolic);
        f.amp_bl_neg_post = bl(f.amp_raw_neg_post);
        f.amp_bl_tP20 = bl(f.tP_20_y); f.amp_bl_tP50 = bl(f.tP_50_y); f.amp_bl_tP80 = bl(f.tP_80_y);
        f.amp_bl_tP20_inv = bl(f.tP_20_y_inv); f.amp_bl_tP50_inv = bl(f.tP_50_y_inv); f.amp_bl_tP80_inv = bl(f.tP_80_y_inv);
        f.amp_bl_tR20 = bl(f.tR_20_y); f.amp_bl_tR50 = bl(f.tR_50_y); f.amp_bl_tR80 = bl(f.tR_80_y);
        f.amp_bl_tR20_inv = bl(f.tR_20_y_inv); f.amp_bl_tR50_inv = bl(f.tR_50_y_inv); f.amp_bl_tR80_inv = bl(f.tR_80_y_inv);

        // Proportional pulse amplitude
        double dia_amp = (f.idx_diastolic >= 0) ? safe_get(beat, dia_l) : kNaN;
        f.proportional_pulse_amp = (!std::isnan(dia_amp) && f.amp_raw_systolic != 0.0)
            ? (f.amp_raw_systolic - dia_amp) / f.amp_raw_systolic
            : kNaN;

        // Area
        {
            std::vector<double> bz(bn);
            double bmin = *std::min_element(beat.begin(), beat.end());
            for (int i = 0; i < bn; ++i) bz[i] = beat[i] - bmin;
            std::vector<double> xv(bn);
            std::iota(xv.begin(), xv.end(), 0.0);
            f.area = trapz(xv, bz);
        }

        // Area baselined & delta systolic
        {
            double slope_val = (beat.back() - beat.front()) / (bn - 1);
            std::vector<double> ys(bn);
            for (int i = 0; i < bn; ++i) ys[i] = beat.front() + slope_val * i;
            std::vector<double> abs_diff(bn);
            for (int i = 0; i < bn; ++i) abs_diff[i] = std::abs(beat[i] - ys[i]);
            std::vector<double> xv(bn);
            std::iota(xv.begin(), xv.end(), 1.0);
            f.area_baselined = trapz(xv, abs_diff);
            f.amp_delta_systolic = std::abs(beat[systolic_local] - ys[systolic_local]);
        }

        // Normalised foot/peak amplitudes
        double bmin_val = *std::min_element(beat.begin(), beat.end());
        if (std::abs(bmin_val) > 1e-15) {
            std::vector<double> beat_norm(bn);
            for (int i = 0; i < bn; ++i) beat_norm[i] = beat[i] / bmin_val;
            f.abs_amp_foot = safe_get(beat_norm, beat_foot);
            f.abs_amp_peak = safe_get(beat_norm, systolic_local);
        }

        // ── R-peak referenced timings ──
        if (ppg_r_idx >= 0) {
            double r_offset = (vs / ppg_fs) * 1000.0 - (ppg_r_idx / ecg_fs) * 1000.0;
            auto rms = [&](double local_x) { return (local_x / ppg_fs) * 1000.0 + r_offset; };
            f.msec_R_to_valley1 = bt_msec[0] * 1000.0 + r_offset;
            f.msec_R_to_foot = rms(static_cast<double>(beat_foot));
            f.msec_R_to_pos_slope = rms(static_cast<double>(pos_slope_local));
            f.msec_R_to_tP20 = rms(f.tP_20_x);
            f.msec_R_to_tP50 = rms(f.tP_50_x);
            f.msec_R_to_tP80 = rms(f.tP_80_x);
            f.msec_R_to_tP20_inv = rms(f.tP_20_x_inv);
            f.msec_R_to_tP50_inv = rms(f.tP_50_x_inv);
            f.msec_R_to_tP80_inv = rms(f.tP_80_x_inv);
            f.msec_R_to_systolic = rms(static_cast<double>(systolic_local));
            f.msec_R_to_neg_pre = (ns_b4_l >= 0) ? rms(static_cast<double>(ns_b4_l)) : kNaN;
            f.msec_R_to_dnotch = (dn_l >= 0) ? rms(static_cast<double>(dn_l)) : kNaN;
            f.msec_R_to_diastolic = (dia_l >= 0) ? rms(static_cast<double>(dia_l)) : kNaN;
            f.msec_R_to_neg_post = (ns_af_l >= 0) ? rms(static_cast<double>(ns_af_l)) : kNaN;
            f.msec_R_to_valley2 = rms(static_cast<double>(bn - 1));
            f.msec_R_to_tR20 = rms(f.tR_20_x);
            f.msec_R_to_tR50 = rms(f.tR_50_x);
            f.msec_R_to_tR80 = rms(f.tR_80_x);
            f.msec_R_to_tR20_inv = rms(f.tR_20_x_inv);
            f.msec_R_to_tR50_inv = rms(f.tR_50_x_inv);
            f.msec_R_to_tR80_inv = rms(f.tR_80_x_inv);
        }

        return f;
    }

} // namespace ppg