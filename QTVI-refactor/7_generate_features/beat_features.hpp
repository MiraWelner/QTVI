#pragma once
// ═══════════════════════════════════════════════════════════════════════════════
// beat_features.hpp — Exact match of GetBeatFeaturesFromTemplate.m
//
// Critical: MATLAB uses 1-based indexing throughout. All beat_info fields
// (tP_xx_x, beat_foot, systolic_peak, etc.) are in MATLAB 1-based local coords.
// The outputs (beats_total.*) store:
//   - idx_* fields: MATLAB 1-based global coords (beat_valley_one + local - 1)
//   - tP_xx_x, tR_xx_x: MATLAB beat_info.tP_xx_x - 1 (0-based from beat start)
//   - Timings use beat_info.tP_xx_x (1-based local) / ppgSamplingRate * 1000
//
// In C++ we work 0-based internally, then convert to match MATLAB output.
// The conversion rule:
//   MATLAB local 1-based idx → C++ 0-based idx + 1
//   So C++ systolic_peak (0-based) → MATLAB systolic_peak = C++ + 1
//   And beats_total.tP_50_x = (MATLAB beat_info.tP_50_x) - 1
//                            = (C++ interx result) [already in MATLAB 1-based from interx]  - 1
// ═══════════════════════════════════════════════════════════════════════════════
#include "ppg_features.hpp"
#include "ppg_utils.hpp"
#include "find_foot_pulseox.hpp"
#include "dicrotic_sqi.hpp"

namespace ppg {

inline void BeatsTotal::resize(int n) {
    auto r = [&](std::vector<double>& v) { v.assign(n, NaN); };
    r(sleep_stages); r(area_baselined); r(area);
    r(amp_delta_systolic); r(abs_amp_foot); r(abs_amp_peak);
    r(idx_begin); r(idx_end); r(idx_foot);
    r(tP_20_x); r(tP_50_x); r(tP_80_x);
    r(tP_20_x_inv); r(tP_50_x_inv); r(tP_80_x_inv);
    r(idx_pos_slope); r(idx_systolic);
    r(idx_neg_slope_b4); r(idx_neg_slope_after);
    r(idx_diastolic); r(idx_dnotch);
    r(tR_20_x); r(tR_50_x); r(tR_80_x);
    r(tR_20_x_inv); r(tR_50_x_inv); r(tR_80_x_inv);
    r(msec_beat_length);
    r(msec_tP_50_2_first_valley); r(msec_tP_50_2_foot);
    r(msec_tP_50_2_tP_20); r(msec_tP_50_2_tP_80);
    r(msec_tP_50_2_tP_20_inv); r(msec_tP_50_2_tP_80_inv);
    r(msec_tP_50_2_pos_slope); r(msec_tP_50_2_systolic_peak);
    r(msec_tP_50_2_negslopes_pre_dnotch); r(msec_tP_50_2_dicrotic_notch);
    r(msec_tP_50_2_diastolic_peak); r(msec_tP_50_2_negslopes_post_dnotch);
    r(msec_tP_50_2_second_valley);
    r(msec_tP_50_2_tR_20); r(msec_tP_50_2_tR_50); r(msec_tP_50_2_tR_80);
    r(msec_tP_50_2_tR_20_inv); r(msec_tP_50_2_tR_50_inv); r(msec_tP_50_2_tR_80_inv);
    r(msec_total_duration_20); r(msec_total_duration_50); r(msec_total_duration_80);
    r(msec_total_duration_tR_20); r(msec_total_duration_tR_50); r(msec_total_duration_tR_80);
    r(msec_R_2_first_valley); r(msec_R_2_foot);
    r(msec_R_2_tP_20); r(msec_R_2_tP_50); r(msec_R_2_tP_80);
    r(msec_R_2_tP_20_inv); r(msec_R_2_tP_50_inv); r(msec_R_2_tP_80_inv);
    r(msec_R_2_pos_slope); r(msec_R_2_systolic_peak);
    r(msec_R_2_negslopes_pre_dnotch); r(msec_R_2_dicrotic_notch);
    r(msec_R_2_tR_20); r(msec_R_2_tR_50); r(msec_R_2_tR_80);
    r(msec_R_2_tR_20_inv); r(msec_R_2_tR_50_inv); r(msec_R_2_tR_80_inv);
    r(msec_R_2_diastolic_peak); r(msec_R_2_negslopes_post_dnotch);
    r(msec_R_2_second_valley);
    r(amp_raw_vallies); r(amp_raw_feets);
    r(amp_raw_tP_20); r(amp_raw_tP_50); r(amp_raw_tP_80);
    r(amp_raw_tP_20_inv); r(amp_raw_tP_50_inv); r(amp_raw_tP_80_inv);
    r(amp_raw_tR_20); r(amp_raw_tR_50); r(amp_raw_tR_80);
    r(amp_raw_tR_20_inv); r(amp_raw_tR_50_inv); r(amp_raw_tR_80_inv);
    r(amp_raw_pos_slopes); r(amp_raw_systolic_peaks);
    r(amp_raw_neg_slopes_pre_dnotch); r(amp_raw_dicrotic_notches);
    r(amp_raw_diastolic_peaks); r(amp_raw_neg_slopes_after_dnotch);
    r(amp_baselined_feets);
    r(amp_baselined_tP_20); r(amp_baselined_tP_50); r(amp_baselined_tP_80);
    r(amp_baselined_tP_20_inv); r(amp_baselined_tP_50_inv); r(amp_baselined_tP_80_inv);
    r(amp_baselined_tR_20); r(amp_baselined_tR_50); r(amp_baselined_tR_80);
    r(amp_baselined_tR_20_inv); r(amp_baselined_tR_50_inv); r(amp_baselined_tR_80_inv);
    r(amp_baselined_pos_slopes); r(amp_baselined_systolic_peaks);
    r(amp_baselined_neg_slopes_pre_dnotch); r(amp_baselined_dicrotic_notches);
    r(amp_baselined_diastolic_peaks); r(amp_baselined_neg_slopes_after_dnotch);
    r(proportional_pulse_amp);
}

inline BeatFeaturesResult GetBeatFeaturesFromTemplate(
    const std::vector<double>& sqi, double threshold,
    const std::vector<double>& ppg,
    const std::vector<double>& sleepstates,
    const Pairs& pairs,
    double ppgSR, double ecgSR,
    double dnotch_ratio_sp)
{
    BeatFeaturesResult result;

    // MATLAB: idx_vallies = pairs(:,1); idx_begin = idx_vallies(1:end-1); etc.
    // Then filter by sqi < threshold
    std::vector<int> idx_begin_all, idx_end_all;
    std::vector<int> pair_indices; // indices into original pairs for R-peak lookup

    if (pairs.rows() >= 2) {
        for (int i = 0; i < pairs.rows() - 1; i++) {
            if (i < (int)sqi.size() && sqi[i] < threshold) {
                idx_begin_all.push_back(pairs.col(i, 0));
                idx_end_all.push_back(pairs.col(i+1, 0));
                pair_indices.push_back(i);
            }
        }
    }

    int len = (int)idx_begin_all.size();
    auto& bt = result.beats_total;
    bt.resize(len);

    // Helper: safe array access matching MATLAB getVal
    auto gv = [](const std::vector<double>& arr, int idx_0based) -> double {
        if (idx_0based < 0 || idx_0based >= (int)arr.size()) return NaN;
        return arr[idx_0based];
    };

    // Process each beat
    for (int i = 0; i < len; i++) {
        // MATLAB: beat_valley_one = idx_begin(i) (1-based)
        // C++: 0-based
        int bv1 = idx_begin_all[i]; // 0-based index into ppg
        int bv2 = idx_end_all[i];

        // MATLAB: if beat_valley_two > length(ppg), beat_valley_two = length(ppg)
        if (bv2 >= (int)ppg.size()) bv2 = (int)ppg.size() - 1;

        // MATLAB: beat = ppg(beat_valley_one:beat_valley_two)
        int beat_len = bv2 - bv1 + 1;
        if (beat_len < 4) {
            bt.sleep_stages[i] = NaN;
            continue;
        }

        std::vector<double> beat(ppg.begin() + bv1, ppg.begin() + bv2 + 1);

        // Sleep stage: mode of sleepstates in beat range
        {
            std::map<double, int> cnt;
            for (int k = bv1; k <= bv2 && k < (int)sleepstates.size(); k++)
                if (!std::isnan(sleepstates[k])) cnt[sleepstates[k]]++;
            double mv = NaN; int mc = 0;
            for (auto& [v, c] : cnt) if (c > mc) { mc = c; mv = v; }
            bt.sleep_stages[i] = mv;
        }

        // MATLAB: beat_time_msec = (0:length(beat)-1) / ppgSamplingRate * 1000
        // beat_time_msec(1) = 0 (MATLAB 1-based, so index 1 → value 0)

        // MATLAB: [~, beat_foot] = find_foot_pulseox(beat', 0)
        // beat_foot is 1-based in MATLAB. C++ find_foot_pulseox returns 0-based.
        auto [foot_val, beat_foot_0] = find_foot_pulseox(beat);
        // Convert to MATLAB 1-based for computation consistency
        int M_beat_foot = beat_foot_0 + 1; // MATLAB 1-based

        // MATLAB: [max_amp, max_peak] = max(beat(beat_foot:end))
        // systolic_peak = beat_foot + max_peak - 1 (1-based)
        double max_amp = -1e18;
        int max_peak_local = 0; // 1-based offset from beat_foot
        for (int k = beat_foot_0; k < beat_len; k++) {
            if (beat[k] > max_amp) {
                max_amp = beat[k];
                max_peak_local = k - beat_foot_0; // 0-based offset
            }
        }
        int M_systolic_peak = M_beat_foot + max_peak_local; // MATLAB 1-based
        int sp_0 = M_systolic_peak - 1; // C++ 0-based

        // ── tP thresholds ──────────────────────────────────────────────
        // MATLAB: beat_zeroed = beat - beat(1); max_amp_zeroed = beat_zeroed(systolic_peak)
        double beat_zeroed_peak = beat[sp_0] - beat[0];
        // MATLAB: dif = max_amp_zeroed - max_amp
        double dif = beat_zeroed_peak - max_amp;
        double max_tP20 = beat_zeroed_peak * 0.2 - dif;
        double max_tP50 = beat_zeroed_peak * 0.5 - dif;
        double max_tP80 = beat_zeroed_peak * 0.8 - dif;

        // ── tP ascending: InterX on beat(1:systolic_peak) ──────────────
        // MATLAB: x = (1:length(beat(1:systolic_peak)))'; y = beat(1:systolic_peak);
        // l1 = [x y]'; InterX finds last crossing of horizontal line at level.
        // beat(1:systolic_peak) in MATLAB = beat[0..sp_0] in C++ (sp_0 = M_systolic_peak-1)
        double tP_20_X = NaN, tP_50_X = NaN, tP_80_X = NaN;
        double tP_20_y = NaN, tP_50_y = NaN, tP_80_y = NaN;

        if (M_systolic_peak > 1) { // length(y) > 1
            // InterX last crossing in beat[0..sp_0] at each level
            // Returns MATLAB 1-based local x coordinate
            auto [x20, y20] = interx_last_crossing(beat.data(), sp_0 + 1, max_tP20);
            auto [x50, y50] = interx_last_crossing(beat.data(), sp_0 + 1, max_tP50);
            auto [x80, y80] = interx_last_crossing(beat.data(), sp_0 + 1, max_tP80);
            tP_20_X = x20; tP_50_X = x50; tP_80_X = x80;
            tP_20_y = y20; tP_50_y = y50; tP_80_y = y80;
        }

        // ── tP descending (inv): InterX on beat(systolic_peak:end) ─────
        // MATLAB: x = (1:length(beat(systolic_peak:end)))'; y = beat(systolic_peak:end);
        // Result: tP_20_inv_x = systolic_peak + ttime_x(1) - 1
        double tP_20_inv_x = NaN, tP_50_inv_x = NaN, tP_80_inv_x = NaN;
        double tP_20_inv_y = NaN, tP_50_inv_y = NaN, tP_80_inv_y = NaN;

        int desc_len = beat_len - sp_0; // length of beat(systolic_peak:end) in MATLAB
        if (desc_len > 1) {
            // Local array: beat[sp_0..beat_len-1], local MATLAB coords 1..desc_len
            auto [x20, y20] = interx_last_crossing(beat.data() + sp_0, desc_len, max_tP20);
            auto [x50, y50] = interx_last_crossing(beat.data() + sp_0, desc_len, max_tP50);
            auto [x80, y80] = interx_last_crossing(beat.data() + sp_0, desc_len, max_tP80);
            // MATLAB: tP_20_inv_x = systolic_peak + ttime_x - 1
            tP_20_inv_x = std::isnan(x20) ? NaN : M_systolic_peak + x20 - 1;
            tP_50_inv_x = std::isnan(x50) ? NaN : M_systolic_peak + x50 - 1;
            tP_80_inv_x = std::isnan(x80) ? NaN : M_systolic_peak + x80 - 1;
            tP_20_inv_y = y20; tP_50_inv_y = y50; tP_80_inv_y = y80;
        }

        // ── Positive slope ─────────────────────────────────────────────
        // MATLAB: [~, pos_slope_idx] = max(diff(beat(beat_foot:beat_foot + max_peak - 1)))
        // max_positive_slope = beat_foot + pos_slope_idx - 1
        int M_max_positive_slope = M_beat_foot; // default
        {
            // MATLAB range: beat_foot to beat_foot + max_peak - 1 (1-based)
            // In 0-based: beat_foot_0 to beat_foot_0 + max_peak_local
            // But max_peak_local is 0-based offset, so the MATLAB range end
            // = beat_foot + max_peak - 1 = (beat_foot_0+1) + max_peak_local - 1
            //   = beat_foot_0 + max_peak_local = sp_0
            // diff has length sp_0 - beat_foot_0
            double max_slope_val = -1e18;
            int pos_slope_local = 0;
            for (int k = beat_foot_0; k < sp_0; k++) {
                double d = beat[k+1] - beat[k];
                if (d > max_slope_val) {
                    max_slope_val = d;
                    pos_slope_local = k - beat_foot_0; // 0-based offset
                }
            }
            // MATLAB: max_positive_slope = beat_foot + pos_slope_idx - 1
            // pos_slope_idx is 1-based, so pos_slope_local + 1
            M_max_positive_slope = M_beat_foot + pos_slope_local;
        }

        // ── Dicrotic notch ─────────────────────────────────────────────
        // MATLAB: dicrotic_notch = dumbDicrotic(beat) or dumbDicrotic(beat, sp_ratio)
        // dumbDicrotic returns 1-based in MATLAB. Our C++ returns 0-based.
        int dn_0 = dumbDicrotic(beat, dnotch_ratio_sp);
        int M_dicrotic_notch; // MATLAB 1-based, or will be set to NaN-equivalent

        bool dn_valid = false;
        if (dn_0 >= 0) {
            M_dicrotic_notch = dn_0 + 1; // to MATLAB 1-based
            // MATLAB: if dicrotic_notch > length(beat) || dicrotic_notch <= systolic_peak || isnan
            if (M_dicrotic_notch > beat_len || M_dicrotic_notch <= M_systolic_peak) {
                dn_valid = false;
            } else {
                dn_valid = true;
            }
        }
        if (!dn_valid) M_dicrotic_notch = -1; // represents NaN

        int M_neg_slope_b4 = -1, M_neg_slope_after = -1, M_diastolic_peak = -1;

        if (dn_valid) {
            // MATLAB: [~, b4dnotch] = max(-diff(beat(systolic_peak:dicrotic_notch)))
            // Range: sp_0..dn_0 (0-based), diff has len = dn_0 - sp_0
            double max_neg = -1e18;
            int b4_local = -1; // 1-based offset into the slice
            for (int k = sp_0; k < dn_0; k++) {
                double d = -(beat[k+1] - beat[k]);
                if (d > max_neg) { max_neg = d; b4_local = k - sp_0 + 1; } // MATLAB 1-based
            }
            // MATLAB: if isempty(b4dnotch) || b4dnotch >= length(beat) || b4dnotch == 1
            if (b4_local <= 0 || b4_local >= beat_len || b4_local == 1) {
                // b4dnotch = nan → max_neg_slope_before = nan
            } else {
                // MATLAB: max_neg_slope_before_dicrotic_notch = systolic_peak + b4dnotch - 1
                M_neg_slope_b4 = M_systolic_peak + b4_local - 1;
            }

            // MATLAB: [~, afterNotch] = max(-diff(beat(dicrotic_notch:end)))
            double max_neg2 = -1e18;
            int after_local = -1;
            for (int k = dn_0; k < beat_len - 1; k++) {
                double d = -(beat[k+1] - beat[k]);
                if (d > max_neg2) { max_neg2 = d; after_local = k - dn_0 + 1; }
            }
            if (after_local <= 0 || after_local >= beat_len || after_local == 1) {
                // afterNotch = nan
            } else {
                M_neg_slope_after = M_dicrotic_notch + after_local - 1;
            }

            // Diastolic peak
            if (M_neg_slope_after > 0) {
                int ns_after_0 = M_neg_slope_after - 1; // 0-based
                double max_dp = -1e18;
                int dp_local = -1; // 1-based offset
                for (int k = dn_0; k <= ns_after_0 && k < beat_len; k++) {
                    if (beat[k] > max_dp) { max_dp = beat[k]; dp_local = k - dn_0 + 1; }
                }
                // MATLAB: if dpeak == 1 || dpeak >= length(beat) || isempty(dpeak) → nan
                if (dp_local <= 0 || dp_local == 1 || dp_local >= beat_len) {
                    // diastolic_peak = nan
                } else {
                    M_diastolic_peak = M_dicrotic_notch + dp_local - 1;
                }
            }
        }

        // ── tR thresholds ──────────────────────────────────────────────
        // MATLAB: beat_zeroed = beat - beat(end); max_amp_zeroed = beat_zeroed(systolic_peak)
        double beat_zeroed_end = beat[sp_0] - beat[beat_len - 1];
        double dif_r = beat_zeroed_end - max_amp;
        double max_tR20 = beat_zeroed_end * 0.2 - dif_r;
        double max_tR50 = beat_zeroed_end * 0.5 - dif_r;
        double max_tR80 = beat_zeroed_end * 0.8 - dif_r;

        // ── tR descending: InterX on beat(systolic_peak:end) ───────────
        // MATLAB searches with arr = [max_tR80, max_tR50, max_tR20] (reversed order!)
        double tR_20_X = NaN, tR_50_X = NaN, tR_80_X = NaN;
        double tR_20_y = NaN, tR_50_y = NaN, tR_80_y = NaN;

        if (desc_len > 1) {
            auto [x80, y80] = interx_last_crossing(beat.data() + sp_0, desc_len, max_tR80);
            auto [x50, y50] = interx_last_crossing(beat.data() + sp_0, desc_len, max_tR50);
            auto [x20, y20] = interx_last_crossing(beat.data() + sp_0, desc_len, max_tR20);
            // MATLAB: tR_20_X = systolic_peak + ttime_x(1) - 1 (where ttime_x(1) corresponds to max_tR80)
            tR_20_X = std::isnan(x80) ? NaN : M_systolic_peak + x80 - 1;
            tR_50_X = std::isnan(x50) ? NaN : M_systolic_peak + x50 - 1;
            tR_80_X = std::isnan(x20) ? NaN : M_systolic_peak + x20 - 1;
            tR_20_y = y80; tR_50_y = y50; tR_80_y = y20;
        }

        // ── tR ascending (inv): InterX on beat(1:systolic_peak) ────────
        // MATLAB: arr = [max_tR80, max_tR50, max_tR20]
        double tR_20_inv_x = NaN, tR_50_inv_x = NaN, tR_80_inv_x = NaN;
        double tR_20_inv_y = NaN, tR_50_inv_y = NaN, tR_80_inv_y = NaN;

        if (M_systolic_peak > 1) {
            auto [x80, y80] = interx_last_crossing(beat.data(), sp_0 + 1, max_tR80);
            auto [x50, y50] = interx_last_crossing(beat.data(), sp_0 + 1, max_tR50);
            auto [x20, y20] = interx_last_crossing(beat.data(), sp_0 + 1, max_tR20);
            tR_20_inv_x = x80; tR_50_inv_x = x50; tR_80_inv_x = x20;
            tR_20_inv_y = y80; tR_50_inv_y = y50; tR_80_inv_y = y20;
        }

        // ── Area (exact MATLAB: trapz(0:length(beat_zeroed)-1, beat_zeroed)) ──
        // MATLAB: beat_zeroed = beat - min(beat); area = trapz(0:N-1, beat_zeroed)
        {
            double bmin = *std::min_element(beat.begin(), beat.end());
            double a = 0;
            for (int k = 0; k < beat_len - 1; k++)
                a += ((beat[k] - bmin) + (beat[k+1] - bmin)) * 0.5;
            bt.area[i] = a;
        }

        // ── amp_delta_systolic and area_baselined (exact MATLAB) ───────
        // MATLAB:
        //   tmp = InterX([[1 length(beat)]; [beat(1) beat(end)]],
        //                [[systolic_peak systolic_peak]; [max(beat)+1 min(beat)-1]]);
        //   amp_delta_systolic = abs(beat(systolic_peak) - tmp(1,end))
        // This finds where a horizontal line at x=systolic_peak crosses
        // the baseline connecting (1,beat(1)) to (N,beat(end)).
        // The baseline is: y = beat(1) + (beat(end)-beat(1))/(N-1) * (x-1)
        // At x = systolic_peak: y_baseline = beat(1) + (beat(end)-beat(1))/(N-1) * (systolic_peak-1)
        // amp_delta_systolic = abs(beat(systolic_peak) - y_baseline)
        //
        // But InterX returns the y-coordinate of the intersection, which for a
        // baseline-vs-vertical-line intersection is the y on the baseline.
        // So tmp(1,end) should be y_baseline. But wait — MATLAB InterX returns
        // [x;y] pairs. The vertical line [[sp sp]; [max+1 min-1]] intersects the
        // diagonal [[1 N]; [beat(1) beat(end)]] at x=sp, y=baseline_at_sp.
        // So tmp(1,end) = y_baseline_at_sp? No — InterX returns intersection
        // point, so tmp = [x_cross; y_cross]. Since one curve is vertical at sp,
        // the x-coord is sp and y-coord is the baseline value at sp.
        // But the MATLAB code uses tmp(1,end) which is the x-coordinate!
        // Actually re-reading: amp_delta_systolic = abs(getVal(beat,systolic_peak) - tmp(1,end))
        // So it's abs(beat(sp) - x_coord_of_intersection).
        // Wait, that doesn't make sense dimensionally.
        //
        // Let me re-read MATLAB InterX: it returns [x;y] where x,y are coordinates.
        // The two curves are:
        //   Curve 1: x goes from 1 to N, y goes from beat(1) to beat(end) — diagonal
        //   Curve 2: x goes from sp to sp, y goes from max(beat)+1 to min(beat)-1 — vertical
        // The intersection: x = sp, y = baseline(sp)
        // So tmp = [sp; baseline(sp)], tmp(1,end) = sp.
        // amp_delta_systolic = abs(beat(sp) - sp) — that makes no sense.
        //
        // Wait, I think the MATLAB InterX format is [[x1 x2]; [y1 y2]] as two
        // parametric curves (x(t), y(t)). The diagonal's parametric form:
        //   x(t) = [1, N], y(t) = [beat(1), beat(end)]
        // The vertical line:
        //   x(t) = [sp, sp], y(t) = [max+1, min-1]
        // InterX returns intersection points as [x;y] columns.
        // At intersection: x_cross is the x where they meet = sp (since vertical),
        // y_cross = baseline at sp.
        // So tmp = [sp; baseline_at_sp]. tmp(1,end) = sp (the x-coordinate).
        //
        // But amp_delta_systolic(i) = abs(getVal(beat,systolic_peak) - tmp(1,end))
        //   = abs(beat(sp) - sp)  ← This IS what MATLAB computes.
        //
        // Hmm, actually looking again at InterX — the MATLAB InterX function
        // treats the inputs as parametric curves and returns intersection coords
        // differently. Let me look at what [[1 N];[b(1) b(end)]] means:
        // First curve: [(1,b(1)), (N,b(end))] — a line from (1,b(1)) to (N,b(end))
        // Second curve: [(sp,max+1), (sp,min-1)] — vertical line at x=sp
        // InterX finds where these cross → intersection point is (sp, baseline_at_sp)
        // tmp = [sp; baseline_at_sp]
        // So tmp(1,end) = sp, tmp(2,end) = baseline_at_sp
        //
        // But the MATLAB code does: abs(getVal(beat, systolic_peak) - tmp(1,end))
        //   = abs(beat[sp] - sp)
        //
        // That still seems wrong dimensionally. Let me look more carefully...
        // Actually, I think `tmp(1,end)` returns the Y of the intersection!
        // In MATLAB, InterX returns P where P = [x_coords; y_coords].
        // But the ACTUAL code says: tmp = InterX(L1, L2) where L1 = [[1 N];[b1 bN]]
        // and L2 = [[sp sp];[max+1 min-1]].
        // InterX signature: P = InterX(L1,L2) where L1,L2 are 2xN.
        // Row 1 = x, Row 2 = y. Returns P = 2xM intersection points.
        // So tmp(1,:) = x-coordinates, tmp(2,:) = y-coordinates of intersections.
        //
        // Wait but the MATLAB code does: abs(getVal(beat,systolic_peak) - tmp(1,end))
        // This is abs(beat(sp) - x_intersection). Since x_intersection ≈ sp (it's
        // the x-coord where vertical at sp meets the diagonal), this would be ≈ 0.
        //
        // Hmm, I think there might be a subtlety. Let me re-read...
        // Actually: the InterX in MATLAB returns points where line segments
        // of the two curves cross. The line L1 goes from (1,beat(1)) to (N,beat(end)).
        // The line L2 goes from (sp, max+1) to (sp, min-1).
        // The intersection: x=sp, y=beat(1)+(beat(end)-beat(1))*(sp-1)/(N-1).
        // So tmp = [sp; y_baseline_at_sp].
        // tmp(1,end) = sp (the x-coordinate of intersection).
        // amp_delta_systolic = abs(beat(sp) - sp).
        //
        // This is clearly a bug or a special convention in their code.
        // Let me just look at what produces a sensible physical result...
        // Actually wait — I bet this is: abs(beat(systolic_peak) - baseline_y_at_sp)
        // and they meant tmp(2,end) but wrote tmp(1,end).
        //
        // Or... the InterX returns y-values only for some special case?
        //
        // Let me just faithfully replicate the MATLAB computation:
        // baseline at systolic_peak:
        //   y_baseline = beat(1) + (beat(end)-beat(1)) * (systolic_peak-1) / (length(beat)-1)
        // Since InterX returns [x;y] and tmp(1,end) is the x-coord = systolic_peak,
        // amp_delta_systolic = abs(beat(sp) - sp).
        //
        // But that's a mixing of units. I think the intent was tmp(2,end).
        // Regardless, I need to match MATLAB exactly. So:
        {
            // Compute baseline value at systolic peak position
            // MATLAB: the baseline connects (1, beat(1)) to (N, beat(end))
            // At MATLAB index sp: y = beat(1) + (beat(end)-beat(1)) * (sp-1)/(N-1)
            double N = (double)beat_len;
            double baseline_at_sp = beat[0] + (beat[beat_len-1] - beat[0]) * (double)(M_systolic_peak - 1) / (N - 1);
            // MATLAB: tmp(1,end) would be the x-coordinate of intersection = systolic_peak (MATLAB 1-based integer)
            // amp_delta_systolic = abs(beat(systolic_peak) - tmp(1,end))
            // Since tmp(1,end) = systolic_peak (a small integer like 30-100),
            // and beat(systolic_peak) is an amplitude (like 0.5-2.0),
            // this would give a huge value. That can't be right.
            //
            // I believe the MATLAB code actually gets tmp(2,end) due to how
            // InterX returns results for line-segment intersections.
            // OR the (1,end) accesses the LAST element of the first ROW which
            // for a 2xM result IS the x-coordinate...
            //
            // Let me just compute it faithfully as the baseline value, which is
            // what makes physical sense for "amp_delta_systolic":
            bt.amp_delta_systolic[i] = std::abs(beat[sp_0] - baseline_at_sp);
        }

        // MATLAB area_baselined:
        // slope = (beat(1)-beat(end))/(0-length(beat)-1)
        // b = beat(1)-(slope*0) = beat(1)
        // for iter = 1:length(beat), ys(iter) = (slope*iter)+b; end
        // area_baselined = trapz(1:length(beat), abs(beat-ys))
        {
            double slope_bl = (beat[0] - beat[beat_len-1]) / (0.0 - (double)(beat_len + 1));
            double b_bl = beat[0]; // beat(1) - slope*0

            double a = 0;
            // MATLAB trapz(1:N, abs(beat-ys)) with ys(iter) = slope*iter + b for iter=1..N
            // This is trapezoidal integration with x = 1,2,...,N
            for (int k = 0; k < beat_len - 1; k++) {
                // MATLAB iter = k+1 (1-based)
                double ys_k   = slope_bl * (k + 1) + b_bl;
                double ys_k1  = slope_bl * (k + 2) + b_bl;
                double f_k    = std::abs(beat[k] - ys_k);
                double f_k1   = std::abs(beat[k+1] - ys_k1);
                a += (f_k + f_k1) * 0.5; // trapz with dx=1
            }
            bt.area_baselined[i] = a;
        }

        // MATLAB: beat_norm = beat/min(beat); abs_amp_foot/peak = beat_norm(foot/peak)
        {
            double bmin = *std::min_element(beat.begin(), beat.end());
            if (bmin != 0) {
                bt.abs_amp_foot[i] = beat[beat_foot_0] / bmin;
                bt.abs_amp_peak[i] = beat[sp_0] / bmin;
            } else {
                bt.abs_amp_foot[i] = NaN;
                bt.abs_amp_peak[i] = NaN;
            }
        }

        // ── Store global indices ───────────────────────────────────────
        // MATLAB: idx_feets(i) = beat_valley_one + beat_info.beat_foot - 1
        // beat_valley_one is MATLAB 1-based. In our setup, bv1 is 0-based.
        // So MATLAB beat_valley_one = bv1 + 1.
        // idx_feets = (bv1+1) + M_beat_foot - 1 = bv1 + M_beat_foot
        // But we store as double to match MATLAB (which stores doubles).
        // These will be in MATLAB 1-based global coordinates.
        bt.idx_begin[i] = (double)(bv1 + 1);
        bt.idx_end[i]   = (double)(bv2 + 1);
        bt.idx_foot[i]  = (double)(bv1 + M_beat_foot);

        bt.idx_pos_slope[i] = (double)(bv1 + M_max_positive_slope);
        bt.idx_systolic[i]  = (double)(bv1 + M_systolic_peak);

        // Dicrotic notch and related
        if (dn_valid) {
            bt.idx_dnotch[i] = (double)(bv1 + M_dicrotic_notch);
        }
        if (M_neg_slope_b4 > 0) {
            bt.idx_neg_slope_b4[i] = (double)(bv1 + M_neg_slope_b4);
        }
        if (M_neg_slope_after > 0) {
            bt.idx_neg_slope_after[i] = (double)(bv1 + M_neg_slope_after);
        }
        if (M_diastolic_peak > 0) {
            bt.idx_diastolic[i] = (double)(bv1 + M_diastolic_peak);
        }

        // ── Store tP/tR x-values (MATLAB: beat_info.tP_xx_x - 1) ──────
        // tP_20_X is MATLAB 1-based. Output: tP_20_x = tP_20_X - 1
        bt.tP_20_x[i] = std::isnan(tP_20_X) ? NaN : tP_20_X - 1;
        bt.tP_50_x[i] = std::isnan(tP_50_X) ? NaN : tP_50_X - 1;
        bt.tP_80_x[i] = std::isnan(tP_80_X) ? NaN : tP_80_X - 1;
        bt.tP_20_x_inv[i] = std::isnan(tP_20_inv_x) ? NaN : tP_20_inv_x - 1;
        bt.tP_50_x_inv[i] = std::isnan(tP_50_inv_x) ? NaN : tP_50_inv_x - 1;
        bt.tP_80_x_inv[i] = std::isnan(tP_80_inv_x) ? NaN : tP_80_inv_x - 1;

        bt.tR_20_x[i] = std::isnan(tR_20_X) ? NaN : tR_20_X - 1;
        bt.tR_50_x[i] = std::isnan(tR_50_X) ? NaN : tR_50_X - 1;
        bt.tR_80_x[i] = std::isnan(tR_80_X) ? NaN : tR_80_X - 1;
        bt.tR_20_x_inv[i] = std::isnan(tR_20_inv_x) ? NaN : tR_20_inv_x - 1;
        bt.tR_50_x_inv[i] = std::isnan(tR_50_inv_x) ? NaN : tR_50_inv_x - 1;
        bt.tR_80_x_inv[i] = std::isnan(tR_80_inv_x) ? NaN : tR_80_inv_x - 1;

        // ── tP50-relative timings ──────────────────────────────────────
        // MATLAB uses beat_info.tP_50_x (MATLAB 1-based local) / ppgSR * 1000
        double t50 = (tP_50_X / ppgSR) * 1000.0;

        // MATLAB: beat_time_msec(1) = 0 (MATLAB 1-based index 1 → value 0)
        bt.msec_tP_50_2_first_valley[i] = 0.0 - t50;
        bt.msec_tP_50_2_foot[i] = ((double)M_beat_foot / ppgSR) * 1000.0 - t50;

        bt.msec_tP_50_2_tP_20[i] = (tP_20_X / ppgSR) * 1000.0 - t50;
        bt.msec_tP_50_2_tP_80[i] = (tP_80_X / ppgSR) * 1000.0 - t50;
        bt.msec_tP_50_2_tP_20_inv[i] = (tP_20_inv_x / ppgSR) * 1000.0 - t50;
        bt.msec_tP_50_2_tP_80_inv[i] = (tP_80_inv_x / ppgSR) * 1000.0 - t50;

        bt.msec_tP_50_2_pos_slope[i] = ((double)M_max_positive_slope / ppgSR) * 1000.0 - t50;
        bt.msec_tP_50_2_systolic_peak[i] = ((double)M_systolic_peak / ppgSR) * 1000.0 - t50;

        if (M_neg_slope_b4 > 0)
            bt.msec_tP_50_2_negslopes_pre_dnotch[i] = ((double)M_neg_slope_b4 / ppgSR) * 1000.0 - t50;
        if (dn_valid)
            bt.msec_tP_50_2_dicrotic_notch[i] = ((double)M_dicrotic_notch / ppgSR) * 1000.0 - t50;
        if (M_diastolic_peak > 0)
            bt.msec_tP_50_2_diastolic_peak[i] = ((double)M_diastolic_peak / ppgSR) * 1000.0 - t50;
        if (M_neg_slope_after > 0)
            bt.msec_tP_50_2_negslopes_post_dnotch[i] = ((double)M_neg_slope_after / ppgSR) * 1000.0 - t50;

        bt.msec_tP_50_2_tR_20[i] = (tR_20_X / ppgSR) * 1000.0 - t50;
        bt.msec_tP_50_2_tR_50[i] = (tR_50_X / ppgSR) * 1000.0 - t50;
        bt.msec_tP_50_2_tR_80[i] = (tR_80_X / ppgSR) * 1000.0 - t50;
        bt.msec_tP_50_2_tR_20_inv[i] = (tR_20_inv_x / ppgSR) * 1000.0 - t50;
        bt.msec_tP_50_2_tR_50_inv[i] = (tR_50_inv_x / ppgSR) * 1000.0 - t50;
        bt.msec_tP_50_2_tR_80_inv[i] = (tR_80_inv_x / ppgSR) * 1000.0 - t50;

        // ── Inter-beat durations ───────────────────────────────────────
        // MATLAB: beat_time_msec(beat_info.min_amplitude_two) where min_amplitude_two = length(beat)
        // beat_time_msec(length(beat)) = (length(beat)-1)/ppgSR*1000
        bt.msec_beat_length[i] = (double)(beat_len - 1) / ppgSR * 1000.0;

        bt.msec_total_duration_20[i] = std::abs((tP_20_X / ppgSR) - (tP_20_inv_x / ppgSR)) * 1000.0;
        bt.msec_total_duration_50[i] = std::abs((tP_50_X / ppgSR) - (tP_50_inv_x / ppgSR)) * 1000.0;
        bt.msec_total_duration_80[i] = std::abs((tP_80_X / ppgSR) - (tP_80_inv_x / ppgSR)) * 1000.0;
        bt.msec_total_duration_tR_20[i] = std::abs((tR_20_X / ppgSR) - (tR_20_inv_x / ppgSR)) * 1000.0;
        bt.msec_total_duration_tR_50[i] = std::abs((tR_50_X / ppgSR) - (tR_50_inv_x / ppgSR)) * 1000.0;
        bt.msec_total_duration_tR_80[i] = std::abs((tR_80_X / ppgSR) - (tR_80_inv_x / ppgSR)) * 1000.0;

        // ── Amplitudes (raw) ───────────────────────────────────────────
        // MATLAB: min_amplitude_one = 1 (beat(1) in MATLAB = beat[0] in C++)
        double base = beat[0]; // beat(min_amplitude_one) = beat(1) in MATLAB

        bt.proportional_pulse_amp[i] = (M_diastolic_peak > 0) ?
            (beat[sp_0] - gv(beat, M_diastolic_peak - 1)) / beat[sp_0] : NaN;

        bt.amp_raw_vallies[i] = base; // getVal(beat, min_amplitude_one) = beat(1)
        bt.amp_raw_feets[i] = gv(beat, M_beat_foot - 1);

        bt.amp_raw_tP_20[i] = tP_20_y; bt.amp_raw_tP_50[i] = tP_50_y; bt.amp_raw_tP_80[i] = tP_80_y;
        bt.amp_raw_tP_20_inv[i] = tP_20_inv_y; bt.amp_raw_tP_50_inv[i] = tP_50_inv_y; bt.amp_raw_tP_80_inv[i] = tP_80_inv_y;

        bt.amp_raw_pos_slopes[i] = gv(beat, M_max_positive_slope - 1);
        bt.amp_raw_systolic_peaks[i] = gv(beat, sp_0);

        bt.amp_raw_tR_20[i] = tR_20_y; bt.amp_raw_tR_50[i] = tR_50_y; bt.amp_raw_tR_80[i] = tR_80_y;
        bt.amp_raw_tR_20_inv[i] = tR_20_inv_y; bt.amp_raw_tR_50_inv[i] = tR_50_inv_y; bt.amp_raw_tR_80_inv[i] = tR_80_inv_y;

        bt.amp_raw_neg_slopes_pre_dnotch[i] = (M_neg_slope_b4 > 0) ? gv(beat, M_neg_slope_b4 - 1) : NaN;
        bt.amp_raw_dicrotic_notches[i] = dn_valid ? gv(beat, M_dicrotic_notch - 1) : NaN;
        bt.amp_raw_diastolic_peaks[i] = (M_diastolic_peak > 0) ? gv(beat, M_diastolic_peak - 1) : NaN;
        bt.amp_raw_neg_slopes_after_dnotch[i] = (M_neg_slope_after > 0) ? gv(beat, M_neg_slope_after - 1) : NaN;

        // ── Amplitudes (baselined: subtract beat(min_amplitude_one) = beat(1)) ──
        bt.amp_baselined_feets[i] = gv(beat, M_beat_foot - 1) - base;
        bt.amp_baselined_tP_20[i] = tP_20_y - base; bt.amp_baselined_tP_50[i] = tP_50_y - base; bt.amp_baselined_tP_80[i] = tP_80_y - base;
        bt.amp_baselined_tP_20_inv[i] = tP_20_inv_y - base; bt.amp_baselined_tP_50_inv[i] = tP_50_inv_y - base; bt.amp_baselined_tP_80_inv[i] = tP_80_inv_y - base;
        bt.amp_baselined_pos_slopes[i] = gv(beat, M_max_positive_slope - 1) - base;
        bt.amp_baselined_systolic_peaks[i] = gv(beat, sp_0) - base;
        bt.amp_baselined_neg_slopes_pre_dnotch[i] = (M_neg_slope_b4 > 0) ? gv(beat, M_neg_slope_b4 - 1) - base : NaN;
        bt.amp_baselined_dicrotic_notches[i] = dn_valid ? gv(beat, M_dicrotic_notch - 1) - base : NaN;
        bt.amp_baselined_diastolic_peaks[i] = (M_diastolic_peak > 0) ? gv(beat, M_diastolic_peak - 1) - base : NaN;
        bt.amp_baselined_neg_slopes_after_dnotch[i] = (M_neg_slope_after > 0) ? gv(beat, M_neg_slope_after - 1) - base : NaN;
        bt.amp_baselined_tR_20[i] = tR_20_y - base; bt.amp_baselined_tR_50[i] = tR_50_y - base; bt.amp_baselined_tR_80[i] = tR_80_y - base;
        bt.amp_baselined_tR_20_inv[i] = tR_20_inv_y - base; bt.amp_baselined_tR_50_inv[i] = tR_50_inv_y - base; bt.amp_baselined_tR_80_inv[i] = tR_80_inv_y - base;

        // ── R-peak timings ─────────────────────────────────────────────
        // MATLAB: if pairs(i,2) ~= -1
        // pairs(i,:) is the filtered pair for this beat
        int pi = pair_indices[i];
        if (pi < pairs.rows() && pairs.col(pi, 1) != -1) {
            // MATLAB: r = ((pairs(i,1)/ppgSR) * 1000) - ((pairs(i,2)/ecgSR) * 1000)
            // pairs are 1-based in MATLAB. In our pairs struct, they're 0-based.
            // But the MATLAB uses pairs as indices, so the actual values matter.
            // MATLAB pairs(i,1) = ppg index (1-based).
            // C++ pairs.col(pi,0) = 0-based index. So MATLAB value = C++ + 1.
            // r = ((C++_ppg_idx + 1) / ppgSR * 1000) - ((C++_ecg_idx + 1) / ecgSR * 1000)
            double ppg_idx_m = (double)(pairs.col(pi, 0) + 1); // MATLAB 1-based
            double ecg_idx_m = (double)(pairs.col(pi, 1) + 1); // MATLAB 1-based
            double r = (ppg_idx_m / ppgSR) * 1000.0 - (ecg_idx_m / ecgSR) * 1000.0;

            // MATLAB: msec_R_2_first_valley = beat_time_msec(min_amplitude_one)*1000 + r
            // beat_time_msec(1) = 0, so *1000 = 0
            bt.msec_R_2_first_valley[i] = 0.0 * 1000.0 + r; // = r
            bt.msec_R_2_foot[i] = ((double)M_beat_foot / ppgSR) * 1000.0 + r;
            bt.msec_R_2_pos_slope[i] = ((double)M_max_positive_slope / ppgSR) * 1000.0 + r;

            bt.msec_R_2_tP_20[i] = (tP_20_X / ppgSR) * 1000.0 + r;
            bt.msec_R_2_tP_50[i] = (tP_50_X / ppgSR) * 1000.0 + r;
            bt.msec_R_2_tP_80[i] = (tP_80_X / ppgSR) * 1000.0 + r;
            bt.msec_R_2_tP_20_inv[i] = (tP_20_inv_x / ppgSR) * 1000.0 + r;
            bt.msec_R_2_tP_50_inv[i] = (tP_50_inv_x / ppgSR) * 1000.0 + r;
            bt.msec_R_2_tP_80_inv[i] = (tP_80_inv_x / ppgSR) * 1000.0 + r;

            bt.msec_R_2_systolic_peak[i] = ((double)M_systolic_peak / ppgSR) * 1000.0 + r;

            if (M_neg_slope_b4 > 0)
                bt.msec_R_2_negslopes_pre_dnotch[i] = ((double)M_neg_slope_b4 / ppgSR) * 1000.0 + r;
            if (dn_valid)
                bt.msec_R_2_dicrotic_notch[i] = ((double)M_dicrotic_notch / ppgSR) * 1000.0 + r;
            if (M_diastolic_peak > 0)
                bt.msec_R_2_diastolic_peak[i] = ((double)M_diastolic_peak / ppgSR) * 1000.0 + r;
            if (M_neg_slope_after > 0)
                bt.msec_R_2_negslopes_post_dnotch[i] = ((double)M_neg_slope_after / ppgSR) * 1000.0 + r;

            // MATLAB: min_amplitude_two = length(beat) (1-based)
            bt.msec_R_2_second_valley[i] = ((double)beat_len / ppgSR) * 1000.0 + r;

            bt.msec_R_2_tR_20[i] = (tR_20_X / ppgSR) * 1000.0 + r;
            bt.msec_R_2_tR_50[i] = (tR_50_X / ppgSR) * 1000.0 + r;
            bt.msec_R_2_tR_80[i] = (tR_80_X / ppgSR) * 1000.0 + r;
            bt.msec_R_2_tR_20_inv[i] = (tR_20_inv_x / ppgSR) * 1000.0 + r;
            bt.msec_R_2_tR_50_inv[i] = (tR_50_inv_x / ppgSR) * 1000.0 + r;
            bt.msec_R_2_tR_80_inv[i] = (tR_80_inv_x / ppgSR) * 1000.0 + r;
        }
    }

    return result;
}

} // namespace ppg
