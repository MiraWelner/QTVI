#pragma once
// ═══════════════════════════════════════════════════════════════════════════════
// flatten_and_generate.hpp — Exact match of flatten_beat_idx.m & GenerateFeatures.m
//
// Key fixes vs old C++:
//  - ALL fields from BeatsTotal are copied to BeatsFlattened (not just a subset)
//  - ALL inter-beat intervals computed (sec_tP_20_2_tP_20, etc.)
//  - No std::map — all flat struct members
//  - Proper sleep state logic matching MATLAB exactly
//  - Correct time correction from ppg_bin_indexs
// ═══════════════════════════════════════════════════════════════════════════════
#include "ppg_features.hpp"
#include "bin_io.hpp"

namespace ppg {

    // ─── Count Bins ─────────────────────────────────────────────────────────────
    inline BinMarksMask countBins(const std::vector<WaveData>& waveData,
        const std::vector<TemplateInfo>* tmplate) {
        BinMarksMask mask;
        int n = (int)waveData.size();
        mask.could_not_identify_PPG.resize(n, false);
        mask.poor_ppg_or_ecg_template_manually_excluded.resize(n, false);
        mask.ecg_template_manually_excluded.resize(n, false);
        mask.ppg_template_manually_excluded.resize(n, false);
        for (int i = 0; i < n; i++)
            mask.could_not_identify_PPG[i] = (bool)waveData[i].bad_segment;
        if (tmplate) {
            for (int i = 0; i < n && i < (int)tmplate->size(); i++) {
                mask.poor_ppg_or_ecg_template_manually_excluded[i] = (bool)(*tmplate)[i].TemplateBad;
                mask.ecg_template_manually_excluded[i] = (bool)(*tmplate)[i].bad_r_templates;
                mask.ppg_template_manually_excluded[i] = (bool)(*tmplate)[i].bad_ppg_templates;
            }
        }
        return mask;
    }

    // ─── Flatten Beat Indices (exact flatten_beat_idx.m) ────────────────────────
    inline BeatsFlattened flatten_beat_idx(
        const std::vector<BeatsTotal>& bins,
        const std::vector<AnnealedSegment>& processSegments)
    {
        BeatsFlattened flat;

        // Count total beats
        int total_beats = 0;
        for (auto& b : bins) total_beats += (int)b.area.size();
        if (total_beats == 0) return flat;

        double ppgRate = processSegments[0].ppgSampleRate;

        // Initialize all vectors to NaN
        auto init = [&](std::vector<double>& v) { v.assign(total_beats, NaN); };

        // Index fields (need start offset added)
        init(flat.idx_begin); init(flat.idx_end); init(flat.idx_foot);
        init(flat.idx_pos_slope); init(flat.idx_systolic);
        init(flat.idx_neg_slope_b4); init(flat.idx_neg_slope_after);
        init(flat.idx_diastolic); init(flat.idx_dnotch);

        // tP/tR x-values (not idx, so no offset)
        init(flat.tP_20_x); init(flat.tP_50_x); init(flat.tP_80_x);
        init(flat.tP_20_x_inv); init(flat.tP_50_x_inv); init(flat.tP_80_x_inv);
        init(flat.tR_20_x); init(flat.tR_50_x); init(flat.tR_80_x);
        init(flat.tR_20_x_inv); init(flat.tR_50_x_inv); init(flat.tR_80_x_inv);

        // Scalar fields
        init(flat.sleep_stages); init(flat.area_baselined); init(flat.area);
        init(flat.amp_delta_systolic); init(flat.abs_amp_foot); init(flat.abs_amp_peak);
        init(flat.msec_beat_length);

        // All timing fields
        init(flat.msec_tP_50_2_first_valley); init(flat.msec_tP_50_2_foot);
        init(flat.msec_tP_50_2_tP_20); init(flat.msec_tP_50_2_tP_80);
        init(flat.msec_tP_50_2_tP_20_inv); init(flat.msec_tP_50_2_tP_80_inv);
        init(flat.msec_tP_50_2_pos_slope); init(flat.msec_tP_50_2_systolic_peak);
        init(flat.msec_tP_50_2_negslopes_pre_dnotch); init(flat.msec_tP_50_2_dicrotic_notch);
        init(flat.msec_tP_50_2_diastolic_peak); init(flat.msec_tP_50_2_negslopes_post_dnotch);
        init(flat.msec_tP_50_2_second_valley);
        init(flat.msec_tP_50_2_tR_20); init(flat.msec_tP_50_2_tR_50); init(flat.msec_tP_50_2_tR_80);
        init(flat.msec_tP_50_2_tR_20_inv); init(flat.msec_tP_50_2_tR_50_inv); init(flat.msec_tP_50_2_tR_80_inv);
        init(flat.msec_total_duration_20); init(flat.msec_total_duration_50); init(flat.msec_total_duration_80);
        init(flat.msec_total_duration_tR_20); init(flat.msec_total_duration_tR_50); init(flat.msec_total_duration_tR_80);

        // R-peak timings
        init(flat.msec_R_2_first_valley); init(flat.msec_R_2_foot);
        init(flat.msec_R_2_tP_20); init(flat.msec_R_2_tP_50); init(flat.msec_R_2_tP_80);
        init(flat.msec_R_2_tP_20_inv); init(flat.msec_R_2_tP_50_inv); init(flat.msec_R_2_tP_80_inv);
        init(flat.msec_R_2_pos_slope); init(flat.msec_R_2_systolic_peak);
        init(flat.msec_R_2_negslopes_pre_dnotch); init(flat.msec_R_2_dicrotic_notch);
        init(flat.msec_R_2_tR_20); init(flat.msec_R_2_tR_50); init(flat.msec_R_2_tR_80);
        init(flat.msec_R_2_tR_20_inv); init(flat.msec_R_2_tR_50_inv); init(flat.msec_R_2_tR_80_inv);
        init(flat.msec_R_2_diastolic_peak); init(flat.msec_R_2_negslopes_post_dnotch);
        init(flat.msec_R_2_second_valley);

        // Amplitudes
        init(flat.amp_raw_vallies); init(flat.amp_raw_feets);
        init(flat.amp_raw_tP_20); init(flat.amp_raw_tP_50); init(flat.amp_raw_tP_80);
        init(flat.amp_raw_tP_20_inv); init(flat.amp_raw_tP_50_inv); init(flat.amp_raw_tP_80_inv);
        init(flat.amp_raw_tR_20); init(flat.amp_raw_tR_50); init(flat.amp_raw_tR_80);
        init(flat.amp_raw_tR_20_inv); init(flat.amp_raw_tR_50_inv); init(flat.amp_raw_tR_80_inv);
        init(flat.amp_raw_pos_slopes); init(flat.amp_raw_systolic_peaks);
        init(flat.amp_raw_neg_slopes_pre_dnotch); init(flat.amp_raw_dicrotic_notches);
        init(flat.amp_raw_diastolic_peaks); init(flat.amp_raw_neg_slopes_after_dnotch);
        init(flat.amp_baselined_feets);
        init(flat.amp_baselined_tP_20); init(flat.amp_baselined_tP_50); init(flat.amp_baselined_tP_80);
        init(flat.amp_baselined_tP_20_inv); init(flat.amp_baselined_tP_50_inv); init(flat.amp_baselined_tP_80_inv);
        init(flat.amp_baselined_tR_20); init(flat.amp_baselined_tR_50); init(flat.amp_baselined_tR_80);
        init(flat.amp_baselined_tR_20_inv); init(flat.amp_baselined_tR_50_inv); init(flat.amp_baselined_tR_80_inv);
        init(flat.amp_baselined_pos_slopes); init(flat.amp_baselined_systolic_peaks);
        init(flat.amp_baselined_neg_slopes_pre_dnotch); init(flat.amp_baselined_dicrotic_notches);
        init(flat.amp_baselined_diastolic_peaks); init(flat.amp_baselined_neg_slopes_after_dnotch);
        init(flat.proportional_pulse_amp);

        // SQI
        init(flat.sqi_mean_corr_dtw); init(flat.sqi_corrcoff_direct);
        init(flat.sqi_corrcoff_interp); init(flat.sqi_dtw); init(flat.sqi_frechet);

        // ── Copy from bins to flat ──────────────────────────────────────────
        int prev = 0;
        double start = 0; // MATLAB: cumulative offset for idx fields

        for (int x = 0; x < (int)bins.size(); x++) {
            int len = (int)bins[x].area.size();
            if (x > 0) start += (double)processSegments[x - 1].po.size() - 1;
            if (len == 0) continue;

            // Helper: copy idx field (add offset)
            auto ci = [&](std::vector<double>& dst, const std::vector<double>& src) {
                for (int i = 0; i < len && i < (int)src.size(); i++)
                    dst[prev + i] = src[i] + start;
                };
            // Helper: copy non-idx field
            auto cd = [&](std::vector<double>& dst, const std::vector<double>& src) {
                for (int i = 0; i < len && i < (int)src.size(); i++)
                    dst[prev + i] = src[i];
                };

            // Index fields (add start offset)
            ci(flat.idx_begin, bins[x].idx_begin);
            ci(flat.idx_end, bins[x].idx_end);
            ci(flat.idx_foot, bins[x].idx_foot);
            ci(flat.idx_pos_slope, bins[x].idx_pos_slope);
            ci(flat.idx_systolic, bins[x].idx_systolic);
            ci(flat.idx_neg_slope_b4, bins[x].idx_neg_slope_b4);
            ci(flat.idx_neg_slope_after, bins[x].idx_neg_slope_after);
            ci(flat.idx_diastolic, bins[x].idx_diastolic);
            ci(flat.idx_dnotch, bins[x].idx_dnotch);

            // Non-idx (copy as-is)
            cd(flat.tP_20_x, bins[x].tP_20_x); cd(flat.tP_50_x, bins[x].tP_50_x); cd(flat.tP_80_x, bins[x].tP_80_x);
            cd(flat.tP_20_x_inv, bins[x].tP_20_x_inv); cd(flat.tP_50_x_inv, bins[x].tP_50_x_inv); cd(flat.tP_80_x_inv, bins[x].tP_80_x_inv);
            cd(flat.tR_20_x, bins[x].tR_20_x); cd(flat.tR_50_x, bins[x].tR_50_x); cd(flat.tR_80_x, bins[x].tR_80_x);
            cd(flat.tR_20_x_inv, bins[x].tR_20_x_inv); cd(flat.tR_50_x_inv, bins[x].tR_50_x_inv); cd(flat.tR_80_x_inv, bins[x].tR_80_x_inv);

            cd(flat.sleep_stages, bins[x].sleep_stages);
            cd(flat.area_baselined, bins[x].area_baselined); cd(flat.area, bins[x].area);
            cd(flat.amp_delta_systolic, bins[x].amp_delta_systolic);
            cd(flat.abs_amp_foot, bins[x].abs_amp_foot); cd(flat.abs_amp_peak, bins[x].abs_amp_peak);
            cd(flat.msec_beat_length, bins[x].msec_beat_length);

            cd(flat.msec_tP_50_2_first_valley, bins[x].msec_tP_50_2_first_valley);
            cd(flat.msec_tP_50_2_foot, bins[x].msec_tP_50_2_foot);
            cd(flat.msec_tP_50_2_tP_20, bins[x].msec_tP_50_2_tP_20);
            cd(flat.msec_tP_50_2_tP_80, bins[x].msec_tP_50_2_tP_80);
            cd(flat.msec_tP_50_2_tP_20_inv, bins[x].msec_tP_50_2_tP_20_inv);
            cd(flat.msec_tP_50_2_tP_80_inv, bins[x].msec_tP_50_2_tP_80_inv);
            cd(flat.msec_tP_50_2_pos_slope, bins[x].msec_tP_50_2_pos_slope);
            cd(flat.msec_tP_50_2_systolic_peak, bins[x].msec_tP_50_2_systolic_peak);
            cd(flat.msec_tP_50_2_negslopes_pre_dnotch, bins[x].msec_tP_50_2_negslopes_pre_dnotch);
            cd(flat.msec_tP_50_2_dicrotic_notch, bins[x].msec_tP_50_2_dicrotic_notch);
            cd(flat.msec_tP_50_2_diastolic_peak, bins[x].msec_tP_50_2_diastolic_peak);
            cd(flat.msec_tP_50_2_negslopes_post_dnotch, bins[x].msec_tP_50_2_negslopes_post_dnotch);
            cd(flat.msec_tP_50_2_second_valley, bins[x].msec_tP_50_2_second_valley);
            cd(flat.msec_tP_50_2_tR_20, bins[x].msec_tP_50_2_tR_20);
            cd(flat.msec_tP_50_2_tR_50, bins[x].msec_tP_50_2_tR_50);
            cd(flat.msec_tP_50_2_tR_80, bins[x].msec_tP_50_2_tR_80);
            cd(flat.msec_tP_50_2_tR_20_inv, bins[x].msec_tP_50_2_tR_20_inv);
            cd(flat.msec_tP_50_2_tR_50_inv, bins[x].msec_tP_50_2_tR_50_inv);
            cd(flat.msec_tP_50_2_tR_80_inv, bins[x].msec_tP_50_2_tR_80_inv);
            cd(flat.msec_total_duration_20, bins[x].msec_total_duration_20);
            cd(flat.msec_total_duration_50, bins[x].msec_total_duration_50);
            cd(flat.msec_total_duration_80, bins[x].msec_total_duration_80);
            cd(flat.msec_total_duration_tR_20, bins[x].msec_total_duration_tR_20);
            cd(flat.msec_total_duration_tR_50, bins[x].msec_total_duration_tR_50);
            cd(flat.msec_total_duration_tR_80, bins[x].msec_total_duration_tR_80);

            cd(flat.msec_R_2_first_valley, bins[x].msec_R_2_first_valley);
            cd(flat.msec_R_2_foot, bins[x].msec_R_2_foot);
            cd(flat.msec_R_2_tP_20, bins[x].msec_R_2_tP_20); cd(flat.msec_R_2_tP_50, bins[x].msec_R_2_tP_50); cd(flat.msec_R_2_tP_80, bins[x].msec_R_2_tP_80);
            cd(flat.msec_R_2_tP_20_inv, bins[x].msec_R_2_tP_20_inv); cd(flat.msec_R_2_tP_50_inv, bins[x].msec_R_2_tP_50_inv); cd(flat.msec_R_2_tP_80_inv, bins[x].msec_R_2_tP_80_inv);
            cd(flat.msec_R_2_pos_slope, bins[x].msec_R_2_pos_slope);
            cd(flat.msec_R_2_systolic_peak, bins[x].msec_R_2_systolic_peak);
            cd(flat.msec_R_2_negslopes_pre_dnotch, bins[x].msec_R_2_negslopes_pre_dnotch);
            cd(flat.msec_R_2_dicrotic_notch, bins[x].msec_R_2_dicrotic_notch);
            cd(flat.msec_R_2_tR_20, bins[x].msec_R_2_tR_20); cd(flat.msec_R_2_tR_50, bins[x].msec_R_2_tR_50); cd(flat.msec_R_2_tR_80, bins[x].msec_R_2_tR_80);
            cd(flat.msec_R_2_tR_20_inv, bins[x].msec_R_2_tR_20_inv); cd(flat.msec_R_2_tR_50_inv, bins[x].msec_R_2_tR_50_inv); cd(flat.msec_R_2_tR_80_inv, bins[x].msec_R_2_tR_80_inv);
            cd(flat.msec_R_2_diastolic_peak, bins[x].msec_R_2_diastolic_peak);
            cd(flat.msec_R_2_negslopes_post_dnotch, bins[x].msec_R_2_negslopes_post_dnotch);
            cd(flat.msec_R_2_second_valley, bins[x].msec_R_2_second_valley);

            cd(flat.amp_raw_vallies, bins[x].amp_raw_vallies); cd(flat.amp_raw_feets, bins[x].amp_raw_feets);
            cd(flat.amp_raw_tP_20, bins[x].amp_raw_tP_20); cd(flat.amp_raw_tP_50, bins[x].amp_raw_tP_50); cd(flat.amp_raw_tP_80, bins[x].amp_raw_tP_80);
            cd(flat.amp_raw_tP_20_inv, bins[x].amp_raw_tP_20_inv); cd(flat.amp_raw_tP_50_inv, bins[x].amp_raw_tP_50_inv); cd(flat.amp_raw_tP_80_inv, bins[x].amp_raw_tP_80_inv);
            cd(flat.amp_raw_tR_20, bins[x].amp_raw_tR_20); cd(flat.amp_raw_tR_50, bins[x].amp_raw_tR_50); cd(flat.amp_raw_tR_80, bins[x].amp_raw_tR_80);
            cd(flat.amp_raw_tR_20_inv, bins[x].amp_raw_tR_20_inv); cd(flat.amp_raw_tR_50_inv, bins[x].amp_raw_tR_50_inv); cd(flat.amp_raw_tR_80_inv, bins[x].amp_raw_tR_80_inv);
            cd(flat.amp_raw_pos_slopes, bins[x].amp_raw_pos_slopes); cd(flat.amp_raw_systolic_peaks, bins[x].amp_raw_systolic_peaks);
            cd(flat.amp_raw_neg_slopes_pre_dnotch, bins[x].amp_raw_neg_slopes_pre_dnotch);
            cd(flat.amp_raw_dicrotic_notches, bins[x].amp_raw_dicrotic_notches);
            cd(flat.amp_raw_diastolic_peaks, bins[x].amp_raw_diastolic_peaks);
            cd(flat.amp_raw_neg_slopes_after_dnotch, bins[x].amp_raw_neg_slopes_after_dnotch);
            cd(flat.amp_baselined_feets, bins[x].amp_baselined_feets);
            cd(flat.amp_baselined_tP_20, bins[x].amp_baselined_tP_20); cd(flat.amp_baselined_tP_50, bins[x].amp_baselined_tP_50); cd(flat.amp_baselined_tP_80, bins[x].amp_baselined_tP_80);
            cd(flat.amp_baselined_tP_20_inv, bins[x].amp_baselined_tP_20_inv); cd(flat.amp_baselined_tP_50_inv, bins[x].amp_baselined_tP_50_inv); cd(flat.amp_baselined_tP_80_inv, bins[x].amp_baselined_tP_80_inv);
            cd(flat.amp_baselined_tR_20, bins[x].amp_baselined_tR_20); cd(flat.amp_baselined_tR_50, bins[x].amp_baselined_tR_50); cd(flat.amp_baselined_tR_80, bins[x].amp_baselined_tR_80);
            cd(flat.amp_baselined_tR_20_inv, bins[x].amp_baselined_tR_20_inv); cd(flat.amp_baselined_tR_50_inv, bins[x].amp_baselined_tR_50_inv); cd(flat.amp_baselined_tR_80_inv, bins[x].amp_baselined_tR_80_inv);
            cd(flat.amp_baselined_pos_slopes, bins[x].amp_baselined_pos_slopes); cd(flat.amp_baselined_systolic_peaks, bins[x].amp_baselined_systolic_peaks);
            cd(flat.amp_baselined_neg_slopes_pre_dnotch, bins[x].amp_baselined_neg_slopes_pre_dnotch);
            cd(flat.amp_baselined_dicrotic_notches, bins[x].amp_baselined_dicrotic_notches);
            cd(flat.amp_baselined_diastolic_peaks, bins[x].amp_baselined_diastolic_peaks);
            cd(flat.amp_baselined_neg_slopes_after_dnotch, bins[x].amp_baselined_neg_slopes_after_dnotch);
            cd(flat.proportional_pulse_amp, bins[x].proportional_pulse_amp);

            // SQI
            for (int i = 0; i < len && i < (int)bins[x].sqi.size(); i++) {
                auto& s = bins[x].sqi[i];
                if (s.size() >= 5) {
                    flat.sqi_mean_corr_dtw[prev + i] = s[0];
                    flat.sqi_corrcoff_direct[prev + i] = s[1];
                    flat.sqi_corrcoff_interp[prev + i] = s[2];
                    flat.sqi_dtw[prev + i] = s[3];
                    flat.sqi_frechet[prev + i] = s[4];
                }
            }

            prev += len;
        }

        // ── Inter-beat intervals (MATLAB: [0; diff(...)]) ───────────────────
        // MATLAB: begin_sec = flattened.idx_begin / ppgSampleRate
        // sec_valley_2_valley = [0; diff(begin_sec)]
        auto make_ibi_idx = [&](const std::vector<double>& idx_src) {
            std::vector<double> result(total_beats, NaN);
            result[0] = 0;
            for (int i = 1; i < total_beats; i++)
                result[i] = (idx_src[i] - idx_src[i - 1]) / ppgRate;
            return result;
            };

        // MATLAB: begin_sec + tP_20_x / ppgSampleRate → then diff
        auto make_ibi_composite = [&](const std::vector<double>& begin_sec,
            const std::vector<double>& offset_x) {
                std::vector<double> result(total_beats, NaN);
                result[0] = 0;
                for (int i = 1; i < total_beats; i++) {
                    double cur = begin_sec[i] + offset_x[i] / ppgRate;
                    double prv = begin_sec[i - 1] + offset_x[i - 1] / ppgRate;
                    result[i] = cur - prv;
                }
                return result;
            };

        // begin_sec for composite IBIs
        std::vector<double> begin_sec(total_beats);
        for (int i = 0; i < total_beats; i++)
            begin_sec[i] = flat.idx_begin[i] / ppgRate;

        flat.sec_valley_2_valley = make_ibi_idx(flat.idx_begin);
        flat.sec_foot_2_foot = make_ibi_idx(flat.idx_foot);

        flat.sec_tP_20_2_tP_20 = make_ibi_composite(begin_sec, flat.tP_20_x);
        flat.sec_tP_50_2_tP_50 = make_ibi_composite(begin_sec, flat.tP_50_x);
        flat.sec_tP_80_2_tP_80 = make_ibi_composite(begin_sec, flat.tP_80_x);
        flat.sec_tP_20_inv_2_tP_20_inv = make_ibi_composite(begin_sec, flat.tP_20_x_inv);
        flat.sec_tP_50_inv_2_tP_50_inv = make_ibi_composite(begin_sec, flat.tP_50_x_inv);
        flat.sec_tP_80_inv_2_tP_80_inv = make_ibi_composite(begin_sec, flat.tP_80_x_inv);

        flat.sec_pos_slope_2_pos_slope = make_ibi_idx(flat.idx_pos_slope);
        flat.sec_systolic_2_systolic = make_ibi_idx(flat.idx_systolic);
        flat.sec_neg_slope_b4_2_neg_slope_b4 = make_ibi_idx(flat.idx_neg_slope_b4);
        flat.sec_neg_slope_after_2_neg_slope_after = make_ibi_idx(flat.idx_neg_slope_after);
        flat.sec_diastolic_2_diastolic = make_ibi_idx(flat.idx_diastolic);
        flat.sec_dnotch_2_dnotch = make_ibi_idx(flat.idx_dnotch);

        flat.sec_tR_20_2_tR_20 = make_ibi_composite(begin_sec, flat.tR_20_x);
        flat.sec_tR_50_2_tR_50 = make_ibi_composite(begin_sec, flat.tR_50_x);
        flat.sec_tR_80_2_tR_80 = make_ibi_composite(begin_sec, flat.tR_80_x);
        flat.sec_tR_20_inv_2_tR_20_inv = make_ibi_composite(begin_sec, flat.tR_20_x_inv);
        flat.sec_tR_50_inv_2_tR_50_inv = make_ibi_composite(begin_sec, flat.tR_50_x_inv);
        flat.sec_tR_80_inv_2_tR_80_inv = make_ibi_composite(begin_sec, flat.tR_80_x_inv);

        // ── Sleep state adjustment (exact MATLAB) ───────────────────────────
        auto rl = RunLength(flat.sleep_stages);
        int before_sleep_wake_count = 0;
        for (size_t i = 0; i < rl.values.size(); i++) {
            if (rl.values[i] == 0) { before_sleep_wake_count = rl.lengths[i]; break; }
        }

        int after_sleep_count = 0;
        // MATLAB: wakeIDX = find(B(:,1) ~= 0 & ~isnan(B(:,1)), 1, 'last') + 1
        int wakeIDX = -1;
        for (int i = (int)rl.values.size() - 1; i >= 0; i--) {
            if (rl.values[i] != 0 && !std::isnan(rl.values[i])) {
                wakeIDX = i + 1;
                break;
            }
        }
        if (wakeIDX >= 0 && wakeIDX < (int)rl.values.size()) {
            after_sleep_count = rl.lengths[wakeIDX];
        }

        flat.adjusted_sleep_state = flat.sleep_stages;
        auto& adj = flat.adjusted_sleep_state;

        // 4 = awake after sleep
        if (after_sleep_count > 0) {
            for (int i = total_beats - after_sleep_count; i < total_beats; i++)
                if (i >= 0) adj[i] = 4;
        }
        // 2 = REM (was -1)
        for (auto& v : adj) if (v == -1) v = 2;
        // 1 = NREM (everything not 0, 2, or 4)
        for (auto& v : adj) if (v != 0 && v != 2 && v != 4 && !std::isnan(v)) v = 1;
        // 3 = awake during sleep (0 → 3)
        for (auto& v : adj) if (v == 0) v = 3;
        // 0 = awake before sleep
        for (int i = 0; i < before_sleep_wake_count && i < total_beats; i++)
            adj[i] = 0;

        // ── Correct times (exact MATLAB correctTimes) ───────────────────────
        flat.correct_idx_begin.assign(total_beats, 0);
        flat.edge_beat_mask.assign(total_beats, 0);
        flat.corrected_time_sec.assign(total_beats, 0);

        int begIdx = 0;
        for (int x = 0; x < (int)processSegments.size(); x++) {
            for (int y = 0; y < (int)processSegments[x].ppg_bin_indexs.size(); y++) {
                int bs = processSegments[x].ppg_bin_indexs[y][0];
                int be = processSegments[x].ppg_bin_indexs[y][1];
                int len_bin = be - bs;
                int endIdx = begIdx + len_bin;

                int first = -1, last = -1;
                for (int i = 0; i < total_beats; i++) {
                    double ib = flat.idx_begin[i];
                    if (ib >= begIdx && ib <= endIdx) {
                        flat.correct_idx_begin[i] = ib - begIdx + bs;
                        if (first < 0) first = i;
                        last = i;
                    }
                }
                if (first >= 0) flat.edge_beat_mask[first] = 1;
                if (last >= 0)  flat.edge_beat_mask[last] = 1;
                begIdx = endIdx;
            }
        }

        // MATLAB: remove first and last edge marks
        for (int i = 0; i < total_beats; i++)
            if (flat.edge_beat_mask[i] == 1) { flat.edge_beat_mask[i] = 0; break; }
        for (int i = total_beats - 1; i >= 0; i--)
            if (flat.edge_beat_mask[i] == 1) { flat.edge_beat_mask[i] = 0; break; }

        // MATLAB: corrected_time_sec = (correct_idx_begin - 1) / samplingRate
        for (int i = 0; i < total_beats; i++)
            flat.corrected_time_sec[i] = (flat.correct_idx_begin[i] - 1) / ppgRate;

        // ppg_flat_time_msec
        int tpl = 0;
        for (auto& s : processSegments) tpl += (int)s.po.size();
        flat.ppg_flat_time_msec.resize(tpl + 1);
        for (int i = 0; i <= tpl; i++)
            flat.ppg_flat_time_msec[i] = (double)i / ppgRate;

        // sec_to_first_onset_of_sleep and sec_from_last_onset_of_sleep
        double sot = (before_sleep_wake_count < total_beats) ?
            flat.corrected_time_sec[before_sleep_wake_count] : 0;
        int last_sleep_idx = total_beats - after_sleep_count - 1;
        double lst = (last_sleep_idx >= 0 && last_sleep_idx < total_beats) ?
            flat.corrected_time_sec[last_sleep_idx] : 0;

        flat.sec_to_first_onset_of_sleep.resize(total_beats);
        flat.sec_from_last_onset_of_sleep.resize(total_beats);
        for (int i = 0; i < total_beats; i++) {
            flat.sec_to_first_onset_of_sleep[i] = (flat.corrected_time_sec[i] - sot) * -1.0;
            flat.sec_from_last_onset_of_sleep[i] = flat.corrected_time_sec[i] - lst;
        }

        return flat;
    }

    // ─── Main Feature Generation (exact GenerateFeatures.m) ────────────────────
    inline int GenerateFeatures(const std::string& anneal_path,
        const std::string& wave_path,
        const std::string& template_info_path,
        const std::string& template_marking_path,
        const std::string& output_path) {
        namespace fs = std::filesystem;

        std::string stem = fs::path(anneal_path).stem().string();
        std::regex re("(\\d+_\\d+)");
        std::smatch match;
        std::string name = std::regex_search(stem, match, re) ? match[1].str() : stem;

        const int window_length = 30;
        const double sqi_threshold = std::numeric_limits<double>::infinity();

        std::cout << "Loading annealed segments..." << std::endl;
        std::vector<AnnealedSegment> annealedSegments;
        try {
            annealedSegments = read_annealed_bin(anneal_path);
            std::cout << "  -> " << annealedSegments.size() << " segments loaded" << std::endl;
        }
        catch (const std::exception& e) {
            std::cerr << "  FAILED reading annealed: " << e.what() << std::endl;
            return 0;
        }

        std::cout << "Loading wave data..." << std::endl;
        std::vector<WaveData> wave_data;
        try {
            wave_data = read_wave_data_bin(wave_path);
            std::cout << "  -> " << wave_data.size() << " bins loaded" << std::endl;
        }
        catch (const std::exception& e) {
            std::cerr << "  FAILED reading wave data: " << e.what() << std::endl;
            return 0;
        }

        bool has_template = !template_info_path.empty() && fs::exists(template_info_path);
        std::vector<TemplateInfo> template_info;

        if (has_template) {
            std::cout << "Loading template info..." << std::endl;
            try {
                template_info = read_template_info_bin(template_info_path);
                std::cout << "  -> " << template_info.size() << " templates loaded" << std::endl;
            }
            catch (const std::exception& e) {
                std::cerr << "  FAILED reading template info: " << e.what() << std::endl;
                return 0;
            }
            if (!template_marking_path.empty() && fs::exists(template_marking_path)) {
                std::cout << "Applying template markings..." << std::endl;
                apply_template_markings_bin(template_marking_path, template_info);
            }

            // Handle segment/template mismatch (exact MATLAB logic)
            if (annealedSegments.size() != template_info.size()) {
                int diff = (int)annealedSegments.size() - (int)template_info.size();
                bool fixed = false;
                if (diff == 1) {
                    TemplateInfo dummy;
                    dummy.TemplateBad = 1; dummy.bad_r_templates = 1; dummy.bad_ppg_templates = 1;
                    if (annealedSegments.front().po.size() <= 258) {
                        template_info.insert(template_info.begin(), dummy);
                        fixed = true;
                    }
                    else if (annealedSegments.back().po.size() <= 258) {
                        template_info.push_back(dummy);
                        fixed = true;
                    }
                }
                if (!fixed) {
                    std::cerr << "Error: segment/template mismatch for " << name << std::endl;
                    return 0;
                }
            }
        }

        std::cout << "Finding beat features..." << std::endl;
        BinMarksMask bmm = has_template ? countBins(wave_data, &template_info) : countBins(wave_data);
        int ns = (int)wave_data.size();
        std::vector<BeatsTotal> bib(ns);

        // Set thread count — use all available cores
        int nthreads = omp_get_max_threads();
        std::cout << "OpenMP threads: " << nthreads << std::endl;

        // Pre-allocate one SqiWorkspace per thread — persists across all segments.
        std::vector<std::unique_ptr<SqiWorkspace>> workspaces(nthreads);
        for (int i = 0; i < nthreads; i++)
            workspaces[i] = std::make_unique<SqiWorkspace>();

        // Segments processed in parallel — each thread handles one segment at a time,
        // processing beats serially within it. This maximizes core utilization
        // (no idle threads when a segment has few beats).
        // Each thread has its own pre-allocated SqiWorkspace.
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int t = 0; t < ns; t++) {
            try {
                int tid = omp_get_thread_num();
                auto& ws = *workspaces[tid];

#pragma omp critical
                std::cout << "Section " << (t + 1) << " of " << ns << std::endl;

                if (has_template) {
                    std::vector<int> ppgidxs;
                    for (int r = 0; r < wave_data[t].pairs.rows(); r++)
                        ppgidxs.push_back(wave_data[t].pairs.col(r, 0));
                    int wl = window_length * (int)annealedSegments[t].ppgSampleRate;

                    // Inline SQI computation using this thread's workspace
                    // (no nested parallelism — beats processed serially per segment)
                    auto sq = PPG_SQI_serial(annealedSegments[t].po, ppgidxs,
                        template_info[t].ppgTemplate, wl,
                        annealedSegments[t].ppgSampleRate, ws);

                    std::vector<double> sc(sq.typeMnemonic.size());
                    for (size_t j = 0; j < sq.typeMnemonic.size(); j++)
                        sc[j] = sq.typeMnemonic[j][0];

                    double sp2e = template_info[t].End - template_info[t].Peak;
                    double rsp = NaN;
                    if (!std::isnan(sp2e) && sp2e >= 1)
                        rsp = std::abs(template_info[t].Dicrotic - template_info[t].Peak) / sp2e;

                    auto fr = GetBeatFeaturesFromTemplate(
                        sc, sqi_threshold,
                        annealedSegments[t].po, annealedSegments[t].sleep_stages,
                        wave_data[t].pairs,
                        annealedSegments[t].ppgSampleRate, annealedSegments[t].ecgSampleRate,
                        rsp);
                    bib[t] = std::move(fr.beats_total);

                    for (auto& row : sq.typeMnemonic)
                        if (row[0] < sqi_threshold) bib[t].sqi.push_back(row);
                    bib[t].sqilabels = sq.labels;
                }
                else {
                    int npm = std::max(0, wave_data[t].pairs.rows() - 1);
                    std::vector<double> zs(npm, 0.0);
                    auto fr = GetBeatFeaturesFromTemplate(
                        zs, sqi_threshold,
                        annealedSegments[t].po, annealedSegments[t].sleep_stages,
                        wave_data[t].pairs,
                        annealedSegments[t].ppgSampleRate, annealedSegments[t].ecgSampleRate);
                    bib[t] = std::move(fr.beats_total);
                    bib[t].sqilabels = { "not_run" };
                }
                bib[t].error_ppg_segmentation = bmm.could_not_identify_PPG[t];
                bib[t].review_bad_ppg_template = bmm.ppg_template_manually_excluded[t];
                bib[t].review_bad_r_template = bmm.ecg_template_manually_excluded[t];
            }
            catch (const std::exception& e) {
#pragma omp critical
                std::cerr << "Section " << (t + 1) << " FAILED: " << e.what() << std::endl;
            }
            catch (...) {
#pragma omp critical
                std::cerr << "Section " << (t + 1) << " FAILED: unknown error" << std::endl;
            }
        }

        std::cout << "Flattening..." << std::endl;
        auto bf = flatten_beat_idx(bib, annealedSegments);

        // Build ppg_wout_noise (concatenate all segment PPGs with overlap-1)
        int tl = 0;
        for (auto& s : annealedSegments) tl += (int)s.po.size();
        bf.ppg_wout_noise.resize(tl);
        int pos = 0;
        for (auto& s : annealedSegments) {
            for (int k = 0; k < (int)s.po.size() && pos < tl; k++)
                bf.ppg_wout_noise[pos++] = s.po[k];
            if (pos > 0) pos--; // overlap by 1
        }

        int nb = 0;
        for (auto& b : bib) nb += (int)b.area.size();
        std::string out = (fs::path(output_path) / (name + "_feature_output.bin")).string();
        std::cout << "Saving to: " << out << std::endl;
        write_feature_output_bin(out, bf, nb);
        return 1;
    }

} // namespace ppg