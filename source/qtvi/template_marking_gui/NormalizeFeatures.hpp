#pragma once
/*
The ECG Normalization algorithm is as follows:

    1. Calculate the total QRS vector magnitude for each beat to control for cardiac axis rotation:
        RS_peak(t) = abs(R_peak(t)) + abs(S_peak(t))
    2. Find the global reference by taking the median across all bins for each individual:
        Global_Ref_person = median(RS_peak(t))
    3. Normalize any amplitude feature like P, R, or T waves using the equation:
        Feature_peak_norm_abs = Feature_peak(t) / Global_Ref_person


The PPG Normalization algorithm is as follows:
    1. First, calculate the local PI for each beat:
        PI(t) = ((systolic_peak(t) - diastolic_trough(t)) / abs(diastolic_trough(t))) * 100.
    2. Find the global reference by taking the median PI across the entire recording for that individual:
        Global_Ref_person = median(PI(t)).
    3. Normalize your amplitude feature by converting it to its local baseline ratio first and then
       dividing by the global reference:
        Feature_peak_norm_abs = Feature_Local_Ratio(t) / Global_Ref_person.

*/

#include "TemplateBinIO.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace normalize_features {

    inline double median_finite(std::vector<double> v) {
		// Remove NaN and Inf values, then compute the median of the remaining finite values.
        v.erase(std::remove_if(v.begin(), v.end(),
            [](double x) { return !std::isfinite(x); }), v.end());
        if (v.empty()) return std::nan("");
        const size_t mid = v.size() / 2;
        std::nth_element(v.begin(), v.begin() + mid, v.end());
        const double a = v[mid];
        if (v.size() % 2 == 1) return a;
        // Even count: pair with max of the lower half.
        double b = -std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < mid; ++i) if (v[i] > b) b = v[i];
        return 0.5 * (a + b);
    }

    inline double sample_y(const std::vector<double>& v, int idx) {
        // Look up sample y-value at marker index, NaN if unavailable.
        if (idx < 0 || idx >= static_cast<int>(v.size())) return std::nan("");
        const double y = v[idx];
        return std::isnan(y) ? std::nan("") : y;
    }

    inline double ecg_norm(double raw, double ref) {
        //ref = median over bins of |R| + |S|
        if (!std::isfinite(ref) || ref == 0.0 || std::isnan(raw)) return raw;
        return raw / ref;
    }

    inline double calculate_perfusion_index(double y, double foot_y) {
        //PI = 100 * (y - foot) / |foot|
        if (std::isnan(y) || std::isnan(foot_y) || std::abs(foot_y) < 1e-12)
            return std::nan("");
        return 100.0 * (y - foot_y) / std::abs(foot_y);
    }

    inline double pulse_norm(double y, double foot_y, double ref) {
		//divide PPG by the global reference (PI)
        if (!std::isfinite(ref) || ref == 0.0) return std::nan("");
        const double lr = calculate_perfusion_index(y, foot_y);
        return std::isnan(lr) ? lr : lr / ref;
    }

    // ------------------------------------------------------------------
    // ECG per-channel Global_Ref = median over bins of |R_peak| + |S_peak|.
    // ------------------------------------------------------------------
    inline double compute_ecg_global_ref(const std::vector<TemplateBin>& bins, int ch, double sampleRateHz)
    {
        std::vector<double> vals;
        vals.reserve(bins.size());
        for (const auto& b : bins) {
            if (b.bad_segment) continue;
            if (b.bad_r_ch[ch]) continue;
            const ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };
            const auto& ecg = chs[ch]->ecgTemplate_raw;
            if (ecg.empty()) continue;

            EcgFeatures f = computeEcgFeatures(ecg,
                b.p_peak_ch[ch], b.q_begin_ch[ch], b.r_peak_ch[ch],
                b.s_end_ch[ch], b.t_end_ch[ch], sampleRateHz);
            const double ry = sample_y(ecg, f.r_idx);
            const double sy = sample_y(ecg, f.s_idx);
            if (std::isnan(ry) || std::isnan(sy)) continue;
            vals.push_back(std::abs(ry) + std::abs(sy));
        }
        return median_finite(std::move(vals));
    }

    // ------------------------------------------------------------------
    // Pulse channel accessor (0=PPG, 1=ABP, 2=ART, 3=ART_PULM).
    // ------------------------------------------------------------------
    struct PulseChannel {
        const std::vector<double>* trace;
        int   foot_idx;
        int   peak_idx;
        uint8_t issue;   // 0 = ok, 1 = user-bad, 2 = absent
    };

    inline PulseChannel pulseChan(const TemplateBin& b, int which) {
        switch (which) {
        case 0: return { &b.ppgTemplate,     b.ppg_onset,    b.ppg_peak,    b.ppg_issue };
        case 1: return { &b.abpTemplate,     b.abp_onset,    b.abp_peak,    b.abp_issue };
        case 2: return { &b.artTemplate,     b.art_onset,    b.art_peak,    b.art_issue };
        default: return { &b.artPulmTemplate, b.art_pulm_onset, b.art_pulm_peak, b.art_pulm_issue };
        }
    }

    inline double compute_pulse_global_ref(const std::vector<TemplateBin>& bins, int which)
    {
        //PI(bin) = 100 * (peak_y - foot_y) / |foot_y|
        std::vector<double> vals;
        vals.reserve(bins.size());
        for (const auto& b : bins) {
            if (b.bad_segment) continue;
            const PulseChannel pc = pulseChan(b, which);
            if (pc.issue != 0) continue;
            if (pc.trace->empty()) continue;
            const double foot_y = sample_y(*pc.trace, pc.foot_idx);
            const double peak_y = sample_y(*pc.trace, pc.peak_idx);
            vals.push_back(calculate_perfusion_index(peak_y, foot_y));
        }
        return median_finite(std::move(vals));
    }

}   // namespace normalize_features