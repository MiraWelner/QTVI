#pragma once
// ═══════════════════════════════════════════════════════════════════════════════
// ppg_features.hpp — Core data structures and forward declarations.
// Faithfully mirrors every MATLAB field. Uses flat arrays for speed.
// ═══════════════════════════════════════════════════════════════════════════════

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <numeric>
#include <limits>
#include <functional>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <optional>
#include <set>
#include <regex>
#include <cstring>
#include <memory>
#include <omp.h>

namespace ppg {

    // ─── Utility ────────────────────────────────────────────────────────────────
    constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

    inline double getVal(const double* arr, int len, int idx) {
        // MATLAB: try c=arr(x); catch c=nan; end
        // MATLAB 1-based idx → C++ 0-based: caller must pass 0-based idx
        if (idx < 0 || idx >= len) return NaN;
        return arr[idx];
    }
    inline double getVal(const std::vector<double>& arr, int idx) {
        return getVal(arr.data(), (int)arr.size(), idx);
    }

    // ─── Data Structures ────────────────────────────────────────────────────────

    struct Pairs {
        std::vector<std::vector<int>> data; // Nx2: [ppgIdx, ecgIdx]
        int rows() const { return (int)data.size(); }
        int col(int r, int c) const { return data[r][c]; }
    };

    struct AnnealedSegment {
        std::vector<double> po;              // PPG signal
        std::vector<double> sleep_stages;
        double ppgSampleRate = 256.0;
        double ecgSampleRate = 256.0;
        std::vector<std::vector<int>> ppg_bin_indexs; // Nx2
    };

    struct TemplateInfo {
        std::vector<double> ecgTemplate;
        std::vector<double> ppgTemplate;
        double alignment_point = NaN;
        int TemplateBad = 0;
        int bad_r_templates = 0;
        int bad_ppg_templates = 0;
        double Dicrotic = NaN;
        double Onset = NaN;
        double Peak = NaN;
        double End = NaN;
    };

    struct WaveData {
        Pairs pairs;
        int bad_segment = 0;
        std::vector<int> ch1_raw_idx, ch1_sq_idx, ch1_abs_idx;
        std::vector<int> ch2_raw_idx, ch2_sq_idx, ch2_abs_idx;
        std::vector<int> ch3_raw_idx, ch3_sq_idx, ch3_abs_idx;
        std::vector<double> ppgMaxAmps;
        std::vector<double> ppgMinAmps;
        uint8_t noise_flags[9] = {};
        std::vector<std::vector<int>> ppg_bin_indexs;
        std::vector<std::vector<int>> ecg_bin_indexs;
    };

    struct BinMarksMask {
        std::vector<bool> could_not_identify_PPG;
        std::vector<bool> poor_ppg_or_ecg_template_manually_excluded;
        std::vector<bool> ecg_template_manually_excluded;
        std::vector<bool> ppg_template_manually_excluded;
    };

    // ─── BeatsTotal: mirrors EVERY beats_total.* field from MATLAB ──────────────
    // All vectors are length=nBeats for this segment, initialized to NaN.
    struct BeatsTotal {
        // Scalars per-beat
        std::vector<double> sleep_stages, area_baselined, area;
        std::vector<double> amp_delta_systolic, abs_amp_foot, abs_amp_peak;

        // Indices (stored as double to match MATLAB, NaN for invalid)
        std::vector<double> idx_begin, idx_end, idx_foot;
        std::vector<double> tP_20_x, tP_50_x, tP_80_x;
        std::vector<double> tP_20_x_inv, tP_50_x_inv, tP_80_x_inv;
        std::vector<double> idx_pos_slope, idx_systolic;
        std::vector<double> idx_neg_slope_b4, idx_neg_slope_after;
        std::vector<double> idx_diastolic, idx_dnotch;
        std::vector<double> tR_20_x, tR_50_x, tR_80_x;
        std::vector<double> tR_20_x_inv, tR_50_x_inv, tR_80_x_inv;

        // Timings (msec)
        std::vector<double> msec_beat_length;
        std::vector<double> msec_tP_50_2_first_valley, msec_tP_50_2_foot;
        std::vector<double> msec_tP_50_2_tP_20, msec_tP_50_2_tP_80;
        std::vector<double> msec_tP_50_2_tP_20_inv, msec_tP_50_2_tP_80_inv;
        std::vector<double> msec_tP_50_2_pos_slope, msec_tP_50_2_systolic_peak;
        std::vector<double> msec_tP_50_2_negslopes_pre_dnotch;
        std::vector<double> msec_tP_50_2_dicrotic_notch;
        std::vector<double> msec_tP_50_2_diastolic_peak;
        std::vector<double> msec_tP_50_2_negslopes_post_dnotch;
        std::vector<double> msec_tP_50_2_second_valley;
        std::vector<double> msec_tP_50_2_tR_20, msec_tP_50_2_tR_50, msec_tP_50_2_tR_80;
        std::vector<double> msec_tP_50_2_tR_20_inv, msec_tP_50_2_tR_50_inv, msec_tP_50_2_tR_80_inv;
        std::vector<double> msec_total_duration_20, msec_total_duration_50, msec_total_duration_80;
        std::vector<double> msec_total_duration_tR_20, msec_total_duration_tR_50, msec_total_duration_tR_80;

        // R-peak timings
        std::vector<double> msec_R_2_first_valley, msec_R_2_foot;
        std::vector<double> msec_R_2_tP_20, msec_R_2_tP_50, msec_R_2_tP_80;
        std::vector<double> msec_R_2_tP_20_inv, msec_R_2_tP_50_inv, msec_R_2_tP_80_inv;
        std::vector<double> msec_R_2_pos_slope, msec_R_2_systolic_peak;
        std::vector<double> msec_R_2_negslopes_pre_dnotch, msec_R_2_dicrotic_notch;
        std::vector<double> msec_R_2_tR_20, msec_R_2_tR_50, msec_R_2_tR_80;
        std::vector<double> msec_R_2_tR_20_inv, msec_R_2_tR_50_inv, msec_R_2_tR_80_inv;
        std::vector<double> msec_R_2_diastolic_peak, msec_R_2_negslopes_post_dnotch;
        std::vector<double> msec_R_2_second_valley;

        // Raw amplitudes
        std::vector<double> amp_raw_vallies, amp_raw_feets;
        std::vector<double> amp_raw_tP_20, amp_raw_tP_50, amp_raw_tP_80;
        std::vector<double> amp_raw_tP_20_inv, amp_raw_tP_50_inv, amp_raw_tP_80_inv;
        std::vector<double> amp_raw_tR_20, amp_raw_tR_50, amp_raw_tR_80;
        std::vector<double> amp_raw_tR_20_inv, amp_raw_tR_50_inv, amp_raw_tR_80_inv;
        std::vector<double> amp_raw_pos_slopes, amp_raw_systolic_peaks;
        std::vector<double> amp_raw_neg_slopes_pre_dnotch, amp_raw_dicrotic_notches;
        std::vector<double> amp_raw_diastolic_peaks, amp_raw_neg_slopes_after_dnotch;

        // Baselined amplitudes
        std::vector<double> amp_baselined_feets;
        std::vector<double> amp_baselined_tP_20, amp_baselined_tP_50, amp_baselined_tP_80;
        std::vector<double> amp_baselined_tP_20_inv, amp_baselined_tP_50_inv, amp_baselined_tP_80_inv;
        std::vector<double> amp_baselined_tR_20, amp_baselined_tR_50, amp_baselined_tR_80;
        std::vector<double> amp_baselined_tR_20_inv, amp_baselined_tR_50_inv, amp_baselined_tR_80_inv;
        std::vector<double> amp_baselined_pos_slopes, amp_baselined_systolic_peaks;
        std::vector<double> amp_baselined_neg_slopes_pre_dnotch, amp_baselined_dicrotic_notches;
        std::vector<double> amp_baselined_diastolic_peaks, amp_baselined_neg_slopes_after_dnotch;

        std::vector<double> proportional_pulse_amp;

        // SQI (per-beat, 5 cols)
        std::vector<std::vector<double>> sqi;
        std::vector<std::string> sqilabels;

        // Flags (per-segment, not per-beat)
        bool error_ppg_segmentation = false;
        bool review_bad_ppg_template = false;
        bool review_bad_r_template = false;

        void resize(int n);
    };

    // ─── BeatsFlattened — flat struct, no std::map ──────────────────────────────
    struct BeatsFlattened {
        // All per-beat double vectors — same names as MATLAB flattened.*
        // Indices
        std::vector<double> idx_begin, idx_end, idx_foot;
        std::vector<double> tP_20_x, tP_50_x, tP_80_x;
        std::vector<double> tP_20_x_inv, tP_50_x_inv, tP_80_x_inv;
        std::vector<double> idx_pos_slope, idx_systolic;
        std::vector<double> idx_neg_slope_b4, idx_neg_slope_after;
        std::vector<double> idx_diastolic, idx_dnotch;
        std::vector<double> tR_20_x, tR_50_x, tR_80_x;
        std::vector<double> tR_20_x_inv, tR_50_x_inv, tR_80_x_inv;

        // Scalars
        std::vector<double> sleep_stages, area_baselined, area;
        std::vector<double> amp_delta_systolic, abs_amp_foot, abs_amp_peak;
        std::vector<double> msec_beat_length;

        // tP50-relative timings
        std::vector<double> msec_tP_50_2_first_valley, msec_tP_50_2_foot;
        std::vector<double> msec_tP_50_2_tP_20, msec_tP_50_2_tP_80;
        std::vector<double> msec_tP_50_2_tP_20_inv, msec_tP_50_2_tP_80_inv;
        std::vector<double> msec_tP_50_2_pos_slope, msec_tP_50_2_systolic_peak;
        std::vector<double> msec_tP_50_2_negslopes_pre_dnotch;
        std::vector<double> msec_tP_50_2_dicrotic_notch;
        std::vector<double> msec_tP_50_2_diastolic_peak;
        std::vector<double> msec_tP_50_2_negslopes_post_dnotch;
        std::vector<double> msec_tP_50_2_second_valley;
        std::vector<double> msec_tP_50_2_tR_20, msec_tP_50_2_tR_50, msec_tP_50_2_tR_80;
        std::vector<double> msec_tP_50_2_tR_20_inv, msec_tP_50_2_tR_50_inv, msec_tP_50_2_tR_80_inv;
        std::vector<double> msec_total_duration_20, msec_total_duration_50, msec_total_duration_80;
        std::vector<double> msec_total_duration_tR_20, msec_total_duration_tR_50, msec_total_duration_tR_80;

        // R-peak timings
        std::vector<double> msec_R_2_first_valley, msec_R_2_foot;
        std::vector<double> msec_R_2_tP_20, msec_R_2_tP_50, msec_R_2_tP_80;
        std::vector<double> msec_R_2_tP_20_inv, msec_R_2_tP_50_inv, msec_R_2_tP_80_inv;
        std::vector<double> msec_R_2_pos_slope, msec_R_2_systolic_peak;
        std::vector<double> msec_R_2_negslopes_pre_dnotch, msec_R_2_dicrotic_notch;
        std::vector<double> msec_R_2_tR_20, msec_R_2_tR_50, msec_R_2_tR_80;
        std::vector<double> msec_R_2_tR_20_inv, msec_R_2_tR_50_inv, msec_R_2_tR_80_inv;
        std::vector<double> msec_R_2_diastolic_peak, msec_R_2_negslopes_post_dnotch;
        std::vector<double> msec_R_2_second_valley;

        // Raw amplitudes
        std::vector<double> amp_raw_vallies, amp_raw_feets;
        std::vector<double> amp_raw_tP_20, amp_raw_tP_50, amp_raw_tP_80;
        std::vector<double> amp_raw_tP_20_inv, amp_raw_tP_50_inv, amp_raw_tP_80_inv;
        std::vector<double> amp_raw_tR_20, amp_raw_tR_50, amp_raw_tR_80;
        std::vector<double> amp_raw_tR_20_inv, amp_raw_tR_50_inv, amp_raw_tR_80_inv;
        std::vector<double> amp_raw_pos_slopes, amp_raw_systolic_peaks;
        std::vector<double> amp_raw_neg_slopes_pre_dnotch, amp_raw_dicrotic_notches;
        std::vector<double> amp_raw_diastolic_peaks, amp_raw_neg_slopes_after_dnotch;

        // Baselined amplitudes
        std::vector<double> amp_baselined_feets;
        std::vector<double> amp_baselined_tP_20, amp_baselined_tP_50, amp_baselined_tP_80;
        std::vector<double> amp_baselined_tP_20_inv, amp_baselined_tP_50_inv, amp_baselined_tP_80_inv;
        std::vector<double> amp_baselined_tR_20, amp_baselined_tR_50, amp_baselined_tR_80;
        std::vector<double> amp_baselined_tR_20_inv, amp_baselined_tR_50_inv, amp_baselined_tR_80_inv;
        std::vector<double> amp_baselined_pos_slopes, amp_baselined_systolic_peaks;
        std::vector<double> amp_baselined_neg_slopes_pre_dnotch, amp_baselined_dicrotic_notches;
        std::vector<double> amp_baselined_diastolic_peaks, amp_baselined_neg_slopes_after_dnotch;

        std::vector<double> proportional_pulse_amp;

        // SQI (5 named columns)
        std::vector<double> sqi_mean_corr_dtw, sqi_corrcoff_direct, sqi_corrcoff_interp;
        std::vector<double> sqi_dtw, sqi_frechet;

        // Inter-beat intervals
        std::vector<double> sec_valley_2_valley, sec_foot_2_foot;
        std::vector<double> sec_tP_20_2_tP_20, sec_tP_50_2_tP_50, sec_tP_80_2_tP_80;
        std::vector<double> sec_tP_20_inv_2_tP_20_inv, sec_tP_50_inv_2_tP_50_inv, sec_tP_80_inv_2_tP_80_inv;
        std::vector<double> sec_pos_slope_2_pos_slope, sec_systolic_2_systolic;
        std::vector<double> sec_neg_slope_b4_2_neg_slope_b4, sec_neg_slope_after_2_neg_slope_after;
        std::vector<double> sec_diastolic_2_diastolic, sec_dnotch_2_dnotch;
        std::vector<double> sec_tR_20_2_tR_20, sec_tR_50_2_tR_50, sec_tR_80_2_tR_80;
        std::vector<double> sec_tR_20_inv_2_tR_20_inv, sec_tR_50_inv_2_tR_50_inv, sec_tR_80_inv_2_tR_80_inv;

        // Sleep / time
        std::vector<double> adjusted_sleep_state;
        std::vector<double> corrected_time_sec;
        std::vector<double> sec_to_first_onset_of_sleep;
        std::vector<double> sec_from_last_onset_of_sleep;
        std::vector<double> ppg_flat_time_msec;
        std::vector<double> correct_idx_begin;
        std::vector<int>    edge_beat_mask;

        // Full PPG
        std::vector<double> ppg_wout_noise;
    };

    // ─── Pipeline Config ────────────────────────────────────────────────────────
    struct PipelineConfig {
        std::string dataType;
        std::string annealedPath;
        std::string wavePath;
        std::string templatePath;
        std::string markingPath;
        std::string featureOutputPath;
    };

    inline std::vector<PipelineConfig> readConfigCsv(const std::string& path) {
        std::vector<PipelineConfig> configs;
        std::ifstream f(path);
        if (!f.is_open()) return configs;
        std::string line;
        std::getline(f, line); // skip header
        while (std::getline(f, line)) {
            if (line.empty()) continue;
            std::stringstream ss(line);
            std::string col;
            std::vector<std::string> cols;
            while (std::getline(ss, col, ',')) cols.push_back(col);
            auto trim = [](std::string& c) {
                c.erase(0, c.find_first_not_of(" \t\r\n"));
                c.erase(c.find_last_not_of(" \t\r\n") + 1);
                };
            for (auto& c : cols) trim(c);
            if (cols.size() < 11) continue;
            PipelineConfig pc;
            pc.dataType = cols[0];
            pc.annealedPath = cols[6];
            pc.wavePath = cols[7];
            pc.templatePath = cols[8];
            pc.markingPath = cols[9];
            pc.featureOutputPath = cols[10];
            configs.push_back(pc);
        }
        return configs;
    }

    // ─── Forward Declarations ───────────────────────────────────────────────────

    struct RunLengthResult { std::vector<double> values; std::vector<int> lengths; };
    struct SqiResult { std::vector<std::vector<double>> typeMnemonic; std::vector<std::string> labels; };
    struct BeatFeaturesResult { BeatsTotal beats_total; };
    struct DtwResult { std::vector<double> w; int rows, cols; std::vector<int> ta, tb; };
    struct DpResult { std::vector<int> p, q; double cost; };

    // Pre-allocated workspace for SQI (defined in ppg_utils.hpp)
    struct SqiWorkspace;

    // ppg_utils.hpp
    inline std::vector<double> interp1_linear(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& xq);
    inline std::pair<std::vector<double>, std::vector<int>> PLA(const std::vector<double>& input, int s, double th);
    inline std::vector<double> nanfastsmooth(const std::vector<double>& data, int width);
    inline std::vector<double> PPGmedianfilter(const std::vector<double>& wave, int order, double Fs);
    inline RunLengthResult RunLength(const std::vector<double>& data);

    // find_foot_pulseox.hpp
    inline std::pair<double, int> find_foot_pulseox(const std::vector<double>& data);

    // dicrotic_sqi.hpp
    inline int dumbDicrotic(const std::vector<double>& beat, double sp_ratio = NaN);
    inline SqiResult PPG_SQI(const std::vector<double>& wave, const std::vector<int>& anntime,
        const std::vector<double>& tmplate, int windowlen, double Fs,
        std::vector<std::unique_ptr<SqiWorkspace>>& workspaces);
    inline SqiResult PPG_SQI_serial(const std::vector<double>& wave, const std::vector<int>& anntime,
        const std::vector<double>& tmplate, int windowlen, double Fs,
        SqiWorkspace& ws);

    // beat_features.hpp
    inline BeatFeaturesResult GetBeatFeaturesFromTemplate(
        const std::vector<double>& sqi, double threshold,
        const std::vector<double>& ppg, const std::vector<double>& sleepstates,
        const Pairs& pairs, double ppgSamplingRate, double ecgSamplingRate,
        double dnotch_ratio_sp = NaN);

    // flatten_and_generate.hpp
    inline BinMarksMask countBins(const std::vector<WaveData>& waveData,
        const std::vector<TemplateInfo>* tmplate = nullptr);
    inline BeatsFlattened flatten_beat_idx(const std::vector<BeatsTotal>& bins,
        const std::vector<AnnealedSegment>& processSegments);
    inline int GenerateFeatures(const std::string& anneal_path, const std::string& wave_path,
        const std::string& template_info_path,
        const std::string& template_marking_path,
        const std::string& output_path);

} // namespace ppg