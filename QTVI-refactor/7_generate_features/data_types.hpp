#pragma once

#include "beat_features.hpp"
#include <vector>
#include <string>
#include <cstdint>

/**
 * @file data_types.hpp
 * @brief Top-level data structures for the PPG feature-extraction pipeline.
 *
 * These correspond to the MATLAB structs loaded from .mat files
 * (annealedSegments, wave_data, template_info).
 */

namespace ppg {

    // ─── Input data ─────────────────────────────────────────────────────────────

    /**
     * @brief One annealed segment (contiguous, noise-removed PPG + metadata).
     */
    struct AnnealedSegment {
        std::vector<double> po;              ///< PPG signal.
        std::vector<double> sleep_stages;    ///< Per-sample sleep stage annotation.
        double ppg_sample_rate = 2000.0;     ///< PPG sampling rate (Hz). Always 2000 by this stage.
        double ecg_sample_rate = 2000.0;     ///< ECG sampling rate (Hz). Always 2000 by this stage.
        std::vector<std::pair<int, int>> ppg_bin_indices; ///< Original bin index ranges.
    };

    /**
     * @brief Wave-detection results for one segment.
     */
    struct WaveData {
        std::vector<std::pair<int, int>> pairs; ///< (ppg_valley_idx, ecg_R_idx). Second = -1 if unmatched.
        bool bad_segment = false;               ///< Flag: segment could not be processed.
    };

    /**
     * @brief Template information for one segment.
     */
    struct TemplateInfo {
        std::vector<double> ppg_template;  ///< Beat template waveform.
        std::vector<double> ecg_template;  ///< ECG beat template (may be empty).
        double alignment_point = kNaN;
        bool template_bad = false;
        bool bad_r_templates = false;
        bool bad_ppg_templates = false;
        double dicrotic = kNaN;
        double onset = kNaN;
        double peak = kNaN;
        double end = kNaN;
    };

    /**
     * @brief Per-segment bin quality flags.
     */
    struct BinMarks {
        bool could_not_identify_ppg = false;
        bool poor_template_excluded = false;
        bool ecg_template_excluded = false;
        bool ppg_template_excluded = false;
    };

    // ─── Output data ────────────────────────────────────────────────────────────

    /**
     * @brief Aggregated output: all beats flattened across segments.
     */
    struct FlattenedBeats {
        std::vector<BeatFeatures> beats;           ///< Per-beat features.
        std::vector<double>       ppg_wout_noise;  ///< Concatenated clean PPG.

        // Derived inter-beat intervals (sec)
        std::vector<double> sec_valley_to_valley;
        std::vector<double> sec_foot_to_foot;

        // Sleep state (adjusted: 0=awake-pre, 1=NREM, 2=REM, 3=awake-during, 4=awake-post)
        std::vector<double> adjusted_sleep_state;
        std::vector<double> sec_to_first_onset;
        std::vector<double> sec_from_last_onset;
        std::vector<double> corrected_time_sec;
        std::vector<int>    edge_beat_mask;

        double ppg_sample_rate = 2000.0;
    };

} // namespace ppg