#pragma once
//
// feature_marks.hpp
//
// Every marker/landmark computation lives in this single class. Three
// categories:
//
//   Fixed     -- auto-detected, NOT user-editable. Currently R peak.
//   Reactive  -- computed from movable markers' current positions;
//                drawn as X glyphs, no draggable bar. Q peak / S peak
//                inside the QRS, plus the PPG glyph bundle used by the
//                overlay.
//   Movable   -- draggable bars in the GUI; auto-detect provides the
//                initial seed, user can then drag.
//
// Also provides seed_all() -- one call that runs every auto-detector
// for one TemplateBin (across all channels + PPG + arterial), copies
// the results into the bin's *_auto_ch fields, and seeds any unset
// user fields with the fresh auto values.
//

#include <utility>
#include <vector>

struct TemplateBin;   // forward-declare -- full definition in TemplateBinIO.hpp

class FeatureMarks {
public:
    // =================================================================
    // Fixed
    // =================================================================

    // (r_peak()/detect_r_peak() removed: R is the deterministic anchor column
    //  r_col_raw carried on the template; QRS polarity is derived from the
    //  known R inside each detector.)

    // =================================================================
    // Reactive
    // =================================================================
    // Given the current positions of the movable markers, return the
    // computed reactive marker's sample index. -1 if inputs invalid.

    // Q peak = extreme (opposite R polarity) inside [q_begin, r_peak].
    static int compute_q_peak(const std::vector<double>& ecg,
        int q_begin, int r_peak);

    // S = first opposite-polarity trough after R (rate-aware window). The
    // single S-trough source, used for s_end detection and QRS-height
    // normalization; needs no s_end bound.
    static int compute_s_peak(const std::vector<double>& ecg, int r_idx, double fs);

    // Whole-beat reactive bundle for the GUI's glyph overlay. Each field is
    // auto-computed but tracks the user's movable markers live (see the
    // compute_* functions below).
    struct EcgGlyphs {
        int p_peak_glyph = -1;    // max within +/-0.05s of user P peak
        int q_begin_glyph = -1;   // cubic-fit knee within +/-0.05s of user Q begin
        int r_peak_glyph = -1;    // argmax|v-baseline| over [q_begin, s_end]
        int s_end_glyph = -1;     // recovery dropoff/knee within +/-0.05s of user S end
        int t_peak_glyph = -1;    // max value between user T begin and T end
        int t_end_glyph = -1;     // = the user's T-end marker (passthrough)
    };
    static EcgGlyphs compute_ecg_glyphs(const std::vector<double>& ecg,
        int p_peak, int q_begin, int s_end, int t_begin, int t_end, double fs);

    // Individual reactive computes (all track user markers live).
    static int compute_p_peak(const std::vector<double>& ecg, int p_onset, int q_begin, int r_peak_idx);
    static int compute_t_peak(const std::vector<double>& ecg, int tBegin, int tEnd);
    static int compute_r_peak(const std::vector<double>& ecg, int qBegin, int sEnd);
    static int compute_s_end(const std::vector<double>& ecg, int sUser, double fs);
    static int compute_q_onset(const std::vector<double>& ecg, int qUser, double fs, int r_idx);
    static int compute_t_end(const std::vector<double>& ecg, int tEndUser, double fs);
    static int compute_p_onset(const std::vector<double>& v, int pUser, double fs, int r_idx);
    static int compute_t_begin(const std::vector<double>& v, int tBeginUser, double fs, int r_idx);


    struct PpgGlyphs {
        int foot = -1;
        int p1 = -1, p1_fallback = -1;
        int p50 = -1, p50_fallback = -1;
        int dic = -1, dic_fallback = -1;
        int p2 = -1, p2_fallback = -1;
        int end = -1;
        bool notch_found = false;
    };
    static PpgGlyphs compute_ppg_glyphs(
        const std::vector<double>& ppg, int foot, int dic, int peak2);

    // =================================================================
    // Movable (auto-detected seeds)
    // =================================================================

    // ECG landmarks.
    static int detect_p_peak(const std::vector<double>& ecg_signal, int r_idx);
    static int detect_q_begin(const std::vector<double>& ecg_signal, int r_idx);
    static int detect_s_end(const std::vector<double>& ecg_signal, int r_idx, double fs);
    static int detect_t_begin(const std::vector<double>& ecg_signal, int r_idx, double fs);
    static int detect_t_end(const std::vector<double>& ecg_signal, int r_idx, double fs);

    // PPG landmarks.
    static int detect_ppg_onset(const std::vector<double>& pulse);
    static int detect_ppg_t80(const std::vector<double>& pulse);
    static int detect_ppg_peak(const std::vector<double>& pulse);
    static int detect_ppg_dicrotic(const std::vector<double>& pulse);
    static int detect_ppg_end(const std::vector<double>& pulse);

    static void seed_all(TemplateBin& bin, double sampleRate);

private:
    // Internal helpers.
    static std::vector<double> first_derivative(const std::vector<double>& v);
};