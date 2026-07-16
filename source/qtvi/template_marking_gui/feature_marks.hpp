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

    // R peak = argmax of |v - baseline| across the first 3/4 of `v`.
    // Baseline is the mean of that same window (DC-tolerant).
    // Returns { index, isPositive }. isPositive tells callers whether
    // the R spike goes UP or DOWN relative to baseline.
    static std::pair<int, bool> r_peak(const std::vector<double>& v);
    static int detect_r_peak(const std::vector<double>& ecg_signal);

    // =================================================================
    // Reactive
    // =================================================================
    // Given the current positions of the movable markers, return the
    // computed reactive marker's sample index. -1 if inputs invalid.

    // Q peak = extreme (opposite R polarity) inside [q_begin, r_peak].
    static int compute_q_peak(const std::vector<double>& ecg,
        int q_begin, int r_peak);

    // S peak = extreme (opposite R polarity) inside [r_peak, s_end].
    static int compute_s_peak(const std::vector<double>& ecg,
        int r_peak, int s_end);

    // Whole-beat reactive bundle for the GUI's glyph overlay. Each field is
    // auto-computed but tracks the user's movable markers live (see the
    // compute_* functions below).
    struct EcgGlyphs {
        int p_wave = -1;    // max within +/-0.05s of user P peak
        int q_onset = -1;   // cubic-fit knee within +/-0.05s of user Q begin
        int r_wave = -1;    // argmax|v-baseline| over [q_begin, s_end]
        int s_end = -1;     // recovery dropoff/knee within +/-0.05s of user S end
        int t_peak = -1;    // max value between user T begin and T end
        int t_end = -1;     // = the user's T-end marker (passthrough)
    };
    static EcgGlyphs compute_ecg_glyphs(const std::vector<double>& ecg,
        int p_peak, int q_begin, int s_end, int t_begin, int t_end, double fs);

    // Individual reactive computes (all track user markers live).
    static int compute_p_wave(const std::vector<double>& ecg, int pUser, double fs);
    static int compute_t_wave(const std::vector<double>& ecg, int tBegin, int tEnd);
    static int compute_r_wave(const std::vector<double>& ecg, int qBegin, int sEnd);
    static int compute_s_end(const std::vector<double>& ecg, int sUser, double fs);
    static int compute_q_onset(const std::vector<double>& ecg, int qUser, double fs);
    static int compute_t_end(const std::vector<double>& ecg, int tEndUser, double fs);

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
    static int detect_p_peak(const std::vector<double>& ecg_signal);
    static int detect_q_begin(const std::vector<double>& ecg_signal);
    static int detect_s_end(const std::vector<double>& ecg_signal);
    static int detect_t_begin(const std::vector<double>& ecg_signal);
    static int detect_t_end(const std::vector<double>& ecg_signal);

    // PPG landmarks.
    static int detect_ppg_onset(const std::vector<double>& pulse);
    static int detect_ppg_t80(const std::vector<double>& pulse);
    static int detect_ppg_peak(const std::vector<double>& pulse);
    static int detect_ppg_dicrotic(const std::vector<double>& pulse);
    static int detect_ppg_end(const std::vector<double>& pulse);

    // =================================================================
    // Bin-level seed
    // =================================================================
    // Run every auto-detector for one TemplateBin (all three ECG
    // channels + PPG + arterial). Populates every *_auto_ch field
    // fresh and seeds unset user-facing fields with the same values.
    // R peak's user field is always overwritten with the fresh auto
    // value (R is auto-only in the GUI).
    static void seed_all(TemplateBin& bin, double sampleRate);

private:
    // Internal helpers.
    static std::vector<double> first_derivative(const std::vector<double>& v);
    static int detect_s(const std::vector<double>& ecg_signal);   // for detect_s_end
};