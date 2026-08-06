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
#include <functional>

struct TemplateBin;   // forward-declare -- full definition in TemplateBinIO.hpp

enum class AnchorType { P_ONSET, P_PEAK, Q_ONSET, R_PEAK, J_POINT, T_PEAK };//J_POINT = S_END

// Returns the landmark sample index for one beat, or -1 if not found.
using AnchorLocator = std::function<int(const std::vector<double>& beat)>;

// Returns the landmark as a sub-sample (floating-point) position, per spec
// I-3. Built on top of make_anchor_locator's already-tested integer result:
// the integer locator supplies the seed (including its own fallback logic,
// reused as-is, not duplicated), and a subsample_refine method appropriate
// to that anchor type -- Gaussian-weighted quadratic for R (symmetric),
// cubic-with-analytic-derivative for P/T (asymmetric), or 4x-upsample-then-
// fit-and-select for the transition anchors (Q-onset, P-onset, J-point) --
// refines it to double precision.
using AnchorLocatorD = std::function<double(const std::vector<double>& beat)>;

// Build the per-beat locator for one anchor. Binds r_col/fs into the detector.
AnchorLocator make_anchor_locator(AnchorType type, int r_col, double fs);
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

    // ---- Reactive glyph bundles ------------------------------------------
    // The ONLY definition of the bracketed glyphs. Each is a pure function of
    // a trace plus the bracketing bar positions, so the result is NEVER
    // stored: the GUI calls it per repaint, the writers call it per row, and
    // no cache exists that could go stale. Pass the *_auto bars to get an
    // autodetect column, the user bars to get a user column -- that choice is
    // the caller's, but the formula is only ever here.
    struct ReactiveEcg { int t_peak = -1; };
    static ReactiveEcg reactive_ecg(const std::vector<double>& ecg,
        int t_begin, int t_end);

    struct ReactivePpg { int t50 = -1, t80 = -1; };
    static ReactivePpg reactive_ppg(const std::vector<double>& ppg,
        int onset, int peak, int end);

    // Individual reactive computes (all track user markers live).
    static int compute_t_peak(const std::vector<double>& ecg, int tBegin, int tEnd);
    // The ONE J-point calculation. Transition detector only (40-sample window,
    // 4x cubic upsample, fit-and-select via transitionAnchor); re-derives the
    // S peak internally, so it needs no seed. Returns a sub-sample double.
    static double compute_j_point(const std::vector<double>& ecg, double fs, int r_col);
    static double compute_q_onset(const std::vector<double>& ecg, int qUser, double fs, int r_idx);
    static double compute_t_end(const std::vector<double>& ecg, int tEndUser, double fs, int r_col);
    static int compute_t_begin(const std::vector<double>& v, int tBeginUser, double fs, int r_idx);
    static double compute_p_begin(const std::vector<double>& v, int pUser, double fs, int r_idx);

    // Single source of truth for EVERY PPG fiducial. Used by seed_all() to
    // populate the *_auto fields AND by the GUI to draw the frozen glyphs --
    // there is no separate/independent glyph recompute anymore, so "peak"
    // (etc.) can't mean two different things in two different places.
    // Dependency order: end is computed FIRST and independently; peak/onset
    // are then bounded by end; the dicrotic notch is bounded by peak and
    // end; t80/p50 are amplitude crossings derived from peak/onset/end.
    // end/peak/notch never exceed W-1, so nothing can land in the invisible
    // tail-overlap region past the displayed window.
    struct PpgFiducials {
        int onset = -1;
        int peak = -1;
        int peak2 = -1;      bool peak2_found = false;
        int end = -1;        bool end_found = false;
        int dicrotic = -1;   bool notch_found = false;
        int t80 = -1;
        int p50 = -1;
    };

    //find all fiducial markers for PPG
    static PpgFiducials detect_ppg_fiducials(const std::vector<double>& v, double tR1, int W, double ppgRate);

    // Sample whose AMPLITUDE is frac of the way from v[a] to v[b] (NOT frac
    // of the sample-index distance). Shared by detect_ppg_fiducials (t80/
    // p50) and the GUI's reactive T80/P50 glyphs, so both always agree.
    static int amplitude_crossing(const std::vector<double>& v, int a, int b, double frac);

    // =================================================================
    // Movable (auto-detected seeds)
    // =================================================================

    // ECG landmarks.
    // QRS polarity from the known R column: positive if the R sample sits
    // above the trace baseline (median), negative otherwise. Exposed
    // publicly so callers outside feature_marks.cpp (e.g. sqi_ecg.hpp's
    // next-beat P-wave search) can search in the same polarity every
    // detector in this file already uses internally.
    static bool qrs_positive_at(const std::vector<double>& ecg_signal, int r_idx);
    static double detect_p_peak(const std::vector<double>& ecg_signal, int r_idx, double fs);
    static int detect_p_end(const std::vector<double>& ecg_signal, int r_idx, double fs);
    static int detect_q_begin(const std::vector<double>& ecg_signal, int r_idx);
    static int detect_t_begin(const std::vector<double>& ecg_signal, int r_idx, double fs);
    static int detect_t_end(const std::vector<double>& ecg_signal, int r_idx, double fs);

    // PPG landmarks.
    static int detect_ppg_onset(const std::vector<double>& pulse);
    static double detect_ppg_peak(const std::vector<double>& pulse);
    static int detect_ppg_dicrotic(const std::vector<double>& pulse);
    static int detect_ppg_peak2(const std::vector<double>& pulse);
    static int detect_ppg_end(const std::vector<double>& pulse);

    static void seed_all(TemplateBin& bin, double sampleRate, double ppgRate, AnchorType anchor);
private:
    // Internal helpers.
    static std::vector<double> first_derivative(const std::vector<double>& v);
};