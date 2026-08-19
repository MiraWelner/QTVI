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
    static int compute_q_peak(const std::vector<double>& ecg, int r_idx, double fs);
    static int compute_s_peak(const std::vector<double>& ecg, int r_idx, double fs);
    struct ReactiveEcg { int t_peak = -1; };
    static ReactiveEcg reactive_ecg(const std::vector<double>& ecg, int t_begin, int t_end);
    struct ReactivePpg { int t50 = -1, t80 = -1; };
    static ReactivePpg reactive_ppg(const std::vector<double>& ppg,
        int onset, int peak, int end);

    // Individual reactive computes (all track user markers live).
    static int compute_t_peak(const std::vector<double>& ecg, int tBegin, int tEnd);
    static double compute_j_point(const std::vector<double>& ecg, double fs, int r_col); //ONE j-point calculation
    static double compute_q_onset(const std::vector<double>& ecg, double fs, int r_idx);
    // T-offset. Window [T-begin + 100 ms, T-begin + 200 ms] -- bounded by
    // T-begin, so no T-end seed is needed and the old circular bound (a T-end
    // estimate bounding the search for T-end) is gone. tBeginIn avoids
    // recomputing T-begin.
    static double compute_t_end(const std::vector<double>& ecg, double fs, int r_col,
        double tBeginIn = -1.0);
    // T-onset. Same onset algorithm as Q-onset / P-onset / the J-point: window
    // [J-point, J-point + 100 ms], baseline at the left edge (the recovered ST
    // level), anchor at 10% of the rise, 4x-upsampled fit-and-select. Has a
    // draggable bar and NO glyph. jPointIn avoids recomputing the J-point.
    static double compute_t_begin(const std::vector<double>& v, double fs, int r_idx,
        double jPointIn = -1.0);
    static double compute_p_begin(const std::vector<double>& v, double fs, int r_idx, double pPeakIn = -1.0);
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

    // Index of the minimum sample on [lo, hi], NaN-skipping; -1 if none.
    static int trough_in(const std::vector<double>& v, int lo, int hi);

    // FIRST crossing of the frac-of-amplitude level on [a, b], linearly
    // interpolated between the bracketing samples and rounded.
    static int first_crossing(const std::vector<double>& v, int a, int b, double frac);

    // Augmentation tally: how many scored pulses had their tallest sample away
    // from the systolic peak, i.e. how many the old argmax detector mislocated.
    // Read after a run; reset per record if you want per-record figures.
    static long long ppg_pulses_scored();
    static long long ppg_tallest_not_first();
    static void       ppg_reset_tally();

    //movable ecg bars
    static bool qrs_positive_at(const std::vector<double>& ecg_signal, int r_idx);
    static double detect_p_peak(const std::vector<double>& ecg_signal, int r_idx, double fs);
    static int detect_p_end(const std::vector<double>& ecg_signal, int r_idx, double fs, double pPeakIn = -1.0);
    static int detect_q_begin(const std::vector<double>& ecg_signal, int r_idx);

    // PPG landmarks
    //
    // Systolic peak located from the UPSTROKE, on [lo, hi) (hi <= 0 => whole
    // vector). Returns -1 when no upstroke is present (flat, all-NaN, too
    // short) -- callers must handle that rather than substituting index 0.

    static int detect_ppg_upstroke_peak(const std::vector<double>& v,
        int lo = 0, int hi = -1);
    static int detect_ppg_onset(const std::vector<double>& pulse);
    static double detect_ppg_peak(const std::vector<double>& pulse);
    static int detect_ppg_dicrotic(const std::vector<double>& pulse, int peak);
    static int detect_ppg_peak2(const std::vector<double>& pulse);
    static int detect_ppg_end(const std::vector<double>& pulse);

    static void seed_all(TemplateBin& bin, double sampleRate, double ppgRate, AnchorType anchor);
private:
    // Internal helpers.
    static std::vector<double> first_derivative(const std::vector<double>& v);
};