#pragma once
/*
* feature_marks.hpp
* This class represents the 3 types of markers on an ECG or PPG beat:
*   Fixed     -- auto-detected, NOT user-editable.
*   Reactive  -- computed from movable markers' current positions  drawn as X glyphs, no draggable bar.
*    Movable   -- draggable bars in the GUI; auto-detect provides the initial seed, user can then drag.
*
*  Also provides seed_all() -- one call that runs every auto-detector
* for one TemplateBin (across all channels + PPG + arterial), copies
* the results into the bin's *_auto_ch fields, and seeds any unset
* user fields with the fresh auto values.
*/

#include <utility>
#include <vector>
#include "template_anchoring\anchor_type.hpp"   // AnchorType (re-exported for existing users)
#include "template_morphology_grouping/template_bank.hpp"
#include <functional>
#include <cmath>
#include <cstdint>
#include <limits>

struct TemplateBin;   // forward-declare -- full definition in TemplateBinIO.hpp

// Returns the landmark subsample index for one beat, or -1 if not found.
using AnchorLocator = std::function<double(const std::vector<double>& beat)>;

using AnchorLocatorD = AnchorLocator;

// Build the per-beat locator for one anchor. Binds r_col/fs into the detector.
AnchorLocator make_anchor_locator(AnchorType type, int r_col, double fs);
class FeatureMarks {
public:
    static double sample_at(const std::vector<double>& v, double p);// Returns the landmark as a sub-sample (floating-point) position
    struct ReactiveEcg { double t_peak = -1.0, p_peak = -1.0; };
    static ReactiveEcg reactive_ecg(const std::vector<double>& ecg, int p_begin, int q_begin, int s_end, int t_end, double sampleRate);
    struct ReactivePpg { double t50 = -1.0, t80 = -1.0, t80_rise = -1.0, pw80 = -1.0, peak2 = -1.0; };
    static ReactivePpg reactive_ppg(const std::vector<double>& ppg, int onset, int peak, int dicrotic, int end);

    static double compute_q_peak(const std::vector<double>& ecg, int r_idx, double fs);
    static double compute_t_end(const std::vector<double>& ecg, double fs, int r_col, double j_point = -1.0);
    static double compute_p_begin(const std::vector<double>& v, double fs, int r_idx, double pPeakIn = -1.0);
    static double compute_t_peak(const std::vector<double>& ecg, double bracketSEnd, double bracketTEnd);
    static double compute_s_peak(const std::vector<double>& ecg, int r_idx, double fs);
    static double compute_j_point(const std::vector<double>& ecg, double fs, int r_col);
    static double compute_q_onset(const std::vector<double>& ecg, double fs, int r_idx, double qPeakIn = -1.0, bool* measured = nullptr);
    static double compute_p_peak(const std::vector<double>& ecg, double pBegin, double qBegin, double fs);
    static double seed_p_peak(const std::vector<double>& ecg_signal, int r_idx, double fs);
    static int detect_p_end(const std::vector<double>& ecg_signal, int r_idx, double fs, double pPeakIn = -1.0);


    struct TemplateLandmarks {
        double r_peak = -1.0;    // refined from nominal_r_col, sub-sample
        double q_begin = -1.0;  bool q_begin_found = false;
        double q_peak = -1.0;
        double s_end = -1.0;    // == J-point
        double t_end = -1.0;
        double p_peak = -1.0;
        double p_begin = -1.0;
        bool   valid = false;    // false => waveform or anchor unusable
    };

    static TemplateLandmarks detect_template_landmarks(const std::vector<double>& tmpl, int nominal_r_col, double sampleRate);

    static void seed_bank_template(const std::vector<double>& tmpl, int r_col,
        double sampleRate, AnchorType anchor, tbank::BankMarkerSet& out);

    // the x, o, |, || or ||| markers for to mark the ppg and to be output in the csv
    struct PpgFiducials {
        double onset = -1.0;
        double peak = -1.0;
        double end = -1.0;
        double peak2 = -1.0;
        double dicrotic = -1.0;   bool notch_found = false;
        double t80 = -1.0, t50 = -1.0;
        double t80_rise = -1.0;
        double t80_rise_y = std::numeric_limits<double>::quiet_NaN();
        double pw80 = -1.0;
        double u = -1.0, v = -1.0, w = -1.0;
        double a = -1.0, b = -1.0, c = -1.0, d = -1.0, e = -1.0, f = -1.0;
        double p1 = -1.0, p2 = -1.0;


        int    dn_tier = 3;         //we have not implemented this yet
        double dn_confidence = 0.0;

        // Derived indices (DeepEntropyX Section 6.3). NaN when the points they
        // depend on are absent; SI stays NaN unless subject height is supplied.
        double ba = NAN, ca = NAN, da = NAN, ea = NAN, fa = NAN;
        double agi = NAN;   // (b - c - d - e)/a on the APG (aging index)
        double ri = NAN;   // amp(p2)/amp(p1) on the pulse (reflection index)
        double si = NAN;   // height / (t_p2 - t_p1) (stiffness index)
        uint16_t foundMask = 0;   // bit k (u,v,w,a,b,c,d,e,f,p1,p2) set when >= 0
    };

    //find all fiducial markers for PPG

    // heightMeters is the subject height from the demographics record, used
    // only for the stiffness index SI. Absent (default NaN) => SI is left NaN.
    static PpgFiducials detect_ppg_fiducials(const std::vector<double>& v, int W, double ppgRate,
        double heightMeters = NAN);

    // Sample whose AMPLITUDE is frac of the way from v[a] to v[b] (NOT frac
    // of the sample-index distance). Shared by detect_ppg_fiducials (t80/
    // t50) and the GUI's reactive T80/t50 glyphs, so both always agree.
    static double amplitude_crossing(const std::vector<double>& v, int a, int b, double frac);

    // Position on [a, b] where v crosses the ABSOLUTE amplitude `target`
    // (interpolated between the straddling samples). Direction inferred from
    // v[a] vs v[b], same as amplitude_crossing. -1 if it never straddles.
    // Used for T80_rise: the upslope point at t80's own absolute level.
    static double crossing_at_level(const std::vector<double>& v, int a, int b, double target);

    // Index of the minimum sample on [lo, hi], NaN-skipping; -1 if none.
    static int trough_in(const std::vector<double>& v, int lo, int hi);

    static double steepest_slope_in(const std::vector<double>& v, int lo, int hi);

    // FIRST crossing of the frac-of-amplitude level on [a, b], linearly
    // interpolated between the bracketing samples and rounded.
    static double first_crossing(const std::vector<double>& v, int a, int b, double frac);

    static bool qrs_positive_at(const std::vector<double>& ecg_signal, int r_idx);
   
    static int detect_ppg_upstroke_peak(const std::vector<double>& v, int lo = 0, int hi = -1);
    static int detect_ppg_onset(const std::vector<double>& pulse);
    static double detect_ppg_peak(const std::vector<double>& pulse);
    static int detect_ppg_dicrotic(const std::vector<double>& pulse, int peak);
    static int detect_ppg_peak2(const std::vector<double>& pulse);
    static int detect_ppg_end(const std::vector<double>& pulse);

    static void seed_all(TemplateBin& bin, double sampleRate, double ppgRate, AnchorType anchor,
        double heightMeters = NAN);

    static void seed_pulse_bank_template(const std::vector<double>& tmpl,
        double ppgRate, tbank::BankPulseMarkerSet& out, double heightMeters = NAN);
};