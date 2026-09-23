#pragma once
/*
* feature_marks.hpp
* This class represents the 3 types of markers on an ECG or PPG beat:
*   Fixed     -- auto-detected, NOT user-editable.
*   Reactive  -- computed from movable markers' current positions  drawn as X glyphs, no draggable bar.
*   Movable   -- draggable bars in the GUI; auto-detect provides the initial seed, user can then drag.
*
*  Also provides seed_all() -- one call that runs every auto-detector
*  for one TemplateBin (across all channels + PPG + arterial), copies
*  the results into the bin's *_auto_ch fields, and seeds any unset
*  user fields with the fresh auto values.
*/

#include <utility>
#include <vector>
#include <functional>
#include <cmath>
#include <cstdint>
#include <limits>

#include "template_marking_gui/anchor_view.hpp"   // AnchorType
#include "template_marking_gui/curve_fit.hpp"
#include "template_morphology_grouping/template_bank.hpp"
#include "subsample_refine.hpp"

struct TemplateBin;   // forward-declare -- full definition in TemplateBinIO.hpp

// Returns the landmark subsample index for one beat, or -1 if not found.
using AnchorLocator = std::function<double(const std::vector<double>& beat)>;

using AnchorLocatorD = AnchorLocator;

// Build the per-beat locator for one anchor. Binds r_col/fs into the detector.
AnchorLocator make_anchor_locator(AnchorType type, int r_col, double fs);

// Lead polarity, from the operator's per-channel "Lead Reversed" checkbox --
// the only polarity authority in the pipeline. There is deliberately no
// detect-from-signal fallback: the old qrs_positive_at decided from one sample
// against the whole-array median, so an R column slightly off the apex could
// flip a template and report the S trough's mirror as a Q peak.
//
// CHANNEL-INDEXED ON PURPOSE. The detectors take a bare sign; this exists for
// the callers that loop leads, where passing ch and the flag separately is the
// obvious way to get them out of step.
struct LeadPolarity {
    bool inverted[3] = { false, false, false };
    bool isInverted(int lead) const { return lead >= 0 && lead < 3 && inverted[lead]; }
    // What the finders take: multiply the trace by this and Q and S are troughs
    // below baseline, whichever way the lead was recorded.
    double sign(int lead) const { return isInverted(lead) ? -1.0 : 1.0; }
};

class FeatureMarks {
public:
    // Amplitude at a sub-sample position, linearly interpolated. Landmarks are
    // fractional doubles throughout, so "the trace's value at this landmark"
    // cannot index a rounded column. NaN if either bracketing sample is NaN --
    // a gap stays a gap.
    static double sample_at(const std::vector<double>& v, double p);

    struct ReactiveEcg { double t_peak = -1.0, p_peak = -1.0; };
    struct ReactivePpg { double t50 = -1.0, t80 = -1.0, t80_rise = -1.0, pw80 = -1.0, peak2 = -1.0; };
    static ReactiveEcg reactive_ecg(const std::vector<double>& ecg, double p_begin, double q_onset, double s_end, double t_end, double sampleRate,
        curve_fit::PeakFitMode peakMode = curve_fit::PeakFitMode::Auto);
    static ReactivePpg reactive_ppg(const std::vector<double>& ppg, double onset, double peak, double dicrotic, double end);

    // ---- ECG LANDMARK FINDERS ---------------------------------------------
    //
    // `sgn` IS NOT DEFAULTED, deliberately. It comes from
    // LeadPolarity::sign(lead), and a default would let a call site silently
    // claim "upright" -- which is worse than the median test it replaced,
    // because it is invisible. Leaving it required makes the compiler list
    // every caller that has not been converted.
    //
    // WHICH FINDERS TAKE IT: the ones whose landmark is defined by polarity.
    // Q is by definition the first NEGATIVE deflection before R, and S the
    // trough after it, so those searches need the lead the right way up.
    //
    // WHICH DO NOT, and this is not an omission: find_t_peak, find_t_end and
    // find_p_peak decide deflection direction LOCALLY, from a bracket baseline
    // and the sample furthest from it in either direction. An inverted T or P
    // on an upright lead is ordinary pathology, not lead reversal -- feeding
    // those the lead sign would make them walk past exactly the finding that
    // matters.
    static double find_q_peak(const std::vector<double>& ecg, int r_idx, double fs, double sgn, curve_fit::PeakFitMode peakMode = curve_fit::PeakFitMode::Auto);
    static double find_s_peak(const std::vector<double>& ecg, int r_idx, double fs, double sgn, curve_fit::PeakFitMode peakMode = curve_fit::PeakFitMode::Auto);
    static double find_j_point(const std::vector<double>& ecg, double fs, int r_col, double sgn, subsample_refine::TransitionCandidates* candOut = nullptr, curve_fit::FitMode mode = curve_fit::FitMode::Auto);
    static double find_q_onset(const std::vector<double>& ecg, double fs, int r_idx, double sgn, double qPeakIn = -1.0, bool* measured = nullptr, subsample_refine::TransitionCandidates* candOut = nullptr, curve_fit::FitMode mode = curve_fit::FitMode::Auto);
    static double find_p_begin(const std::vector<double>& v, double fs, int r_idx, double sgn, double pPeakIn = -1.0, subsample_refine::TransitionCandidates* candOut = nullptr, curve_fit::FitMode mode = curve_fit::FitMode::Auto);
    static int    find_p_end(const std::vector<double>& ecg_signal, int r_idx, double fs, double sgn, double pPeakIn = -1.0);

    // Locally-polarised finders -- no sgn, see the note above.
    static double find_t_peak(const std::vector<double>& ecg, double bracketSEnd, double bracketTEnd, curve_fit::PeakFitMode peakMode = curve_fit::PeakFitMode::Auto);
    static double find_t_end(const std::vector<double>& ecg, double fs, int r_col, double j_point = -1.0, subsample_refine::TransitionCandidates* candOut = nullptr, curve_fit::FitMode mode = curve_fit::FitMode::Auto);
    static double find_p_peak(const std::vector<double>& v, double loIn, double hiIn, double fs, curve_fit::PeakFitMode peakMode = curve_fit::PeakFitMode::Auto);


    struct TemplateLandmarks {
        double r_peak = -1.0;    // refined from nominal_r_col, sub-sample
        double q_onset = -1.0;  bool q_onset_found = false;
        double q_peak = -1.0;
        double s_end = -1.0;    // == J-point
        double t_end = -1.0;
        double p_peak = -1.0;
        double p_begin = -1.0;
        bool   valid = false;    // false => waveform or anchor unusable

        // The candidate curves the DETECTOR fitted, carried out so the focus
        // panel draws the same fits that placed the mark rather than fitting
        // its own. That is the whole invariant in focus_panel_widget.hpp: a
        // curve the panel produced itself would answer a different question
        // from the one that placed the mark.
        subsample_refine::TransitionCandidates q_onset_cand;
        subsample_refine::TransitionCandidates s_end_cand;
        subsample_refine::TransitionCandidates t_end_cand;
        subsample_refine::TransitionCandidates p_begin_cand;
    };

    static TemplateLandmarks detect_template_landmarks(const std::vector<double>& tmpl, int nominal_r_col, double sampleRate, double sgn,
        curve_fit::FitMode fitMode = curve_fit::FitMode::Auto,
        curve_fit::PeakFitMode peakMode = curve_fit::PeakFitMode::Auto);

    static void seed_bank_template(const std::vector<double>& tmpl, int r_col,
        double sampleRate, double sgn, AnchorType anchor, tbank::BankMarkerSet& out,
        curve_fit::FitMode fitMode = curve_fit::FitMode::Auto,
        curve_fit::PeakFitMode peakMode = curve_fit::PeakFitMode::Auto);

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

        subsample_refine::PeakCandidates peak_cand;    // systolic peak
        subsample_refine::PeakCandidates onset_cand;   // foot
        subsample_refine::PeakCandidates end_cand;     // end of cycle

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
    //
    // NO sgn: a pulse waveform has a physical sign, so there is no polarity
    // question to answer. Same for the arterial channels.
    static PpgFiducials detect_ppg_fiducials(const std::vector<double>& v, int W, double ppgRate, double heightMeters = NAN);
    static double amplitude_crossing(const std::vector<double>& v, int a, int b, double frac);
    static double crossing_at_level(const std::vector<double>& v, int a, int b, double target);
    static int trough_in(const std::vector<double>& v, int lo, int hi);
    static double steepest_slope_in(const std::vector<double>& v, int lo, int hi);
    static double first_crossing(const std::vector<double>& v, int a, int b, double frac);
    static int detect_ppg_upstroke_peak(const std::vector<double>& v, int lo = 0, int hi = -1);
    static int detect_ppg_onset(const std::vector<double>& pulse);
    static double detect_ppg_peak(const std::vector<double>& pulse);
    static int detect_ppg_dicrotic(const std::vector<double>& pulse, int peak);
    static double detect_ppg_peak2(const std::vector<double>& v, int sysPeak, double t80, int end);
    static int detect_ppg_end(const std::vector<double>& pulse);

    // pol is indexed by channel inside, so the sign and the lead cannot get out
    // of step across the three-channel loop.
    static void seed_all(TemplateBin& bin, double sampleRate, double ppgRate, AnchorType anchor,
        const LeadPolarity& pol, double heightMeters = NAN,
        curve_fit::FitMode fitMode = curve_fit::FitMode::Auto,
        curve_fit::PeakFitMode peakMode = curve_fit::PeakFitMode::Auto);

    static void seed_pulse_bank_template(const std::vector<double>& tmpl,
        double ppgRate, tbank::BankPulseMarkerSet& out, double heightMeters = NAN);
};