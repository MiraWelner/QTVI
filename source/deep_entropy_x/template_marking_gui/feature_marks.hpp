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
#include "template_morphology_grouping/template_bank.hpp"
#include <functional>
#include <cmath>
#include <cstdint>
#include <limits>

struct TemplateBin;   // forward-declare -- full definition in TemplateBinIO.hpp

enum class AnchorType { P_ONSET, P_PEAK, Q_ONSET, R_PEAK, J_POINT, T_PEAK };//J_POINT = S_END

// Returns the landmark subsample index for one beat, or -1 if not found.
using AnchorLocator = std::function<double(const std::vector<double>& beat)>;

// Returns the landmark as a sub-sample (floating-point) position
using AnchorLocatorD = AnchorLocator;

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
    // Sub-sample positions. These three used to return int, rounding away a
    // refinement they had already computed; they now return the refined double
    // and nothing downstream rounds it. Callers needing the trace value at the
    // landmark use sample_at() rather than indexing a rounded column.
    // Linear interpolation of a trace at a fractional position. NaN outside
    // the array or across a gap. This is how an amplitude is read at a
    // sub-sample landmark: a mark at 104.37 has a value, and it is not
    // ecg[104].
    static double sample_at(const std::vector<double>& v, double p);

    static double compute_q_peak(const std::vector<double>& ecg, int r_idx, double fs);
    static double compute_s_peak(const std::vector<double>& ecg, int r_idx, double fs);
    struct ReactiveEcg { double t_peak = -1.0; };
    static ReactiveEcg reactive_ecg(const std::vector<double>& ecg, int t_begin, int t_end);
    struct ReactivePpg { double t50 = -1.0, t80 = -1.0, t80_rise = -1.0, pw80 = -1.0, peak2 = -1.0;  };
    static ReactivePpg reactive_ppg(const std::vector<double>& ppg, int onset, int peak, int dicrotic, int end);

    // Individual reactive computes (all track user markers live).
    static double compute_j_point(const std::vector<double>& ecg, double fs, int r_col); //ONE j-point calculation
    static double compute_q_onset(const std::vector<double>& ecg, double fs, int r_idx);
    static double compute_t_peak(const std::vector<double>& ecg, double tBegin, double tEnd);
    static double compute_t_end(const std::vector<double>& ecg, double fs, int r_col,  double tBeginIn = -1.0);
    static double compute_p_begin(const std::vector<double>& v, double fs, int r_idx, double pPeakIn = -1.0);

    // ---------------------------------------------------------------
    // Per-BANK-TEMPLATE landmark seeding
    // ---------------------------------------------------------------
    //
    // Fills one BankMarkerSet from a template's OWN waveform. Until this
    // existed, only slot 0 got bars: TemplateViewerWindow applied
    // TemplateBin::marks() -- which describes the bin's sinus template -- and
    // skipped every other column, because a PVC's Q-onset sits at a different
    // column than sinus's and drawing sinus's bars there would be simply
    // wrong, with a drag writing the wrong value back.
    //
    // The fix is not to copy the bin's marks but to RUN THE SAME DETECTORS on
    // the template's own median. Every landmark function below already takes
    // (waveform, fs, r_col) and holds no per-bin state, so a bank template is
    // just another waveform to them. That is why this is a dozen lines rather
    // than a second implementation.
    //
    // r_col is the template's own R column. A template whose r_col is < 0 gets
    // an all -1 marker set: no anchor means no landmark is locatable, and -1 is
    // the established "absent" value the drawing layer already understands.
    // ---------------------------------------------------------------------
    // THE canonical ECG landmark detector. One waveform in, six landmarks out.
    //
    // EVERY TEMPLATE IS ITS OWN TEMPLATE. There used to be three copies of this
    // logic -- seed_all's per-channel block, seed_bank_template, and
    // bin_archive's buildChannelArchive -- and they disagreed. seed_all refined
    // the R anchor with subsample_refine::symmetricExtremum before running the
    // finders; the other two passed the stored r_col through untouched. Since
    // the R anchor is the search origin for all six finders, a one-sample
    // difference there walked through every landmark, so the archive's numbers
    // did not match what the viewer displayed for the same bin. Worse, a bank
    // template with no r_col of its own borrowed the BIN's refined R -- a value
    // refined against a different waveform, which is most wrong exactly where
    // the morphologies differ most.
    //
    // The R anchor is refined HERE, against the waveform being measured. A
    // PVC's peak genuinely sits elsewhere than the sinus median's, so a bank
    // template's refined R may differ from its bin's. That is intended: the
    // bank exists because these are different morphologies.
    //
    // CONVENTION: -1 MEANS ABSENT, never clamped to an edge column. A
    // ventricular template has no P wave, and template_bank.hpp's contract on
    // markers_by_anchor requires every P-dependent feature to come out NaN
    // rather than 0. Callers must treat -1 as a valid state.
    //
    // Positions are sub-sample doubles. Callers that store integers round at
    // the point of storage rather than here, so the sub-sample values stay
    // available to whoever wants them.
    // ---------------------------------------------------------------------
    struct TemplateLandmarks {
        double r_peak = -1.0;    // refined from nominal_r_col, sub-sample
        double q_begin = -1.0;
        double s_end = -1.0;    // == J-point
        double t_end = -1.0;
        double p_peak = -1.0;
        double p_begin = -1.0;
        bool   valid = false;    // false => waveform or anchor unusable
    };

    static TemplateLandmarks detect_template_landmarks(
        const std::vector<double>& tmpl, int nominal_r_col, double sampleRate);

    static void seed_bank_template(const std::vector<double>& tmpl, int r_col,
        double sampleRate, tbank::BankMarkerSet& out);

    // the x, o, |, || or ||| markers for to mark the ppg and to be output in the csv
    struct PpgFiducials {
        double onset = -1.0;
        double peak = -1.0;
        double end = -1.0;
        double peak2 = -1.0;      bool peak2_found = false;
        double dicrotic = -1.0;   bool notch_found = false;
        double t80 = -1.0, t50 = -1.0;
        double t80_rise = -1.0;
        double t80_rise_y = std::numeric_limits<double>::quiet_NaN();
        double pw80 = -1.0;
        double u = -1.0, v = -1.0, w = -1.0;
        double a = -1.0, b = -1.0, c = -1.0, d = -1.0, e = -1.0, f = -1.0;
        double p1 = -1.0, p2 = -1.0;

        // Three-tier dicrotic-notch provenance (E-5). Which tier produced the
        // notch, and its normalized confidence. dn_tier: 1=IEM, 2=Windkessel,
        // 3=absent (matches ppg_dicrotic::DnResult::Tier).
        int    dn_tier = 3;         // ABSENT until a tier resolves it
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

    static double steepest_rise_in(const std::vector<double>& v, int lo, int hi);

    // FIRST crossing of the frac-of-amplitude level on [a, b], linearly
    // interpolated between the bracketing samples and rounded.
    static double first_crossing(const std::vector<double>& v, int a, int b, double frac);

    static bool qrs_positive_at(const std::vector<double>& ecg_signal, int r_idx);
    static double detect_p_peak(const std::vector<double>& ecg_signal, int r_idx, double fs);
    static int detect_p_end(const std::vector<double>& ecg_signal, int r_idx, double fs, double pPeakIn = -1.0);

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