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
#include "template_morphology_grouping/template_bank.hpp"
#include <functional>
#include <cmath>
#include <cstdint>
#include <limits>

struct TemplateBin;   // forward-declare -- full definition in TemplateBinIO.hpp

enum class AnchorType { P_ONSET, P_PEAK, Q_ONSET, R_PEAK, J_POINT, T_PEAK };//J_POINT = S_END

// Returns the landmark sample index for one beat, or -1 if not found.
// Anchor locators return SUB-SAMPLE positions. This was
// std::function<int(...)>, which forced every locator lambda to lround the
// landmark its finder had just refined -- the note that used to sit in
// make_anchor_locator ("Rounded only because AnchorLocator is declared
// int-returning; AnchorLocatorD is the sub-sample form and is still
// unimplemented") described exactly that. AnchorLocatorD was the intended
// replacement and is now what AnchorLocator is.
using AnchorLocator = std::function<double(const std::vector<double>& beat)>;

// Returns the landmark as a sub-sample (floating-point) position, per spec
// I-3. Built on top of make_anchor_locator's already-tested integer result:
// the integer locator supplies the seed (including its own fallback logic,
// reused as-is, not duplicated), and a subsample_refine method appropriate
// to that anchor type -- Gaussian-weighted quadratic for R (symmetric),
// cubic-with-analytic-derivative for P/T (asymmetric), or 4x-upsample-then-
// fit-and-select for the transition anchors (Q-onset, P-onset, J-point) --
// refines it to double precision.
// Retained as an alias so existing references keep compiling; AnchorLocator
// is now the same type.
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
    struct ReactivePpg { double t50 = -1.0, t80 = -1.0, t80_rise = -1.0, pw80 = -1.0; };   // sub-sample
    static ReactivePpg reactive_ppg(const std::vector<double>& ppg,
        int onset, int peak, int end);

    // Individual reactive computes (all track user markers live).
    static double compute_t_peak(const std::vector<double>& ecg, double tBegin, double tEnd);
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
    static void seed_bank_template(const std::vector<double>& tmpl, int r_col,
        double sampleRate, tbank::BankMarkerSet& out);

    // the x, o, |, || or ||| markers for to mark the ppg and to be output in the csv
    struct PpgFiducials {
        // SUB-SAMPLE POSITIONS. Every one of these is produced by a refinement
        // -- symmetricExtremum, asymmetricExtremum, a spline notch, a
        // fractional-crossing search -- so they are doubles and nothing inside
        // FeatureMarks rounds them. Rounding happens once, at the boundary
        // where a position enters TemplateBin's int fields (a versioned
        // on-disk format); see the note at that call site.
        double onset = -1.0;
        double peak = -1.0;
        double end = -1.0;
        double peak2 = -1.0;      bool peak2_found = false;
        double dicrotic = -1.0;   bool notch_found = false;
        double t80 = -1.0, p50 = -1.0;
        // T80_rise: upslope position at the SAME absolute amplitude as t80
        // (the 80%-downslope level), NOT an 80%-of-onset->peak crossing.
        // t80_rise_y is that amplitude; pw80 = t80 - t80_rise (the width).
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

    // ---------------------------------------------------------------
    // ORDERING INVARIANT: the notch precedes the diastolic peak
    // ---------------------------------------------------------------
    //
    // The dicrotic notch is aortic valve closure; the diastolic peak is the
    // reflected wave that arrives after it. notch < peak2, always, on every
    // pulse channel. Nothing enforced it, and the two are found by SEPARATE
    // detectors, so nothing made it true either:
    //
    //   PPG (detect_ppg_fiducials): peak2_found required splineDiastolic to be
    //   past the notch -- but when that test FAILED, peak2 fell back to the
    //   midpoint 0.5*(peak + hi) with no comparison against the notch at all.
    //   A late notch and an early midpoint put peak2 before the notch, marked
    //   not-found, and still written to ppg_peak2_auto and still drawn.
    //
    //   ABP / ART / ART_PULM (seed_pulse_channel): no test in either direction.
    //   The notch comes from detect_ppg_dicrotic on the cycle and peak2 from
    //   detect_ppg_peak2 on the same cycle, independently, and whichever the
    //   two searches land on is what gets stored.
    //
    // Enforced in ONE place so the channels cannot disagree about it. The notch
    // wins ties and is never moved: it has a three-tier detector with recorded
    // provenance and a confidence (E-5), while peak2 is a local-maximum search
    // with neither. Given a contradiction, the measurement that can say how it
    // was obtained is the one to keep.
    //
    // A violated pair is NOT silently swapped. Swapping would produce two
    // plausible landmarks in the right order and destroy the evidence that a
    // detector failed. peak2 is re-sought strictly after the notch, and if
    // nothing is there it is reported ABSENT -- which is a true statement, and
    // one the archive already knows how to carry (peak2_found).
    //
    // hi is the search bound (pulse end); returns false when the invariant
    // could not be satisfied and peak2 was cleared.
    static bool order_notch_before_peak2(const std::vector<double>& pulse,
        double notch, double& peak2, bool& peak2_found, double hi);

    // Integer overload for the arterial seeding path, which works in whole
    // samples. Same policy; -1 out means absent.
    static bool order_notch_before_peak2(const std::vector<double>& pulse,
        int notch, int& peak2, int hi);

    // Sample whose AMPLITUDE is frac of the way from v[a] to v[b] (NOT frac
    // of the sample-index distance). Shared by detect_ppg_fiducials (t80/
    // p50) and the GUI's reactive T80/P50 glyphs, so both always agree.
    static double amplitude_crossing(const std::vector<double>& v, int a, int b, double frac);

    // Position on [a, b] where v crosses the ABSOLUTE amplitude `target`
    // (interpolated between the straddling samples). Direction inferred from
    // v[a] vs v[b], same as amplitude_crossing. -1 if it never straddles.
    // Used for T80_rise: the upslope point at t80's own absolute level.
    static double crossing_at_level(const std::vector<double>& v, int a, int b, double target);

    // Index of the minimum sample on [lo, hi], NaN-skipping; -1 if none.
    static int trough_in(const std::vector<double>& v, int lo, int hi);

    // FIRST crossing of the frac-of-amplitude level on [a, b], linearly
    // interpolated between the bracketing samples and rounded.
    static double first_crossing(const std::vector<double>& v, int a, int b, double frac);

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

    static void seed_all(TemplateBin& bin, double sampleRate, double ppgRate, AnchorType anchor,
        double heightMeters = NAN);
private:
    // Internal helpers.
    static std::vector<double> first_derivative(const std::vector<double>& v);
};