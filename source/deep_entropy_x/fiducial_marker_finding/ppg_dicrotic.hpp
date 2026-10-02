#pragma once
//
// ppg_dicrotic.hpp
//
// DeepEntropyX Phase 2, E-5: dicrotic-notch (DN) detection.
//
// THE ONE DN DETECTOR. Every notch in the pipeline -- the viewer's notch bar,
// the markings files, and notch_ms in the transit logs -- comes from this
// module, so replacing a tier here replaces the notch everywhere.
//
// ===========================================================================
// WHAT IS IMPLEMENTED, AND WHAT IS NOT
// ===========================================================================
//
//   E-5.1  artifact gate              implemented, with ONE SUBSTITUTED
//                                     criterion -- see artifactGate()
//   E-5.3  envelope knots             implemented (extremaOfFirstDerivative)
//   E-5.4  IEM iteration              ppg_pipeline::iemEnvelope (NOT here --\n//                                     one IEM for the pipeline)
//   E-5.5  Tier 1 localization        implemented (tier1Dn)
//   E-5.6  Tier 2 flow reconstruction implemented (flowReconstructionDn),
//                                     with the four divergences from Hoeksel
//                                     named there -- NOT YET VALIDATED
//   E-5.7  windowed pass / join       implemented (detectDicroticNotchWindowed)
//                                     -- 4-second windows, half-window
//                                     advance, nearest-center wins,
//                                     disagreement archived as confidence,
//                                     boundary-only cycles discarded. This is
//                                     the correct entry point and the transit
//                                     log uses it.
//                                     detectDicroticNotchPulse remains for the
//                                     ONE caller that has no window to give:
//                                     see its note.
//   E-5.8  prominence descriptor      implemented, archived, never a gate
//   E-5.9  Tier 3                     implemented (UNDETECTED, index -1,
//                                     distinct from windowRejected)
//   E-5.10 configuration              implemented; the three deleted fields
//                                     (dn_min_prominence, iem_tol,
//                                     dn_enhance_gain) are ABSENT, not
//                                     neutralized -- along with the
//                                     multi-scale Gaussian pre-enhancement
//                                     they belonged to.
//
// ===========================================================================
// THE TWO SILENT FAILURES THIS MODULE IS BUILT TO AVOID
// ===========================================================================
//
// Both produce a zero residual; a zero residual produces no valley; and no
// valley is indistinguishable from "this pulse has no notch". Neither throws.
//
//   1. Locating envelope knots by the SIGN of the second derivative instead of
//      its SIGN CHANGE. Returns 0 knots on real data (verified: 0 and 0
//      against 41 and 40 on the same window), so no envelope can be built.
//      -> extremaOfFirstDerivative, and the acceptance test named on it.
//
//   2. Running the iteration on a single pulse rather than a multi-cycle
//      window: too few knots to constrain an envelope.
//      -> detectDicroticNotchPulse documents that this is what both current
//         callers do.
//
// Assert on knot counts, never on a DN position: a position test passes
// vacuously when the residual is zero.
//

#include <cmath>
#include <limits>
#include <vector>

namespace ppg_dicrotic {

    // =========================================================================
    // NO SHAPE OR PROMINENCE GATE ON TIER 1
    // =========================================================================
    //
    // An earlier revision of THIS file carried kDnMaxValleyWidthFrac, a
    // valley-width condition added on top of E-5.5's two. It is removed.
    //
    // E-5.5 is explicit: "There is no prominence gate... Prominence is still
    // computed, but as an archived per-pulse descriptor, never as an
    // acceptance criterion." E-5.1 is equally explicit for the gate's peak
    // extraction: the source uses scipy.signal.find_peaks with DEFAULT
    // settings -- a plain local-maximum test with no prominence, width or
    // distance constraint -- and "a C++ equivalent must match that behaviour,
    // not add constraints".
    //
    // It was also fixing a non-problem. The threshold was tuned on synthetic
    // pulses to reject a valley at 144 ms after the systolic peak in favour of
    // one at 200 ms. The reference results give systolic-peak-to-DN as
    // 159.5 ms with SD 35.8 ms, so 144 ms is WELL INSIDE one standard
    // deviation and was the better answer all along.

    // =========================================================================
    // E-5.10 Configuration
    // =========================================================================
    //
    // DELETED FROM THE PREVIOUS REVISION AND NOT PRESENT HERE, deliberately,
    // rather than set to a neutral value:
    //   dn_min_prominence   gating Tier 1 on prominence is what made Tier 1
    //                       return nothing. Prominence is still computed, as
    //                       an archived descriptor only (E-5.8).
    //   iem_tol             replaced by iem_beta, which is a DIFFERENT
    //                       criterion: absolute difference of successive
    //                       residual mean-square energies, not a normalized
    //                       norm ratio.
    //   dn_enhance_gain     together with the multi-scale Gaussian
    //                       pre-enhancement it drove -- a local extension
    //                       carried under the Pal citation, not part of the
    //                       published method.
    struct DnConfig {
        double iem_window_s = 4.0;      ///< multi-cycle window (E-5.7 pass)
        double iem_beta = 0.1;          ///< stop when |dE| < beta
        int    iem_max_iter = 12;       ///< SAFETY CEILING. Expect 1 to 3; a
        //                                   count at the ceiling is a problem.
        double ppg_prefilter_hz = 16.0; ///< 4th-order Butterworth, zero-phase
        int    sg_poly_order = 4;       ///< Savitzky-Golay order
        int    sg_nc = 0;               ///< SG window length, odd. 0 = scale
        //                                   with rate: 25 at 256 Hz, 49 at
        //                                   500 Hz (both ~0.098 s).
        double dn_window_lo_ms = 100.0; ///< earliest valley after the systolic
        //                                   peak (25 samples at 256 Hz)
        int    dn_gate_peak_min = 3;    ///< reject window AT OR BELOW this
        int    dn_gate_peak_max = 10;   ///< reject window ABOVE this
        double dn_window_overlap = 0.5; ///< window advance fraction (E-5.7)
        double tier2_rate_review_pct = 10.0;  ///< median per-record Tier 2
        //                                   rate triggering model-order review

        // Local additions, not from the source.
        /// The substituted first artifact criterion. See artifactGate().
        double dropout_run_ms = 100.0;

        /// E-5.2: nc is a WINDOW LENGTH, not a half-width, and must be odd and
        /// greater than sg_poly_order. 25 at 256 Hz (97.7 ms), 49 at 500 Hz
        /// (98.0 ms) -- about one half-width of the shortest cardiac feature,
        /// per Vannuccini. The authors name non-adaptive SG parameters as
        /// their first limitation and say explicitly that nc must be updated
        /// for a different sampling frequency; our rate is nearly double
        /// theirs, so this is not optional.
        int sgWindow(double fs) const;
        /// The same thing as a HALF-width, which is what ppg_deriv takes.
        /// Passing ppg_deriv's own default (12) would give nc = 25 at 500 Hz --
        /// half what E-5.2 requires.
        int sgHalfWidth(double fs) const { return (sgWindow(fs) - 1) / 2; }
    };

    // =========================================================================
    // Shared primitives
    // =========================================================================

    /// 4th-order Butterworth low-pass at `cutoffHz`, applied zero-phase
    /// (forward-backward). Built on the real 4th-order cascade, not a
    /// 2nd-order section run through filtfilt -- that trick doubles the order
    /// rather than preserving it.
    ///
    /// NaN-safe: NaNs are filled by linear interpolation between their finite
    /// neighbours (edge runs held at the nearest finite value) for the filter
    /// run, then restored, so the caller sees the same NaN positions it passed
    /// in. Returns the input unfiltered when there are fewer than 4 samples,
    /// the rate cannot support the cutoff, or every sample is NaN.
    std::vector<double> bandLimit(const std::vector<double>& pulse, double fs,
        double cutoffHz = 16.0);

    /// Min-max scale to [0, 1]: (x - min) / (max - min), over the finite
    /// samples. NaN in, NaN out at that position. Returns the input unchanged
    /// when there are fewer than 2 finite samples or the window is flat.
    std::vector<double> minMaxScale(const std::vector<double>& window);

    // NO DERIVATIVE BANK AND NO SPLINE HERE. Both already exist and there must
    // be one of each:
    //
    //   ppg_deriv::buildDerivatives     the Savitzky-Golay bank (d0 smoothed,
    //                                   d1/d2 in per-second units). E-5.4's
    //                                   pseudocode calls this.
    //   ppg_pipeline::naturalCubicSpline  E-5.3's `splineThrough`.
    //
    // This module duplicated both for a while. The copies differed -- a 49-
    // against 25-sample SG window at 500 Hz, and flat-held against
    // endpoint-anchored envelope tails -- which means they produced different
    // knots and a different residual from the pipeline's own. Deleted.

    // =========================================================================
    // E-5.3 Envelope knots
    // =========================================================================
    //
    // An envelope knot is a local extremum of the FIRST derivative, located by
    // a SIGN CHANGE of the SECOND derivative -- not by the sign of the second
    // derivative.
    //
    // At an extremum of d1, d2 crosses zero, so d2 is approximately zero
    // there. A test on the sign of d2 asks "is d2 positive here", which is
    // true almost everywhere and says nothing about where it changes. The
    // failure is total, not partial: on a real 4-second window the sign test
    // returns 0 maxima and 0 minima, against 41 and 40 for the sign-change
    // test.
    //
    // MANDATORY ACCEPTANCE TEST, not yet written: on record 3010104 at 256 Hz,
    // a 4-second window must yield at least 20 maxima and at least 20 minima
    // (measured median 38 of each, roughly nine per cardiac cycle, minimum 31
    // across 120 windows). Assert on THESE COUNTS, never on the resulting DN
    // position -- a position test passes vacuously when the residual is zero.

    /// `maxima` and `minima` are cleared and filled; the return value is both
    /// combined, in index order. `d1` is not read -- it is in the signature
    /// because it is the series these are extrema OF, and a caller validating
    /// that d1 and d2 describe the same window does so against it.
    std::vector<int> extremaOfFirstDerivative(const std::vector<double>& d1,
        const std::vector<double>& d2,
        std::vector<int>& maxima,
        std::vector<int>& minima);

    // =========================================================================
    // E-5.4 The IEM iteration -- ppg_pipeline::iemEnvelope
    // =========================================================================
    //
    // NOT IMPLEMENTED HERE. ppg_pipeline::iemEnvelope is the one IEM; see its
    // declaration for the three correctness points E-5.4 sets out and which of
    // them that function previously had wrong. Tier 1 reads its
    // `nonStationary` field -- the final residual, which is where the notch
    // lives.

    // =========================================================================
    // E-5.5 / E-5.6 / E-5.9 The result
    // =========================================================================
    struct DnResult {
        int    index = -1;                               ///< sample within the window
        double subSample = std::numeric_limits<double>::quiet_NaN();
        enum Tier { IEM = 1, FLOW = 2, UNDETECTED = 3 } tier = UNDETECTED;
        double prominence = std::numeric_limits<double>::quiet_NaN();  ///< archived descriptor, not a gate
        double confidence = std::numeric_limits<double>::quiet_NaN();  ///< between-window agreement (E-5.7)
        bool   windowRejected = false;                   ///< DISTINCT from UNDETECTED

        /// A real detection, from either tier. False for both UNDETECTED and
        /// windowRejected -- a caller needing to tell those apart reads the
        /// fields.
        bool found() const { return tier != UNDETECTED && index >= 0 && !windowRejected; }
        /// Sub-sample position where one was refined, else the integer index,
        /// else -1.
        double position() const {
            if (index < 0) return -1.0;
            return std::isfinite(subSample) ? subSample : static_cast<double>(index);
        }
    };

    // =========================================================================
    // E-5.1 Window preparation
    // =========================================================================
    //
    // FOUR STEPS, IN THIS ORDER. The order is not incidental -- steps 1 and 2
    // both change what step 3 sees:
    //
    //   1. BAND-LIMIT. 4th-order Butterworth low-pass at ppg_prefilter_hz
    //      (16 Hz), zero-phase. Removes high-frequency noise that would create
    //      spurious envelope knots -- and spurious GATE PEAKS: counted on the
    //      unfiltered window the peak count runs 6 to 17 on real MESA PPG and
    //      rejects 9 of 14 windows, against 5 to 11 and 1 of 14 after
    //      filtering.
    //   2. NORMALIZE. Min-max scale the window to [0, 1]. REQUIRED, not
    //      cosmetic: the IEM's stopping criterion is an ABSOLUTE energy
    //      difference (iem_beta), which is only meaningful on a fixed
    //      amplitude scale. Our PPG spans about 0.1 units, so without this the
    //      residual energies are ~100x smaller than the ones iem_beta was set
    //      against.
    //   3. ARTIFACT GATE. See artifactGate(); one criterion is substituted.
    //   4. RECORD THE OUTCOME. A rejected window yields no DN for any pulse it
    //      covers, and those pulses are marked WINDOW-REJECTED rather than
    //      notch-undetected. Different conditions, different flags.
    //
    // The ordering hazard in step 3 is observed, not predicted: applying the
    // gate to the signal as stored rejects 449 of 449 windows on record
    // 3010104, because most samples are negative.
    struct DnWindow {
        std::vector<double> raw;         ///< fs * iem_window_s samples of PPG
        std::vector<int>    pulseFeet;   ///< foot indices within the window
        std::vector<int>    sysPeaks;    ///< systolic peak indices within it
        double fs = 0.0;
        int centerSample = -1;
    };

    /// Steps 1 and 2 of E-5.1, in order. The gate and the IEM both run on the
    /// result, never on the raw window.
    std::vector<double> prepareWindow(const std::vector<double>& raw, double fs,
        const DnConfig& cfg = {});

    // =========================================================================
    // E-5.1 step 3: the artifact gate
    // =========================================================================
    //
    // Reject the window before any notch search runs on it. Three criteria,
    // any one of which rejects. The second and third are UNCHANGED from the
    // source; the first is NOT APPLIED AS PUBLISHED.
    //
    // ---- SUBSTITUTED CRITERION, recorded as a deviation ----
    //
    // Published criterion 1 is "the window contains zero or negative values".
    // It assumes a natively positive signal: the source's perioperative
    // monitors, where ABP in mmHg and PPG both sit on a positive scale, so a
    // non-positive sample indicates sensor failure.
    //
    // OUR PPG IS SIGNED AT SOURCE. Measured on MESA record 3010104, the pleth
    // channel has median -0.031 with 68.2% of samples negative. This is native
    // to the recording, not a product of our preprocessing: the coarse DC
    // removal in pipeline Section 2 applies to ECG channels only. Applying the
    // published criterion to our data rejects every window on every record --
    // confirmed at 449 of 449.
    //
    // The criterion's INTENT is to exclude windows containing invalid samples.
    // On a signed, zero-centred signal that condition appears as flat runs and
    // rail values rather than as negativity. Substituted with: reject if the
    // window contains a constant run longer than dropout_run_ms (default
    // 100 ms), if any sample sits at either ADC rail, or if the window is
    // constant throughout.
    //
    // ---- UNCHANGED CRITERIA ----
    //   * at or below dn_gate_peak_min (3) peaks exceeding the 75th percentile
    //     of the window amplitude
    //   * above dn_gate_peak_max (10) such peaks
    //
    // "75th percentile of the window amplitude" is read as 75% of the window's
    // own peak-to-trough range -- minMaxScale(window) > 0.75 -- not a
    // statistical quantile of the sample distribution. That makes the
    // threshold self-normalizing per window, which is what a per-window
    // amplitude test must be to mean the same thing at any absolute signal
    // level.

    /// ADC saturation bounds. Unknown (the default) SKIPS the rail check
    /// rather than guessing: a window's natural peak and a true saturated rail
    /// are not reliably distinguishable from samples alone. Rail values are
    /// meant to come from the Phase 1 signal-validity pass, which also supplies
    /// the clipping term of the PPG SQI in Task A -- the two share one source
    /// of truth. Until that pass supplies bounds, pass RailBounds{}.
    struct RailBounds {
        double lo = std::numeric_limits<double>::quiet_NaN();
        double hi = std::numeric_limits<double>::quiet_NaN();
        bool known() const { return std::isfinite(lo) || std::isfinite(hi); }
    };

    /// Which criterion fired -- `reject` alone does not say why, and the
    /// substituted criterion in particular has to be auditable per window.
    struct ArtifactGateResult {
        bool reject = false;
        bool constantRun = false;      ///< run of one repeated finite value,
        //                                  longer than dropout_run_ms
        bool atRail = false;           ///< sample at a supplied RailBounds
        bool constantWindow = false;   ///< every finite sample equal; checked
        //                                  independently of run length, so a
        //                                  NaN-interrupted flat window cannot
        //                                  dodge constantRun
        bool tooFewPeaks = false;      ///< <= dn_gate_peak_min (published)
        bool tooManyPeaks = false;     ///< >  dn_gate_peak_max (published)
        /// E-5.1's resolved ambiguity: BOTH bounds govern the ABOVE-PERCENTILE
        /// count, not all peaks. Both are logged so the choice stays auditable
        /// and reversible -- at 500 Hz the all-peaks reading rejects nearly
        /// every window.
        int  nPeaksAboveP75 = 0;
        int  nPeaksAll = 0;
        int  longestRunSamples = 0;
    };

    ArtifactGateResult artifactGate(const std::vector<double>& window, double fs,
        const DnConfig& cfg = {}, const RailBounds& rail = {});

    // =========================================================================
    // E-5.5 Tier 1 localization
    // =========================================================================
    //
    // Within the non-stationary component, searching between the systolic peak
    // and the DIASTOLIC ENDPOINT of that cardiac cycle, the DN is the FIRST
    // valley satisfying BOTH source conditions:
    //
    //   * it lies at least dn_window_lo_ms after the systolic peak
    //     (25 samples at 256 Hz, 50 at 500 Hz), and
    //   * its value in the non-stationary component is LESS THAN ZERO.
    //
    // The FIRST valley meeting both is taken -- not the deepest. Both
    // conditions exist because secondary reflected waves and non-physiological
    // oscillations each produce additional valleys in this interval.
    //
    // THERE IS NO THIRD CONDITION, and none may be added: E-5.5 states that
    // there is no prominence gate.
    //
    // The upper bound is the diastolic endpoint, the end of the cardiac cycle,
    // NOT a fraction of the RR interval. An earlier revision used 70% of RR.
    //
    // THERE IS NO PROMINENCE GATE. Prominence is computed and archived
    // (E-5.8), never used for acceptance.
    DnResult tier1Dn(const std::vector<double>& nonStationary, double fs,
        int sysPeak, int diastolicEnd, const DnConfig& cfg = {});

    // =========================================================================
    // E-5.6 Tier 2: flow reconstruction
    // =========================================================================
    //
    // Hoeksel et al. 1997, J Clin Monit 13(5):309-316. Arterial flow is
    // reconstructed from the pressure waveform and valve closure is taken at
    // the MINIMUM OF THE FIRST NEGATIVE DIP in reconstructed flow.
    //
    // DO NOT USE A ZERO-CROSSING. Under a two-element model with tau fitted to
    // the diastolic decay, the flow-proportional signal is identically zero
    // throughout ideal diastole by construction, so no isolated zero-crossing
    // exists and the landmark is set by noise. On synthetic pulses the
    // zero-crossing came out 230 to 260 ms late on clean data and returned no
    // answer at all on the no-notch case -- which is the case Tier 2 exists to
    // handle.
    //
    // FOUR DOCUMENTED DIVERGENCES FROM HOEKSEL, ALL REQUIRING VALIDATION
    // BEFORE TIER 2 OUTPUT IS TRUSTED:
    //   1. a two-element rather than three-element Windkessel model;
    //   2. tau fitted across the window rather than per epoch as published;
    //   3. PPG used as a pressure surrogate where Hoeksel validated on radial
    //      artery pressure;
    //   4. a post-systolic rather than full-epoch search.
    // Hoeksel found notch position less sensitive to parameter variation under
    // the nonlinear model; this is the simplest linear form.
    //
    // FIRING-RATE GATE, for whoever runs the cohort: report the Tier 2 rate as
    // a MEDIAN PER RECORD across the cohort, not pooled across pulses. Above
    // tier2_rate_review_pct (10%), revisit the model order. Any individual
    // record above 25% is flagged a data-quality outlier and investigated
    // separately, NOT used as evidence about model order. A median below 1%
    // means Tier 2 is effectively unreachable and triggers the dead-path
    // review in E-5.9, rather than silence.

    /// tau FITTED ACROSS THE WHOLE WINDOW, not per pulse. One pulse gives
    /// roughly 340 ms of diastole against a time constant near 550 ms -- less
    /// than one time constant. Hoeksel used 6-second epochs.
    ///
    /// Log-linear regression on the pooled diastolic segments seeds it, then
    /// two or three Gauss-Newton iterations refine it on the UNTRANSFORMED
    /// exponential: the log transform makes the noise multiplicative and
    /// over-weights the low-amplitude tail.
    ///
    /// `sysPeaks` and `cycleEnds` are parallel, one entry per cardiac cycle.
    /// Returns NaN when there is not enough diastole to fit.
    double fitTauWindow(const std::vector<double>& raw, double fs,
        const std::vector<int>& sysPeaks,
        const std::vector<int>& cycleEnds);

    /// g = dP/dt + P/tau is proportional to flow. Compliance scales g but
    /// cannot move its extremum, so only tau is needed.
    DnResult flowReconstructionDn(const std::vector<double>& raw, double fs,
        int sysPeak, double tauWindow, const DnConfig& cfg = {});

    // =========================================================================
    // E-5.8 Prominence, as an archived descriptor
    // =========================================================================
    //
    // The prominence of the detected valley in the non-stationary component,
    // archived per pulse alongside the position. NEVER an acceptance
    // criterion.
    //
    // This replaces the per-bin Tier 3 fraction as the arterial stiffening
    // marker. Pal reports 100% detectability, so under a correct Tier 1 the
    // Tier 3 fraction is near zero on every record and measures nothing but
    // defect rate. Prominence is continuous, is defined on every detected
    // pulse, and declines as the notch damps -- which is the physiological
    // quantity intended. ANY TIER 3 FRACTION ALREADY COMPUTED ON ANY DATASET
    // IS AN ARTIFACT OF THE BROKEN IMPLEMENTATION AND MUST BE DISCARDED, NOT
    // REINTERPRETED.
    //
    // IN RESIDUAL AMPLITUDE UNITS, not normalized. Whether to normalize, and
    // against what (the residual's range, the pulse's range), is a decision
    // the source does not make -- and normalizing against a quantity that
    // itself varies with notch damping would partly cancel the effect being
    // measured. Left raw and flagged here rather than chosen silently.
    double valleyProminence(const std::vector<double>& nonStationary, int valley);

    /// Width of a valley in SAMPLES, measured at half its prominence on each
    /// side. Exposed as a diagnostic only -- it is NOT an acceptance
    /// criterion, for the same reason prominence is not (E-5.5).
    /// 0 when `valley` is not a valley or its neighbourhood is not finite.
    int valleyWidthSamples(const std::vector<double>& nonStationary, int valley);

    // =========================================================================
    // E-5.7 The windowed pass and the join contract
    // =========================================================================
    //
    // E-5 RUNS ON WINDOWS. E-1 through E-4 run on pulses. The join is
    // specified, not implied:
    //
    //   * ADVANCE BY HALF A WINDOW (dn_window_overlap). Every interior cycle is
    //     covered centrally by exactly one window and peripherally by its
    //     neighbour.
    //   * NEAREST-CENTER WINS. When two windows both produce a DN for a pulse,
    //     keep the result from the window whose center is nearest that pulse.
    //   * ARCHIVE THE DISAGREEMENT. |t_A - t_B| between the two windows'
    //     estimates is stored as DnResult::confidence. Both are already
    //     computed, so it is free, and it directly measures the boundary
    //     effect: verification found the first cycle in each window differing
    //     from the rest by 2 to 5 samples, because spline envelopes are least
    //     constrained at window edges.
    //   * DISCARD BOUNDARY-ONLY CYCLES. A cycle covered by no window centrally
    //     -- at the start and end of a record -- yields no DN.
    //
    // confidence IS IN SAMPLES, not ms, to match the units the boundary effect
    // was verified in (2 to 5 samples). NaN where only one window covered the
    // cycle, which is not the same as zero: zero would read as two windows in
    // perfect agreement.
    //
    // "COVERED CENTRALLY" is the middle half of a window, [w0 + L/4,
    // w0 + 3L/4]. At a half-window advance that tiles the record exactly once,
    // which is what makes "exactly one window centrally" true. The source
    // states the property but not the interval; this is the interval that
    // yields it.
    //
    // A SIGNAL SHORTER THAN ABOUT TWO WINDOWS loses its edge cycles to the
    // boundary rule and may yield nothing at all. That is the specified
    // behaviour, not a failure, and nCyclesBoundaryOnly reports it rather than
    // leaving the caller to infer it from empty results.
    struct WindowedPass {
        /// Parallel to the sysPeaks/cycleEnds passed in, one entry per cycle.
        std::vector<DnResult> perCycle;
        /// 0 where the cycle was boundary-only, so its UNDETECTED is "never
        /// examined centrally" rather than "examined, no qualifying valley".
        /// The spec's DnResult has no field for that distinction and this
        /// keeps it without inventing one.
        std::vector<char> coveredCentrally;

        int nWindows = 0;
        int nWindowsRejected = 0;        ///< windows that failed the artifact gate
        int nTier1 = 0, nTier2 = 0;
        int nUndetected = 0;             ///< examined, no qualifying valley
        int nCyclesWindowRejected = 0;   ///< KEPT DISTINCT from nUndetected
        //                                   (E-5.9): their only central window
        //                                   failed the gate, so they were
        //                                   never examined at all
        int nBoundaryOnly = 0;           ///< no central coverage; discarded

        /// Tier 2 firing rate for THIS record, in percent of cycles that got a
        /// DN. E-5.6's gate is a MEDIAN PER RECORD across the cohort, not this
        /// number pooled across pulses -- so this is the per-record input to
        /// that median, not the gate itself.
        double tier2RatePct() const {
            const int det = nTier1 + nTier2;
            return (det > 0) ? 100.0 * nTier2 / det
                : std::numeric_limits<double>::quiet_NaN();
        }
    };

    /// THE CORRECT ENTRY POINT. `sysPeaks` and `cycleEnds` are parallel, one
    /// entry per cardiac cycle, as indices into `signal`, strictly increasing.
    WindowedPass detectDicroticNotchWindowed(const std::vector<double>& signal,
        double fs,
        const std::vector<int>& sysPeaks,
        const std::vector<int>& cycleEnds,
        const DnConfig& cfg = {},
        const RailBounds& rail = {});

    // =========================================================================
    // THERE IS NO PER-PULSE OR PER-TEMPLATE ENTRY POINT
    // =========================================================================
    //
    // detectDicroticNotchPulse IS DELETED. It took a single pulse -- in
    // practice a 1.9 RR template -- and ran the IEM on it. Nothing in E-5
    // describes that. E-5.1's DnWindow is "fs samples of PPG, iem_window_s
    // seconds": a contiguous slice of the RECORDING. sysPeaks and pulseFeet
    // are indices saying where the cycles fall WITHIN that slice; the pulses
    // are located in the window, not fed to it.
    //
    // Feeding it a template put the NEXT cardiac cycle inside the window --
    // a 1.9 RR frame contains most of one -- and that cycle's foot and
    // upstroke are the largest features present, so they dominated the spline
    // envelopes and distorted the residual inside the pulse of interest.
    // Measured: residual -0.13 to -0.20 at the next pulse's upstroke against
    // -0.002 for a graze in the first pulse's own diastole. Every notch defect
    // traced during implementation came from this, and each attempted repair
    // -- a valley-width gate, an undetected-on-missing-end rule, a
    // require-an-end rule -- was treating a symptom of the input rather than
    // the input. The width gate in particular added an acceptance criterion
    // E-5.5 forbids, to suppress a valley the second cycle had created.
    //
    // THE ONE ENTRY POINT IS detectDicroticNotchWindowed, over the continuous
    // signal, joined back per cycle as E-5.7 specifies. A notch is a property
    // of a cycle in the recording, measured once there. Anything that needs to
    // display one reads that measurement; it does not re-detect.

}  // namespace ppg_dicrotic