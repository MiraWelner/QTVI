#pragma once
//
// ppg_pipeline.hpp
//
// Continuous-recording PPG signal conditioning. Three stages:
//
//   iemEnvelope    Iterative Envelope Mean (E-5.4). THE one IEM: ppg_dicrotic's
//                  Tier 1 calls this, it is not reimplemented there.
//   naturalCubicSpline
//                  the cubic spline E-5.3 cites as `splineThrough`. Exposed
//                  because the IEM's envelopes and E-1's DC baseline both need
//                  it and there must be exactly one.
//   dcEnvelope     E-1 (Section 6.1). Cubic spline through the diastolic troughs
//                  of a continuous trace, then a 2nd-order Butterworth low-pass
//                  at 0.1 Hz applied zero-phase (forward-backward).
//   perfusionIndex E-2 (Section 6.2). (sys - dc)/dc * 100, with a 1st-percentile
//                  denominator floor and sub-0.1% beats excluded.
//
// The IEM operates on ONE template pulse; the DC envelope and perfusion index
// operate on a continuous recording (many beats) -- a 0.1 Hz low-pass needs
// several seconds of signal and is meaningless on a single-pulse template.
//

#include <vector>

namespace ppg_pipeline {

    // ---- Natural cubic spline (E-5.3's `splineThrough`) -----------------
    //
    // Through knots (xs, ys), xs strictly increasing, sampled at every integer
    // in [0, N). Tails outside the knot span are left NaN; callers needing full
    // coverage anchor the endpoints themselves (see dcEnvelope, and
    // envelopeThrough inside the IEM).
    std::vector<double> naturalCubicSpline(const std::vector<int>& xs,
        const std::vector<double>& ys, int N);

    // ---- Iterative Envelope Mean (E-5.4) --------------------------------
    //
    // ONE IEM FOR THE WHOLE PIPELINE. ppg_dicrotic's Tier 1 calls this.
    //
    // THREE CORRECTNESS POINTS, each of which independently broke an earlier
    // revision of this function. All three are now as E-5.4 specifies:
    //
    //   * KNOTS ARE EXTREMA OF THE FIRST DERIVATIVE, located by SIGN CHANGES
    //     OF THE SECOND (ppg_dicrotic::extremaOfFirstDerivative). This
    //     previously used a plain 3-point extremum test on the signal itself,
    //     which is a different and smaller knot set.
    //   * THE NON-STATIONARY COMPONENT IS THE FINAL RESIDUAL, not the
    //     converged mean. `nonStationary` is where the dicrotic notch lives.
    //     The old field names had this inverted: `imf` was the residual and
    //     `envelope` the trend, which reads as the opposite of E-5.4.
    //   * THE STOP CRITERION is the ABSOLUTE difference of successive residual
    //     mean-square energies, |E{R_{i-1}^2} - E{R_i^2}| < beta. It was
    //     Huang's normalised SD ratio, which E-5.4 names as the wrong form.
    //
    // Envelopes are built on the SG-SMOOTHED signal; subtraction is from the
    // UNSMOOTHED current iterate. "Unfiltered" in the source distinguishes the
    // iterate from its smoothed version -- it does NOT mean the original input,
    // and subtracting from the original every round does not converge.
    //
    // Convergence takes 1 to 3 iterations. maxIter = 12 is a safety ceiling; a
    // count at the ceiling indicates a problem.
    struct IemEnvelope {
        std::vector<double> nonStationary;   // final residual: the DN is here
        std::vector<double> stationary;      // accumulated envelope means
        std::vector<double> upper;           // last-iteration upper envelope
        std::vector<double> lower;           // last-iteration lower envelope
        int  iterations = 0;
        bool converged = false;
        bool ok = false;
    };

    /// fs is needed for the Savitzky-Golay derivative bank the knots come from.
    /// beta is E-5.4's stopping threshold on the residual mean-square energy.
    IemEnvelope iemEnvelope(const std::vector<double>& pulse, double fs,
        int maxIter = 12, double beta = 0.1);

    // ---- E-1: DC envelope (continuous trace) ----------------------------
    // troughs are the diastolic-trough sample indices; fs is the PPG rate (Hz).
    // Returns a full-length baseline (endpoints held flat beyond the outermost
    // troughs so the envelope is defined across the whole record).
    std::vector<double> dcEnvelope(const std::vector<double>& ppg,
        const std::vector<int>& troughs, double fs);

    // ---- E-2: Perfusion index (%) with denominator floor ----------------
    // sys and dc are parallel arrays (per-beat or per-sample). The denominator
    // is floored at the 1st percentile of the finite dc values so a near-zero
    // baseline cannot blow the ratio up; beats with PI < 0.1% are excluded
    // (returned as NaN, preserving index alignment with the inputs).
    std::vector<double> perfusionIndex(const std::vector<double>& sys,
        const std::vector<double>& dc);

} // namespace ppg_pipeline