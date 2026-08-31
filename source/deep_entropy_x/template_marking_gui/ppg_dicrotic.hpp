#pragma once
//
// ppg_dicrotic.hpp
//
// E-5: three-tier dicrotic-notch (DN) detection (DeepEntropyX Phase 2, Section 6.3).
//
// The DN marks aortic valve closure. It is prominent in young compliant vessels
// and progressively damped with age/stiffness/vasodilation/low perfusion, to the
// point of vanishing. Simple derivative-zero methods fail on damped pulses, so
// the pipeline runs three tiers in order and records which one produced the answer:
//
//   Tier 1  Iterative Envelope Mean (IEM)  -- primary; detrend the pulse so the
//           notch is a clear local minimum of the residual even when it is only
//           an inflection in the raw pulse. (Pal et al. 2024.)
//   Tier 2  Windkessel-augmented fallback  -- physics-aware; reconstruct a
//           flow-proportional signal and take its zero-crossing. (Saffarpour 2023.)
//   Tier 3  ABSENT                         -- notch genuinely absent; recorded as
//           a feature, never substituted with a nominal position.
//
// Search window (all tiers): 120 ms after the systolic peak to 70% of the current
// RR after the systolic peak.
//
// NOTE ON PHASE-1 DEPENDENCIES. This module references two Phase 1 routines that
// are not in this repo snapshot: splineThrough (Task N, cubic spline) and
// refineAsymmetric (Task I-3, sub-sample refinement). Local stand-ins are provided
// in the .cpp under those exact names; replace them with the Phase 1 versions when
// integrating.
//

#include <cmath>
#include <limits>
#include <vector>

namespace ppg_dicrotic {

    // DN-related config fields (Section 6.3). May be merged into a broader PpgConfig.
    struct PpgConfig {
        double dnEnhanceGain = 0.15;   // dn_enhance_gain (0 disables enhancement)
        // CAUTION: testing shows the default 0.15 *fills* a clear notch (a notch is
        // concave-up, so the added +d2-of-Gaussian response raises it), collapsing
        // its prominence below the accept threshold. The spec flags this step as an
        // unvalidated default and says 0 is the correct setting for validation.
        // Recommend 0 until validated on real records; see notes in the .cpp.
        double dnMinProminence = 0.02;   // dn_min_prominence (Tier-1 accept threshold)
        int    iemMaxIter = 12;     // iem_max_iter
        double iemTol = 1e-4;   // iem_tol
        double dnWindowLoMs = 120.0;  // dn_window_lo_ms
        double dnWindowHiRrFrac = 0.70;   // dn_window_hi_rr_frac

        // UPPER BOUND FROM THE DIASTOLIC PEAK, when the caller knows it.
        //
        // dnWindowHiRrFrac puts the window's high edge at a fixed 70% of the RR
        // after the systolic peak, which is a proxy for "before the diastolic
        // peak" and a loose one: at a long RR it reaches well past the peak, so
        // the notch search sees the diastolic upstroke, its shoulder, and the
        // trough beyond it. Any local minimum of the IEM residual in that region
        // is a candidate, and one past the diastolic peak scores on prominence
        // just as well as the real notch -- so the detector could return a notch
        // AFTER the peak it is supposed to precede, and did.
        //
        // A caller that has already located the diastolic peak should pass it
        // here. The window then ends where the physiology says it ends, and the
        // notch-before-peak2 invariant holds by construction rather than by
        // correction afterwards. -1 (the default) keeps the RR-fraction
        // behaviour, so the detector is unchanged for callers that have no peak2
        // -- which includes anything calling it before the spline fit runs.
        //
        // Both bounds still apply: the effective high edge is the MINIMUM of
        // this and the RR fraction. Passing a peak2 can only narrow the window,
        // never widen it past what the RR allows.
        int dnWindowHiSample = -1;
    };

    struct IemResult {
        std::vector<double> profile;    // extracted non-stationary component
        std::vector<double> residual;   // input minus profile
        int  iterations = 0;
        bool converged = false;
    };

    struct DnResult {
        int    index = -1;                                  // sample in pulse; -1 absent
        double subSample = std::numeric_limits<double>::quiet_NaN();  // refined position
        enum Tier { IEM = 1, WINDKESSEL = 2, ABSENT = 3 } tier = ABSENT;
        double confidence = 0.0;                                 // normalized, 0..1
    };

    // Tier 1 primitive: iterate cur := mean(upper_env, lower_env) until the change
    // falls below tol (relative to the signal norm) or maxIter is reached. The
    // converged iterate is the slow profile; residual = input - profile.
    IemResult iterativeEnvelopeMean(const std::vector<double>& x,
        int maxIter = 12, double tol = 1e-4);

    // Multi-scale second-derivative-of-Gaussian enhancement of the descending limb,
    // to amplify a faint notch before IEM. gain = 0 disables (use when validating
    // against published error figures).
    std::vector<double> enhanceDescendingLimb(const std::vector<double>& pulse,
        double fs, int sysPeak, double gain = 0.15);

    // Tier 2: Windkessel flow-zero fallback within [winLo, winHi].
    DnResult windkesselDn(const std::vector<double>& pulse, double fs,
        int sysPeak, int winLo, int winHi);

    // Full three-tier detector. rrSeconds is the current RR interval.
    DnResult detectDicroticNotch(const std::vector<double>& pulse, double fs,
        int sysPeak, double rrSeconds, const PpgConfig& cfg);

} // namespace ppg_dicrotic