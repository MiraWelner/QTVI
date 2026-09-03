#pragma once
//
// ppg_dicrotic.hpp
//
// E-5: three-tier dicrotic-notch (DN) detection (DeepEntropyX Phase 2, Section 6.3).
//
// ============================================================================
// CURRENT STATE: PLACEHOLDER. NO TIER IS IMPLEMENTED.
// ============================================================================
//
// detectDicroticNotch returns a fixed offset after the systolic peak, flagged
// PLACEHOLDER, with confidence 0. It is NOT a measurement and must not be read
// as one -- DnResult::measured() is the only correct test, and it returns false
// for every result this module currently produces.
//
// The offset is deliberately unphysiological (10 ms; the real notch is aortic
// valve closure at end-systole, 150-350 ms after the peak) so that a placeholder
// reads as obviously provisional on screen rather than as a plausible wrong
// answer. It sits BEFORE the tier search window on purpose: a placeholder that
// fell inside the window would look like something a tier had found.
//
// THE POINT OF THIS FILE IS THE SEAM, NOT THE ALGORITHM. The signature, the
// config, the result struct, and the window arithmetic are all final; only the
// three tier bodies are missing. Filling them in should touch nothing outside
// this module.
//
// ---------------------------------------------------------------------------
// THE PLAN
// ---------------------------------------------------------------------------
//
// The DN marks aortic valve closure. It is prominent in young compliant vessels
// and progressively damped with age/stiffness/vasodilation/low perfusion, to the
// point of vanishing. Simple derivative-zero methods fail on damped pulses, so
// the pipeline runs three tiers in order and records which one produced the
// answer:
//
//   Tier 1  Iterative Envelope Mean (IEM)  -- primary; detrend the pulse so the
//           notch is a clear local minimum of the residual even when it is only
//           an inflection in the raw pulse. (Pal et al. 2024.)
//   Tier 2  Windkessel-augmented fallback  -- physics-aware; reconstruct a
//           flow-proportional signal and take its zero-crossing. (Saffarpour 2023.)
//   Tier 3  ABSENT                         -- notch genuinely absent; recorded as
//           a feature, never substituted with a nominal position.
//
// Search window (all tiers): dnWindowLoMs after the systolic peak, to
// dnWindowHiRrFrac of the current RR after the systolic peak, further narrowed
// by dnWindowHiSample when the caller knows the diastolic peak.
//
// PHASE-1 DEPENDENCIES for the eventual implementation: splineThrough (Task N,
// cubic spline) and refineAsymmetric (Task I-3, sub-sample refinement). Neither
// is needed by the placeholder.
//

#include <cmath>
#include <limits>
#include <vector>

namespace ppg_dicrotic {

    // DN-related config fields (Section 6.3). May be merged into a broader PpgConfig.
    //
    // Every field except dnPlaceholderMs is unused at the placeholder stage and
    // is retained because it is the tiers' contract with the caller. Keeping the
    // struct stable means enabling a tier later is a change to this file alone.
    struct PpgConfig {
        // ---- placeholder ----
        // Offset after the systolic peak for the placeholder position. Not a
        // physiological estimate; see the header note.
        double dnPlaceholderMs = 10.0;

        // ---- Tier 1 (IEM), unused until implemented ----
        double dnEnhanceGain = 0.0;
        // CAUTION, carried forward from the E-5 notes: gain 0.15 *fills* a clear
        // notch (a notch is concave-up, so the added +d2-of-Gaussian response
        // raises it), collapsing its prominence below the accept threshold. The
        // spec flags this as an unvalidated default and says 0 is correct for
        // validation. Default is 0 here rather than 0.15 so the first run of a
        // real Tier 1 is the validatable one.
        double dnMinProminence = 0.02;   // Tier-1 accept threshold
        int    iemMaxIter = 12;
        double iemTol = 1e-4;

        // ---- search window, shared by all tiers ----
        double dnWindowLoMs = 120.0;
        double dnWindowHiRrFrac = 0.70;

        // UPPER BOUND FROM THE DIASTOLIC PEAK, when the caller knows it.
        //
        // dnWindowHiRrFrac puts the high edge at a fixed 70% of the RR after the
        // systolic peak, which is a proxy for "before the diastolic peak" and a
        // loose one: at a long RR it reaches well past that peak, so the search
        // would see the diastolic upstroke, its shoulder, and the trough beyond.
        // A minimum of the residual past the diastolic peak scores on prominence
        // just as well as the real notch, so a tier could return a notch AFTER
        // the peak it is supposed to precede.
        //
        // A caller that has already located the diastolic peak should pass it
        // here; the window then ends where the physiology says it ends and
        // notch < peak2 holds by construction rather than by correction. -1 keeps
        // the RR-fraction behaviour.
        //
        // Both bounds apply: the effective high edge is the MINIMUM of this and
        // the RR fraction. Passing a peak2 can only narrow the window.
        //
        // PASS ONLY A MEASURED PEAK2, never a placeholder. Bounding a real
        // detector by a fabricated position lets the fabrication suppress a
        // genuine notch.
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
        // PLACEHOLDER is distinct from ABSENT on purpose. ABSENT means a tier ran
        // and found nothing, which is a finding. PLACEHOLDER means no tier ran at
        // all. Collapsing the two would make the placeholder stage invisible in
        // any archive written during it.
        enum Tier { PLACEHOLDER = 0, IEM = 1, WINDKESSEL = 2, ABSENT = 3 } tier = PLACEHOLDER;
        double confidence = 0.0;                                 // normalized, 0..1

        // THE ONLY CORRECT TEST for "is this a measurement". Callers must not
        // write `tier != ABSENT` -- that is true for PLACEHOLDER, which is how a
        // provisional position would get drawn as a solid X and averaged into
        // downstream statistics.
        bool measured() const { return tier == IEM || tier == WINDKESSEL; }
    };

    // Resolved search window for the tiers, in samples. Exposed because it is
    // testable on its own and because a caller may want to know whether there was
    // any window at all: valid() is false when the pulse is too short, the RR too
    // fast, or a supplied peak2 too early for a notch to exist between the two.
    struct DnWindow {
        int  lo = -1;
        int  hi = -1;
        bool valid() const { return lo >= 0 && hi > lo; }
    };

    DnWindow dnSearchWindow(int nSamples, double fs, int sysPeak,
        double rrSeconds, const PpgConfig& cfg);

    // ---- Tier slots. NONE IMPLEMENTED; each returns an explicit not-implemented
    // ---- result rather than a guess. See the .cpp.
    IemResult iterativeEnvelopeMean(const std::vector<double>& x,
        int maxIter = 12, double tol = 1e-4);

    std::vector<double> enhanceDescendingLimb(const std::vector<double>& pulse,
        double fs, int sysPeak, double gain = 0.0);

    DnResult windkesselDn(const std::vector<double>& pulse, double fs,
        int sysPeak, int winLo, int winHi);

    // Full detector. rrSeconds is the current RR interval.
    //
    // AT THE PLACEHOLDER STAGE this ignores the pulse shape entirely and returns
    // sysPeak + cfg.dnPlaceholderMs, tier PLACEHOLDER, confidence 0. It still
    // returns index -1 when there is no room for the placeholder in the array, so
    // "absent" remains representable.
    DnResult detectDicroticNotch(const std::vector<double>& pulse, double fs,
        int sysPeak, double rrSeconds, const PpgConfig& cfg);

} // namespace ppg_dicrotic