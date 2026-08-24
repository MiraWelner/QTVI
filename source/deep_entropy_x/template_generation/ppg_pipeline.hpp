#pragma once
//
// ppg_pipeline.hpp
//
// Continuous-recording PPG signal conditioning. Three stages:
//
//   iemEnvelope    per-template Iterative Envelope Mean baseline (data-adaptive;
//                  splines an upper/lower envelope through a single pulse's own
//                  extrema and iterates the mean-subtraction to leave the slow
//                  profile underneath).
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

    // ---- Iterative Envelope Mean (per template pulse) -------------------
    struct IemEnvelope {
        std::vector<double> envelope;   // slow trend underneath (pulse - imf)
        std::vector<double> imf;        // extracted pulsatile component
        std::vector<double> upper;      // last-iteration upper envelope
        std::vector<double> lower;      // last-iteration lower envelope
        int  iterations = 0;
        bool ok = false;
    };

    // maxIter caps the sift; sdThresh is the Huang stop criterion
    // (sum (h_prev - h)^2 / sum h_prev^2 < sdThresh; 0.2-0.3 is classic).
    IemEnvelope iemEnvelope(const std::vector<double>& pulse,
        int maxIter = 12, double sdThresh = 0.2);

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