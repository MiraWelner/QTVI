#pragma once
//
// landmark_admissibility.hpp
//
// WHICH ECG LANDMARKS A GIVEN ALIGNMENT MAY REPORT.
//
// Alignment is not free. Shifting every beat so its P-onsets coincide leaves the
// R spikes scattered by the PR-interval spread, so the column-wise median smears
// what lies past that spread: the Q notch fills in, amplitudes drop, and a
// finder run on that template returns a position for a feature the averaging
// destroyed. The number is not wrong because the finder failed -- it is wrong
// because the waveform it measured is an artifact of the alignment.
//
// So each pass reports only what it may and stays silent about the rest. -1
// already means absent everywhere downstream (the paint loop skips it, emitPair
// blanks the CSV cell, BankMarkerSet defaults to it), so masking at the point of
// detection suppresses the bar, the glyph and the column together.
//
// ---------------------------------------------------------------------------
// WHY THIS IS A TABLE AND NOT A RULE
// ---------------------------------------------------------------------------
//
// The tempting design derives admissibility from interval rigidity: declare
// which inter-landmark intervals hold their shape under averaging, then admit a
// landmark iff every interval between it and the anchor is rigid. One
// physiological claim per interval, and a new anchor needs no new data.
//
// It does not fit. The set below has P_ONSET admitting R_PEAK across the PR
// segment while Q_ONSET does not admit P_PEAK across the same segment, and has
// Q_ONSET admitting R_PEAK but not S_END while J_POINT admits T_END -- i.e. the
// QRS is not treated as rigid end to end but the ST-T is. No single
// rigidity chain produces that, so a chain implementation would have to be bent
// with exceptions until the exceptions were the content.
//
// The table is therefore the primary source, stated once, per anchor. It is
// EDITORIAL: it records which landmarks this project chooses to trust from each
// pass, which is a decision about the measurement protocol rather than a
// derivable fact. Each row says what it admits; a row is the whole story for
// that anchor.
//
// ADDING AN ANCHOR IS ONE ROW.
//

// AnchorType. feature_marks.hpp must NOT include this header back -- include it
// from feature_marks.cpp only, or the graph cycles. If AnchorType is ever split
// into its own dependency-free header, switch this to that instead: it would let
// template_io.hpp type its raw_anchors map as AnchorType rather than int (see
// the note at its declaration explaining why it settled for int).
#include "feature_marks.hpp"

namespace landmark_admit {

    // Every landmark a pass can be asked about, in cardiac timeline order.
    // Order is presentational only here -- nothing walks the sequence -- but
    // keeping it chronological makes a row readable at a glance.
    enum class Landmark : int {
        P_ONSET = 0,
        P_PEAK  = 1,
        Q_ONSET = 2,
        Q_PEAK  = 3,
        R_PEAK  = 4,
        S_END   = 5,   // == J-point
        T_PEAK  = 6,
        T_END   = 7,
        Count   = 8
    };

    struct Mask {
        // Named to mirror FeatureMarks::TemplateLandmarks field for field, so a
        // call site reads as a one-to-one masking of it.
        bool p_begin = false;
        bool p_peak  = false;
        bool q_begin = false;
        bool q_peak  = false;
        bool r_peak  = false;
        bool s_end   = false;
        bool t_peak  = false;
        bool t_end   = false;

        constexpr bool operator()(Landmark lm) const {
            switch (lm) {
            case Landmark::P_ONSET: return p_begin;
            case Landmark::P_PEAK:  return p_peak;
            case Landmark::Q_ONSET: return q_begin;
            case Landmark::Q_PEAK:  return q_peak;
            case Landmark::R_PEAK:  return r_peak;
            case Landmark::S_END:   return s_end;
            case Landmark::T_PEAK:  return t_peak;
            case Landmark::T_END:   return t_end;
            default:                return false;
            }
        }
    };

    // Q-PEAK IS ALWAYS ADMITTED, ON EVERY PASS.
    //
    // Not because averaging cannot smear it -- it can -- but because it is
    // auto-only and never editable, and compute_q_peak already refuses to
    // fabricate one: it requires a STRICT INTERIOR minimum in [R-50ms, R] and a
    // depth of at least Q_MIN_DEPTH_MV below the PQ isoelectric median, and
    // returns -1 otherwise. A smeared Q fails the depth test on its own, so the
    // detector is its own mask here and a second one would only hide the cases
    // where it correctly succeeded.
    inline constexpr bool kQPeakAlwaysAdmitted = true;

    // R-PEAK IS ALWAYS ADMITTED, ON EVERY PASS.
    //
    // R is not measured from the aligned median at all. alignment.hpp carries
    // the passed-in R column through the median of the applied shifts and then
    // snaps it within +/-5 ms, so it is a bookkept reference column, not a
    // detection -- and it is the search origin every other finder takes. seed_all
    // already refuses to let it be absent, falling back to the unrefined
    // r_col_raw rather than -1, on the grounds that it is what all the other
    // landmarks are expressed against.
    //
    // NOTE this differs from the spec as given, which listed R under P_ONSET and
    // Q_ONSET but not under J_POINT. Suppressing it on the J pass alone would
    // contradict the fallback above and leave that pass's landmarks expressed
    // against a column the CSV does not report. Flip this to false and set the
    // rows individually if the omission was deliberate.
    inline constexpr bool kRPeakAlwaysAdmitted = true;

    inline constexpr Mask maskFor(AnchorType a) {
        Mask m;
        switch (a) {

        // Aligned on the P onset. The P wave is sharp; everything from the PR
        // segment onward is smeared by AV conduction variability.
        case AnchorType::P_ONSET:
            m.p_begin = true;
            m.p_peak  = true;
            break;

        // Aligned on the Q onset. The QRS upstroke is sharp; the J point and
        // beyond are not trusted from this pass.
        case AnchorType::Q_ONSET:
            m.q_begin = true;
            break;

        // Aligned on the J point. The whole ST-T segment is sharp here, which is
        // the reason this pass exists: T-end is unmeasurable from any earlier
        // anchor.
        //
        // T_PEAK needs no detector of its own -- it is the reactive glyph,
        // recomputed from the S_END and T_END bars by FeatureMarks::reactive_ecg,
        // so it follows a drag of either bracket. Admitting it here means
        // admitting the bracket it is derived from, which this row already does.
        case AnchorType::J_POINT:
            m.s_end  = true;
            m.t_peak = true;
            m.t_end  = true;
            break;

        // Aligned on R -- the base alignment, and the template every other pass
        // aligns FROM. Reports everything: it is what a plain R-only file
        // contains, and every consumer of the *_auto_ch fields has always been
        // able to expect a full set from it.
        //
        // This is a convention, not a claim that R-alignment smears nothing. It
        // does smear P and T, which is precisely why the other passes exist. But
        // blanking most of the base pass's columns is a far larger change than
        // the one this file is for, so the exemption is stated here rather than
        // hidden.
        case AnchorType::R_PEAK:
            m.p_begin = m.p_peak = m.q_begin = true;
            m.s_end = m.t_peak = m.t_end = true;
            break;
        }

        if (kQPeakAlwaysAdmitted) m.q_peak = true;
        if (kRPeakAlwaysAdmitted) m.r_peak = true;
        return m;
    }

    // The intended rows, pinned. A careless edit above fails here rather than in
    // a CSV weeks later.
    static_assert(maskFor(AnchorType::P_ONSET).p_begin, "");
    static_assert(maskFor(AnchorType::P_ONSET).p_peak, "");
    static_assert(maskFor(AnchorType::P_ONSET).r_peak, "");
    static_assert(!maskFor(AnchorType::P_ONSET).q_begin, "");
    static_assert(!maskFor(AnchorType::P_ONSET).s_end, "");
    static_assert(!maskFor(AnchorType::P_ONSET).t_end, "");

    static_assert(maskFor(AnchorType::Q_ONSET).q_begin, "");
    static_assert(maskFor(AnchorType::Q_ONSET).q_peak, "");
    static_assert(maskFor(AnchorType::Q_ONSET).r_peak, "");
    static_assert(!maskFor(AnchorType::Q_ONSET).p_peak, "");
    static_assert(!maskFor(AnchorType::Q_ONSET).s_end, "");

    static_assert(maskFor(AnchorType::J_POINT).s_end, "");
    static_assert(maskFor(AnchorType::J_POINT).t_peak, "");
    static_assert(maskFor(AnchorType::J_POINT).t_end, "");
    static_assert(!maskFor(AnchorType::J_POINT).q_begin, "");
    static_assert(!maskFor(AnchorType::J_POINT).p_begin, "");

}  // namespace landmark_admit
