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
// So each alignment asks the operator to place only the landmark it can show
// honestly, and stays silent about the rest. -1 already means absent everywhere
// downstream (the paint loop skips it, emitPair blanks the CSV cell,
// BankMarkerSet defaults to it), so masking at the point of seeding suppresses
// the bar and its *_user column together.
//
// IT SUPPRESSES THE BAR ONLY. The glyphs -- the *_auto_ch fields and their
// *_auto columns -- are NOT masked: every alignment detects every landmark on
// its own average and all four are reported, because the smearing described
// above is a thing worth measuring rather than hiding, and comparing one
// landmark's four *_auto columns is how it becomes visible. See the note above
// maskFor.
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

// AnchorType now comes from anchor_type.hpp, the dependency-free header this
// comment used to ask for. feature_marks.hpp is no longer included here at all,
// so the "must NOT include this header back or the graph cycles" hazard is
// gone: there is nothing left in this file's include closure for feature_marks
// to collide with. (template_io.hpp can now type its raw_anchors map as
// AnchorType rather than int -- see the note at its declaration explaining why
// it settled for int. Left alone here; it is a format-adjacent change.)
#include "anchor_type.hpp"
#include "anchor_view.hpp"

namespace landmark_admit {

    // Every landmark a pass can be asked about, in cardiac timeline order.
    // Order is presentational only here -- nothing walks the sequence -- but
    // keeping it chronological makes a row readable at a glance.
    enum class Landmark : int {
        P_ONSET = 0,
        P_PEAK = 1,
        Q_ONSET = 2,
        Q_PEAK = 3,
        R_PEAK = 4,
        S_END = 5,   // == J-point
        T_PEAK = 6,
        T_END = 7,
        Count = 8
    };

    struct Mask {
        // Named to mirror FeatureMarks::TemplateLandmarks field for field, so a
        // call site reads as a one-to-one masking of it.
        bool p_begin = false;
        bool p_peak = false;
        bool q_begin = false;
        bool q_peak = false;
        bool r_peak = false;
        bool s_end = false;
        bool t_peak = false;
        bool t_end = false;

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
    // (kQPeakAlwaysAdmitted DELETED with kRPeakAlwaysAdmitted: q_peak is a
    //  glyph, computed inside the QRS, and glyphs are no longer masked at all,
    //  so there is nothing left to exempt it from.)

    // ---- THE MASK GOVERNS BARS, NOT GLYPHS ------------------------------
    //
    // It used to govern both, and that was the larger of the two mistakes in
    // this file. A glyph is a MEASUREMENT: every alignment detects every
    // landmark on its own average and all four are reported, because comparing
    // <landmark>_auto_P against ..._auto_Q is precisely how the effect of an
    // alignment on a landmark becomes visible. Masking the *_auto_ch fields
    // blanked three quarters of that comparison in the name of not trusting a
    // smeared position -- but the smearing is the thing being measured.
    //
    // A BAR is a JUDGEMENT, and that is what this mask is for: it decides
    // which alignment the operator is asked to place a landmark on, and
    // therefore which alignment's *_user column reports it. Exactly one
    // alignment per bar, four bars, four alignments.
    //
    // (kRPeakAlwaysAdmitted DELETED. Its note recorded that the spec listed R
    // under P_ONSET and Q_ONSET but not J_POINT, and that suppressing R on the
    // J pass alone would leave that pass's landmarks expressed against a
    // column the CSV never reported. R is a glyph. With the mask off the glyph
    // side, every alignment reports it, the J pass included, and there is
    // nothing left to reconcile.)
    //
    // ONE TABLE, NOT TWO. The rows below are derived from anchor_view::owns so
    // that the CSV writer, the drag handler, the focus panel and this mask
    // cannot disagree about which alignment a bar belongs to. Change
    // anchor_view::anchorFor, not this function.
    inline constexpr Mask maskFor(AnchorType a) {
        Mask m;
        m.p_begin = anchor_view::owns(a, anchor_view::kPBegin);
        m.q_begin = anchor_view::owns(a, anchor_view::kQBegin);
        m.s_end = anchor_view::owns(a, anchor_view::kSEnd);
        m.t_end = anchor_view::owns(a, anchor_view::kTEnd);

        // p_peak / q_peak / t_peak stay FALSE on every alignment. They are
        // glyphs -- p_peak and t_peak reactive from the bars above, q_peak
        // computed inside the QRS -- so nothing seeds a stored copy of one,
        // and the only p_peak that exists is the one
        // TemplateBin::syncReactiveGlyphs recomputes from the bars. r_peak is
        // likewise false here: it has no BankMarkerSet field and no *_user
        // column on any alignment.
        return m;
    }

    // The intended rows, pinned: one bar each, and no glyph admitted anywhere.
    // A careless edit above fails here rather than in a CSV weeks later.
    static_assert(maskFor(AnchorType::P_ONSET).p_begin, "");
    static_assert(!maskFor(AnchorType::P_ONSET).q_begin, "");
    static_assert(!maskFor(AnchorType::P_ONSET).s_end, "");
    static_assert(!maskFor(AnchorType::P_ONSET).t_end, "");

    static_assert(maskFor(AnchorType::Q_ONSET).q_begin, "");
    static_assert(!maskFor(AnchorType::Q_ONSET).p_begin, "");
    static_assert(!maskFor(AnchorType::Q_ONSET).s_end, "");
    static_assert(!maskFor(AnchorType::Q_ONSET).t_end, "");

    // The J point is PLACED against the R-aligned average, where the QRS is
    // sharp, so R owns s_end -- and R owns nothing else. This is the row that
    // changed: R used to report everything (it had to, when an R-only file was
    // all a first pass produced) and J_POINT used to own s_end alongside t_end.
    static_assert(maskFor(AnchorType::R_PEAK).s_end, "");
    static_assert(!maskFor(AnchorType::R_PEAK).p_begin, "");
    static_assert(!maskFor(AnchorType::R_PEAK).q_begin, "");
    static_assert(!maskFor(AnchorType::R_PEAK).t_end, "");

    // The J alignment exists to make T-end measurable, and T-end is what it
    // reports.
    static_assert(maskFor(AnchorType::J_POINT).t_end, "");
    static_assert(!maskFor(AnchorType::J_POINT).s_end, "");
    static_assert(!maskFor(AnchorType::J_POINT).q_begin, "");
    static_assert(!maskFor(AnchorType::J_POINT).p_begin, "");

    // No glyph is admitted on any alignment.
    static_assert(!maskFor(AnchorType::P_ONSET).p_peak, "");
    static_assert(!maskFor(AnchorType::J_POINT).t_peak, "");
    static_assert(!maskFor(AnchorType::R_PEAK).r_peak, "");
    static_assert(!maskFor(AnchorType::Q_ONSET).q_peak, "");

}  // namespace landmark_admit