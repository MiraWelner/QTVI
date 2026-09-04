#pragma once
//
// anchor_view.hpp
//
// ONE SCREEN, FOUR ALIGNMENTS.
//
// The marking session used to run four sequential passes -- one window per
// alignment, "Finish and Next" between them -- and each pass could only see
// its own aligned template. This header replaces that with one rule, stated in
// the vocabulary BinPlotWidget.cpp already uses:
//
//   A GLYPH IS REPORTED ON EVERY ALIGNMENT. A BAR IS REPORTED ON ONE.
//
// A BAR is a draggable vertical line the operator positions. Its column is a
// column OF the aligned average it was placed against -- a P onset measured on
// the P-aligned average is not the same number as one measured on the R-aligned
// average -- so it is reported under that alignment's suffix and no other.
// There are four bars and four alignments, one each:
//
//   bar              close-up shows       column reported in
//   -------------------------------------------------------------
//   P onset          P-aligned            _P
//   Q onset          Q-aligned            _Q
//   J point (S end)  R-aligned            _R
//   T end            T-aligned            _T
//
// A GLYPH is a mark the widget draws and the operator cannot touch --
// markerAtX() never hit-tests one, so a click can neither select nor move it.
// It is a measurement, not a judgement, so it is measured independently on all
// four aligned averages and all four are reported. <landmark>_auto_P through
// _auto_T are four measurements of the same landmark on four waveforms, and
// comparing them is how the effect of an alignment on a landmark becomes
// visible. That is why the admissibility mask governs bars only.
//
// GLYPHS COME IN THE TWO FLAVOURS BinPlotWidget ALREADY NAMES, and the rule
// covers both the same way:
//
//   FROZEN   -- from detection on that alignment's own average (P peak from
//               the detector, R peak from r_col, Q onset, T end, ...). Four
//               alignments, four detections, four columns.
//   REACTIVE -- recomputed from the current bars: P peak between the P-onset
//               and Q-onset bars, T peak between the S-end and T-end bars,
//               Q peak and S peak inside the QRS. Their brackets STRADDLE
//               ALIGNMENTS -- P peak's two bars live on _P and _Q -- so a
//               reactive glyph has no single alignment home either, and is
//               likewise computed once per alignment from the assembled bar
//               set translated into that alignment's frame.
//
// R PEAK IS A GLYPH, NOT A BAR: the alignment anchor itself, re-derived from
// each template's own r_col at every load, never placed by hand, holding no
// BankMarkerSet field. Four auto columns, no user column anywhere. So is
// P PEAK -- reactive, per above -- which is why neither appears in the table.
//
// WHY "T-ALIGNED" IS AnchorType::J_POINT. The enum names each alignment after
// the landmark beats are shifted ONTO; the operator names it after the segment
// it is FOR. Aligning on the J point is what makes the whole ST-T segment
// sharp, and that pass exists for one reason: T-end is not measurable from any
// earlier anchor. It is labelled "T" on screen and in the CSV.
//
// MARKER IDS ARE DUPLICATED HERE ON PURPOSE, and they were RENUMBERED when
// T begin was removed -- 4 is S end, 5 is T end. Nothing persists a marker id,
// so the shift is invisible outside the process; the static_asserts in
// BinPlotWidget.cpp are what keep the two lists honest.
//
// BinPlotWidget.hpp includes
// template_marking_bin_io.hpp, which includes this file, so this header cannot
// include BinPlotWidget.hpp back. The literals below mirror
// BinPlotWidget::Marker and BinPlotWidget.cpp carries static_asserts that fail
// the build if the two ever drift.
//
// AnchorType ONLY. This header is included by template_marking_bin_io.hpp,
// landmark_admissibility.hpp and BinPlotWidget.cpp; pulling feature_marks.hpp
// in here would have made all of them depend on template_bank.hpp and
// annotation_types.hpp just to ask which alignment a bar belongs to.
#include "anchor_type.hpp"

#include <array>
#include <cstring>   // std::strcmp, for hasUserColumn

namespace anchor_view {

    // Every alignment the session holds, in the order sidecar CSVs are merged.
    inline constexpr std::array<AnchorType, 4> kAllAnchors = {
        AnchorType::R_PEAK,     // "R" -- the base; the grid always draws this one
        AnchorType::P_ONSET,    // "P"
        AnchorType::Q_ONSET,    // "Q"
        AnchorType::J_POINT     // "T"
    };

    // Operator-facing name. THIS IS THE CSV COLUMN SUFFIX -- the merged
    // markings CSV gets <col>_R / _P / _Q / _T -- so changing a string here
    // renames columns in every downstream analysis script.
    inline constexpr const char* label(AnchorType a) {
        switch (a) {
        case AnchorType::R_PEAK:  return "R";
        case AnchorType::P_ONSET: return "P";
        case AnchorType::Q_ONSET: return "Q";
        case AnchorType::J_POINT: return "T";
        }
        return "R";
    }

    // ---- ECG marker ids, mirroring BinPlotWidget::Marker -------------------
    enum EcgMarker : int {
        kPBegin = 0,
        kPPeak = 1,
        kQBegin = 2,
        kRPeak = 3,
        kSEnd = 4,
        kTEnd = 5
    };
    inline constexpr bool isEcgMarker(int m) { return m >= kPBegin && m <= kTEnd; }

    // ---- BAR OR GLYPH ------------------------------------------------------
    //
    // THE FOUR BARS, and the same four markerAtX() will hit-test. Everything
    // else in the ECG marker range is a glyph: P peak and T peak are reactive,
    // R peak is the alignment anchor.
    //
    // Keep this in step with markerAtX's skip list. A marker the widget lets
    // the operator drag but this function calls a glyph would take an edit and
    // report it on all four alignments, none of which the operator was looking
    // at.
    inline constexpr bool isBar(int marker) {
        return marker == kPBegin || marker == kQBegin
            || marker == kSEnd || marker == kTEnd;
    }
    inline constexpr bool isGlyph(int marker) {
        return isEcgMarker(marker) && !isBar(marker);
    }

    // Which alignment this BAR is displayed on, stored against, and reported
    // under -- one question, one answer. Meaningless for a glyph; returns
    // R_PEAK there, and `owns` still answers false, so a glyph can never be
    // pinned to one alignment by accident.
    inline constexpr AnchorType anchorFor(int marker) {
        switch (marker) {
        case kPBegin: return AnchorType::P_ONSET;
        case kQBegin: return AnchorType::Q_ONSET;
        case kSEnd:   return AnchorType::R_PEAK;
        case kTEnd:   return AnchorType::J_POINT;
        }
        return AnchorType::R_PEAK;
    }

    // Does alignment `a` report a *_user value for this marker? True for
    // exactly one alignment per bar, and for no alignment for any glyph.
    inline constexpr bool owns(AnchorType a, int marker) {
        return isBar(marker) && anchorFor(marker) == a;
    }

    // ---- THE CSV COLUMN RULE -----------------------------------------------
    //
    // Does this point emit a *_user column at all? Keyed by the CSV column
    // name rather than a marker id because two of the points -- q_peak and
    // s_peak -- are computed inside the QRS and have no marker id to key on.
    //
    // Only r_peak answers false. It is the alignment anchor: re-derived from
    // each template's own r_col at every load, never placed by hand, holding no
    // BankMarkerSet field. Every other glyph gets a user column, because a
    // reactive glyph has a real second value -- the same landmark measured
    // between the operator's bars instead of between the detector's -- and that
    // second value is the one drawn on screen.
    //
    // writeTemplateMarkingsCsv's header emitter and its row loop BOTH call
    // this. They used to encode the rule twice, as a strcmp against "r_peak"
    // in one and a bare `k != 3` index test in the other, so adding a point
    // column meant finding both or writing a header that did not match its
    // rows.
    inline bool hasUserColumn(const char* pointName) {
        return std::strcmp(pointName, "r_peak") != 0;
    }

} // namespace anchor_view