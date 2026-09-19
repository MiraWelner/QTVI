#pragma once
//
// anchor_view.hpp
// this handles the anchored alignments and bars in the grid


#include "anchor_type.hpp"

#include <array>
#include <cstring>   // std::strcmp, for markerForPoint / hasUserColumn

namespace anchor_view {

    // Every alignment the session holds, in the order sidecar CSVs are merged.
    inline constexpr std::array<AnchorType, 4> kAllAnchors = {
        AnchorType::R_PEAK,     // "R" -- the base; the grid always draws this one
        AnchorType::P_ONSET,    // "P"
        AnchorType::Q_ONSET,    // "Q"
        AnchorType::J_POINT     // "T"
    };

    // Operator-facing name. THIS IS THE CSV COLUMN SUFFIX -- the merged
    // markings CSV gets <col>_R / _P / _Q / _J -- so changing a string here
    // renames columns in every downstream analysis script.
    inline constexpr const char* label(AnchorType a) {
        switch (a) {
        case AnchorType::R_PEAK:  return "R";
        case AnchorType::P_ONSET: return "P";
        case AnchorType::Q_ONSET: return "Q";
        case AnchorType::J_POINT: return "J";
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

    // THE CANONICAL ALIGNMENT for a bar: which cell the AUTOMATIC view draws.
    // Automatic means the operator has not chosen an alignment, so the grid
    // assembles one bar per landmark and this says whose copy that is
    // (TemplateBin::userMarks). It is NOT the storage rule, NOT the seeding
    // rule and NOT the CSV rule any more -- all three are showsBar.
    // Meaningless for a glyph; returns R_PEAK there, and both predicates below
    // answer false, so a glyph can never be pinned to one alignment.
    inline constexpr AnchorType anchorFor(int marker) {
        switch (marker) {
        case kPBegin: return AnchorType::P_ONSET;
        case kQBegin: return AnchorType::Q_ONSET;
        case kSEnd:   return AnchorType::Q_ONSET;
        case kTEnd:   return AnchorType::Q_ONSET;
        }
        return AnchorType::R_PEAK;
    }

    // Is `a` the canonical alignment for this bar? Exactly one per bar, none
    // for any glyph. Read only by the Automatic assembly.
    inline constexpr bool owns(AnchorType a, int marker) {
        return isBar(marker) && anchorFor(marker) == a;
    }

    // ---- TWO KINDS OF ALIGNMENT ------------------------------------------
    //
    // P_ONSET and Q_ONSET OWN the canonical bars: p_begin lives in P, and
    // q_onset / s_end / t_end in Q (anchorFor). Automatic shows that set, and
    // forcing P or Q shows the part of it that alignment carries -- the SAME
    // bars, so one value, one colour, one set of output columns.
    //
    // R_PEAK and J_POINT are separately aligned averages. A P onset measured
    // on the R-aligned average is not the P-aligned one moved sideways; it is
    // a different measurement of a different waveform. So those two get their
    // own cells, drawn in their own style, reported in their own columns.
    inline constexpr bool ownsCanonicalBar(AnchorType a) {
        return a == anchorFor(kPBegin) || a == anchorFor(kQBegin);
    }

    // The complement: this alignment's bars are its own, not the canonical set.
    inline constexpr bool hasOwnBars(AnchorType a) {
        return !ownsCanonicalBar(a);
    }

    inline constexpr bool showsBar(AnchorType a, int marker) {
        if (!isBar(marker)) return false;
        switch (a) {
        case AnchorType::R_PEAK:  return true;
        case AnchorType::P_ONSET: return marker == kPBegin;
        case AnchorType::Q_ONSET: return marker == kQBegin
            || marker == kSEnd || marker == kTEnd;
        case AnchorType::J_POINT: return marker == kTEnd;
        }
        return false;
    }

    // A bar's canonical alignment must be one that offers it, or Automatic
    // would assemble from a cell nothing ever seeds.
    static_assert(showsBar(anchorFor(kPBegin), kPBegin), "");
    static_assert(showsBar(anchorFor(kQBegin), kQBegin), "");
    static_assert(showsBar(anchorFor(kSEnd), kSEnd), "");
    static_assert(showsBar(anchorFor(kTEnd), kTEnd), "");

    // The rows, pinned -- the guard landmark_admissibility.hpp used to carry.
    static_assert(showsBar(AnchorType::P_ONSET, kPBegin), "");
    static_assert(!showsBar(AnchorType::P_ONSET, kQBegin), "");
    static_assert(!showsBar(AnchorType::P_ONSET, kTEnd), "");
    static_assert(showsBar(AnchorType::Q_ONSET, kQBegin), "");
    static_assert(showsBar(AnchorType::Q_ONSET, kSEnd), "");
    static_assert(showsBar(AnchorType::Q_ONSET, kTEnd), "");
    static_assert(!showsBar(AnchorType::Q_ONSET, kPBegin), "");
    static_assert(showsBar(AnchorType::J_POINT, kTEnd), "");
    static_assert(!showsBar(AnchorType::J_POINT, kSEnd), "");
    static_assert(showsBar(AnchorType::R_PEAK, kPBegin), "");
    static_assert(showsBar(AnchorType::R_PEAK, kTEnd), "");
    static_assert(!showsBar(AnchorType::R_PEAK, kPPeak), "");

    // ---- THE CSV COLUMN RULE -----------------------------------------------
    //
    // Which BAR a markings-CSV point name refers to, or -1 when the name is a
    // glyph. Keyed by NAME rather than by id because two of the points --
    // q_peak and s_peak -- are computed inside the QRS and have no marker id
    // to key on, so the id-based predicates above cannot be asked about them.
    //
    // The four names here are exactly the four fields writeTemplateMarkingsBin
    // persists per (lead, slot, anchor). That is the cross-check: a landmark
    // absent from the .bin was never placed by hand, so it must not have a
    // _user column, and any name added here without a matching field in that
    // record is reporting an operator value that is not stored anywhere.
    inline int markerForPoint(const char* pointName) {
        if (std::strcmp(pointName, "p_begin") == 0) return kPBegin;
        if (std::strcmp(pointName, "q_onset") == 0) return kQBegin;
        if (std::strcmp(pointName, "s_end") == 0) return kSEnd;
        if (std::strcmp(pointName, "t_end") == 0) return kTEnd;
        return -1;
    }
    inline bool hasUserColumn(const char* pointName, AnchorType a) {
        const int m = markerForPoint(pointName);
        return m >= 0 && showsBar(a, m);
    }

} // namespace anchor_view