#pragma once
/*
* anchor_view.hpp
* @brief handle the
*/

#include <array>
#include <cstring>   // std::strcmp, for markerForPoint / hasUserColumn

// T_END IS LAST ON PURPOSE. The numeric value is the tag every archive and
// markings file stores (raw_anchors, bank_anchors, markers_by_anchor), so the
// existing four keep the values they were written with.
enum class AnchorType { P_ONSET, Q_ONSET, R_PEAK, J_POINT, T_END };

namespace anchor_view {

    // Every alignment the session holds, in the order sidecar CSVs are merged.
    inline constexpr std::array<AnchorType, 5> anchor_array = { AnchorType::R_PEAK, AnchorType::P_ONSET, AnchorType::Q_ONSET, AnchorType::J_POINT, AnchorType::T_END };

    //defines the suffix for the column headers printed to the csv
    inline constexpr const char* label(AnchorType a) {
        switch (a) {
        case AnchorType::R_PEAK:  return "R";
        case AnchorType::P_ONSET: return "P";
        case AnchorType::Q_ONSET: return "Q";
        case AnchorType::J_POINT: return "J";
        case AnchorType::T_END:   return "T";
        }
        return "R";
    }

    // ---- HOW FAR EACH ALIGNMENT MAY MOVE A BEAT, AND WHAT IT LOOKS AT --------
    //
    // align_beat_matrix cross-correlates each beat against the bin template in
    // a window about that alignment's landmark, and takes the best lag within
    // +-max_lag. One pair for every landmark used to be the rule (+-30 ms, lag
    // up to the half-window), and the landmarks are not alike:
    //
    //   the window must hold the ONE transition being aligned and nothing else
    //   -- a window holding two locks onto the larger -- so it scales with how
    //   wide that transition is; and
    //
    //   the lag must cover how far that landmark really wanders beat to beat
    //   RELATIVE TO R, since every beat arrives R-aligned.
    //
    //   R   30 / 10  already the frame's origin; this only cleans sub-sample
    //                residue, so a small lag cannot walk it onto Q or S
    //   P   40 / 30  a slow, low onset; PR varies more than the QRS does
    //   Q   20 / 15  sharp, and 30 ms either side reaches the R upstroke
    //   J   30 / 20  the S-to-ST corner
    //   T  100 / 60  the T downslope is a slow transition tens of ms long, and
    //                QT tracks the preceding RR, so T end moves the furthest
    //
    // Milliseconds, converted at the call. Tune here; nothing else carries a
    // copy.
    struct ShiftWindow { double half_ms; double max_lag_ms; };
    inline constexpr ShiftWindow shiftWindow(AnchorType a) {
        switch (a) {
        case AnchorType::R_PEAK:  return { 30.0, 10.0 };
        case AnchorType::P_ONSET: return { 40.0, 30.0 };
        case AnchorType::Q_ONSET: return { 20.0, 15.0 };
        case AnchorType::J_POINT: return { 30.0, 20.0 };
        case AnchorType::T_END:   return { 100.0, 60.0 };
        }
        return { 30.0, 30.0 };
    }

    //this enum makes it easier to list the fiducial marker locations in order
    enum EcgMarker : int { p_begin = 0, p_peak = 1, q_begin = 2, r_peak = 3, j_point = 4, t_end = 5 };

    inline constexpr bool isEcgMarker(int m) { return m >= p_begin && m <= t_end; }

    inline constexpr bool isBar(int marker) {
        //only these 4 locations have user movable bars
        return marker == p_begin || marker == q_begin || marker == j_point || marker == t_end;
    }

    inline constexpr AnchorType anchorFor(int marker) {
        //right now we are focusing on PQ distance so p begin is p onset aligned while the rest are q aligned
        switch (marker) {
        case p_begin: return AnchorType::P_ONSET;
        case q_begin: return AnchorType::Q_ONSET;
        case j_point:   return AnchorType::Q_ONSET;
        case t_end:   return AnchorType::Q_ONSET;
        }
        return AnchorType::R_PEAK;
    }

    inline constexpr bool owns(AnchorType a, int marker) {
        //is a given algnment the canonical one right now for given marker
        return isBar(marker) && anchorFor(marker) == a;
    }

    inline constexpr bool hasOwnBars(AnchorType a) {
        //p and q alignments show the same bars that are shown in the auto alignment so they don't get their own bars. all others do.
        return !(a == anchorFor(p_begin) || a == anchorFor(q_begin));
    }

    inline constexpr bool showsBar(AnchorType a, int marker) {
        //different alignments show different bars - this shows which ones show which
        if (!isBar(marker)) return false;
        switch (a) {
        case AnchorType::R_PEAK:  return true;
        case AnchorType::P_ONSET: return marker == p_begin;
        case AnchorType::Q_ONSET: return marker == q_begin || marker == j_point || marker == t_end;
        case AnchorType::J_POINT: return marker == t_end;
            // T-aligned, T end is the landmark every beat is stacked on, so this
            // is where it is sharpest to place. Its own cell (hasOwnBars), like
            // J's; the Automatic view still reads t_end from its owner.
        case AnchorType::T_END:   return marker == t_end;
        }
        return false;
    }
    inline int markerForPoint(const char* pointName) {
        if (std::strcmp(pointName, "p_begin") == 0) return p_begin;
        if (std::strcmp(pointName, "q_onset") == 0) return q_begin;
        if (std::strcmp(pointName, "s_end") == 0) return j_point;
        if (std::strcmp(pointName, "t_end") == 0) return t_end;
        return -1;
    }

    // THE BARS THIS ALIGNMENT HAS -- showsBar, and nothing else:
    //
    //   _R : p_begin, q_onset, s_end, t_end
    //   _P : p_begin
    //   _Q : q_onset, s_end, t_end
    //   _J : t_end
    //   _T : t_end
    //
    // A user column exists where a bar exists. No column is emitted for a bar
    // the alignment does not have.
    //
    // WHAT THE COLUMN CONTAINS is where that bar IS: the operator's edit if
    // they moved it, the detection they left it on if they did not. See
    // finalBarsFor -- barsForPanel's rule, so the number in the file is the
    // number that was on screen.
    inline constexpr bool hasUserColumnFor(AnchorType a, int marker) {
        return showsBar(a, marker);
    }
    inline bool hasUserColumn(const char* pointName, AnchorType a) {
        const int m = markerForPoint(pointName);
        return m >= 0 && hasUserColumnFor(a, m);
    }

} // namespace anchor_view