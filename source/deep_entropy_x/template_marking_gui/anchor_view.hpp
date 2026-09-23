#pragma once
/*
* anchor_view.hpp
* @brief handle the
*/

#include <array>
#include <cstring>   // std::strcmp, for markerForPoint / hasUserColumn

enum class AnchorType { P_ONSET, Q_ONSET, R_PEAK, J_POINT };

namespace anchor_view {

    // Every alignment the session holds, in the order sidecar CSVs are merged.
    inline constexpr std::array<AnchorType, 4> anchor_array = { AnchorType::R_PEAK, AnchorType::P_ONSET, AnchorType::Q_ONSET, AnchorType::J_POINT };

    //defines the suffix for the column headers printed to the csv
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
        p_begin = 0,
        p_peak = 1,
        q_begin = 2,
        r_peak = 3,
        s_end = 4,
        t_end = 5,
        j_point = s_end   // the J point and the end of S are one landmark
    };
    inline constexpr bool isEcgMarker(int m) { return m >= p_begin && m <= t_end; }

    inline constexpr bool isBar(int marker) {
        //only these 4 locations have user movable bars
        return marker == p_begin || marker == q_begin || marker == s_end || marker == t_end;
    }
    inline constexpr bool isGlyph(int marker) {
        //the rest are just glyphs
        return isEcgMarker(marker) && !isBar(marker);
    }

    inline constexpr AnchorType anchorFor(int marker) {
        //right now we are focusing on PQ distance so p begin is p onset aligned while the rest are q aligned
        switch (marker) {
        case p_begin: return AnchorType::P_ONSET;
        case q_begin: return AnchorType::Q_ONSET;
        case s_end:   return AnchorType::Q_ONSET;
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
        case AnchorType::Q_ONSET: return marker == q_begin || marker == s_end || marker == t_end;
        case AnchorType::J_POINT: return marker == t_end;
        }
        return false;
    }
    inline int markerForPoint(const char* pointName) {
        if (std::strcmp(pointName, "p_begin") == 0) return p_begin;
        if (std::strcmp(pointName, "q_onset") == 0) return q_begin;
        if (std::strcmp(pointName, "s_end") == 0) return s_end;
        if (std::strcmp(pointName, "t_end") == 0) return t_end;
        return -1;
    }
    inline bool hasUserColumn(const char* pointName, AnchorType a) {
        const int m = markerForPoint(pointName);
        return m >= 0 && showsBar(a, m);
    }

} // namespace anchor_view