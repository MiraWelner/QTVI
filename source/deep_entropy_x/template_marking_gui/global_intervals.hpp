#pragma once
//
// global_intervals.hpp
//
// GLOBAL (across-lead) ECG intervals, reduced from the per-channel landmarks
// already carried on a TemplateBin.
//
// Why global rather than per-lead: a single ECG lead is a projection of the
// heart's electrical vector onto one axis, so a wavefront moving nearly
// perpendicular to that axis produces almost no deflection in it. The first
// few milliseconds of depolarization are therefore invisible in some leads and
// the last few are invisible in others, and EVERY lead reports a QRS that is
// too short by a lead-specific amount. The standard measurement takes the
// union: QRS onset is the EARLIEST onset seen in any lead, QRS offset is the
// LATEST offset seen in any lead. No single lead is trusted to contain the
// whole complex. Per-lead values are kept beside the global ones so a lead
// that disagrees wildly (usually a bad fiducial, not physiology) is visible.
//
// ---------------------------------------------------------------------------
// ALIGNMENT: the reason this is not a min/max over raw columns
// ---------------------------------------------------------------------------
// Each channel is aligned independently -- ch1/ch2/ch3 each get their own
// r_col_raw out of their own alignment pass -- so a sample column means
// something different in each channel's template. Column 100 in ch1 and
// column 100 in ch2 are NOT the same instant. Taking min/max of q_begin_ch[]
// directly would be comparing three different clocks, and the answer would
// silently absorb the inter-channel R-column offset.
//
// Everything here is therefore reduced in R-RELATIVE sample offsets
// (landmark - that channel's own R column). R is the one landmark all three
// channels share by construction, so it is the common origin. Onsets come out
// negative (before R), offsets positive.
//
// ---------------------------------------------------------------------------
// Conventions, matching the rest of the marking code
// ---------------------------------------------------------------------------
//   Positions : sample columns / offsets, sub-sample doubles. -1 = not found,
//               so a negative RAW column must never be arithmetic'd on before
//               checking. (R-relative offsets are legitimately negative --
//               those use the `valid` flags, not a sign test.)
//   Intervals : milliseconds, NAN when not computable. Same as
//               EcgFeatures::qrs_ms / qt_ms.
//   Channels  : 0..2 == ch1/ch2/ch3, the same indexing as the *_auto_ch[]
//               arrays and MarkerSet::*_ch[].
//
// Landmark names follow this codebase: q_begin is the QRS onset, s_end is the
// J point (AnchorType::J_POINT == S_END), p_begin is the P onset, t_end is the
// T offset.
//

#include "template_marking_bin_io.hpp"   // TemplateBin, AnchorType

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

namespace global_intervals {

    inline constexpr int    kNumEcgCh = 3;
    inline constexpr double kNotFound = -1.0;

    /// Which set of landmark positions to reduce. The auto values live in
    /// TemplateBin::*_auto_ch[] as sub-sample doubles; the user values live in
    /// the per-anchor MarkerSet as ints, since a dragged bar lands on a column.
    enum class MarkerSource { AUTO, USER };

    /**
     * @brief The landmarks of ONE channel, pulled off the bin as raw sample
     *        columns (-1 = not found) plus that channel's own R column.
     *
     *        rPeak is carried alongside deliberately: it is what converts the
     *        rest of this struct out of the channel's private clock and into
     *        the shared R-relative one. A LeadMarkers with rPeak < 0 cannot
     *        contribute to a global interval at all, no matter how many other
     *        landmarks it has, because there is no way to place them in time
     *        relative to the other channels.
     */
    struct LeadMarkers {
        double pOnset = kNotFound;   ///< p_begin -- P wave onset
        double qOnset = kNotFound;   ///< q_begin -- QRS onset (Q onset, or R onset if no Q)
        double rPeak = kNotFound;   ///< R column for THIS channel; the common origin
        double jPoint = kNotFound;   ///< s_end -- QRS offset / J point
        double tOffset = kNotFound;  ///< t_end -- T wave offset

        bool   hasR()      const { return rPeak >= 0.0; }
        /// Landmark as an offset from this channel's R. Caller must have
        /// checked both are found.
        double rel(double landmark) const { return landmark - rPeak; }
    };

    /**
     * @brief Extract one channel's landmarks from a bin.
     *
     * @param ch      0..2.
     * @param anchor  Which per-anchor user MarkerSet to read. Marker columns
     *                are only meaningful relative to the alignment they were
     *                placed on, so the anchor is not optional even for AUTO
     *                (it selects nothing there, but keeping one signature
     *                stops the two paths from drifting apart).
     */
    inline LeadMarkers leadMarkersFor(const TemplateBin& bin, int ch,
        AnchorType anchor, MarkerSource src) {
        LeadMarkers m;
        if (ch < 0 || ch >= kNumEcgCh) return m;

        // R is auto-only and not per-anchor. Prefer the sub-sample value, fall
        // back to the integer column.
        m.rPeak = (bin.r_peak_auto_ch[ch] >= 0.0)
            ? bin.r_peak_auto_ch[ch]
            : static_cast<double>(bin.r_peak_ch[ch]);

        if (src == MarkerSource::AUTO) {
            m.pOnset = bin.p_begin_auto_ch[ch];
            m.qOnset = bin.q_begin_auto_ch[ch];
            m.jPoint = bin.s_end_auto_ch[ch];
            m.tOffset = bin.t_end_auto_ch[ch];
        }
        else {
            const TemplateBin::MarkerSet& u = bin.marks(anchor);
            m.pOnset = static_cast<double>(u.p_begin_ch[ch]);
            m.qOnset = static_cast<double>(u.q_begin_ch[ch]);
            m.jPoint = static_cast<double>(u.s_end_ch[ch]);
            m.tOffset = static_cast<double>(u.t_end_ch[ch]);
        }
        return m;
    }

    /// All three channels at once. Channels with no R are still returned (so
    /// the caller can see they were considered) but contribute nothing.
    inline std::map<int, LeadMarkers> allLeadMarkers(const TemplateBin& bin,
        AnchorType anchor, MarkerSource src) {
        std::map<int, LeadMarkers> out;
        for (int c = 0; c < kNumEcgCh; ++c) out[c] = leadMarkersFor(bin, c, anchor, src);
        return out;
    }

    /**
     * @brief Intervals measured across ALL leads at once.
     *
     *        qrsOnset / qrsOffset are R-RELATIVE sample offsets (see the
     *        alignment note at the top), not columns in any one channel's
     *        template. onsetLead / offsetLead name the channel that supplied
     *        each extreme, which is what makes the result auditable: if one
     *        channel always wins both, its fiducials are worth a look.
     *
     *        `valid` is false when the QRS window could not be established;
     *        the ms fields are then NAN rather than a difference of sentinels.
     */
    struct GlobalIntervals {
        double qrsOnset = std::numeric_limits<double>::quiet_NaN();   ///< earliest q_begin, samples from R
        double qrsOffset = std::numeric_limits<double>::quiet_NaN();  ///< latest s_end, samples from R
        double qrsDuration_ms = std::numeric_limits<double>::quiet_NaN();
        double prInterval_ms = std::numeric_limits<double>::quiet_NaN();
        double qtInterval_ms = std::numeric_limits<double>::quiet_NaN();

        std::map<int, double> perLeadQrsOnset;    ///< per-lead, samples from R
        std::map<int, double> perLeadQrsOffset;   ///< per-lead, samples from R

        int  onsetLead = -1;   ///< channel that supplied the global onset
        int  offsetLead = -1;  ///< channel that supplied the global offset
        bool valid = false;    ///< QRS onset AND offset both established
    };

    /**
     * @brief Reduce per-channel landmarks to the global intervals.
     *
     * @param leadMarkers  Per-channel landmarks, from allLeadMarkers().
     *                     Channels may be missing, may lack R, or may have all
     *                     landmarks unset -- all three are simply "no
     *                     evidence" and are skipped identically.
     * @param rateHz       ECG sample rate, for the ms conversions. <= 0 leaves
     *                     every ms field NAN while the sample-offset fields
     *                     still populate.
     */
    inline GlobalIntervals computeGlobalIntervals(
        const std::map<int, LeadMarkers>& leadMarkers, double rateHz) {
        GlobalIntervals g;

        // Extremes tracked in locals seeded with the opposite infinity, copied
        // out only if a channel actually moved them. This is what keeps a
        // "nothing was found" call from returning lowest() - max() as a
        // duration -- the failure mode of reading the struct's own sentinels
        // back as if they were measurements.
        double onset = std::numeric_limits<double>::max();
        double offset = std::numeric_limits<double>::lowest();
        double pOnset = std::numeric_limits<double>::max();
        double tOffset = std::numeric_limits<double>::lowest();

        for (const auto& [lead, marks] : leadMarkers) {
            // No R column => no shared clock => this channel cannot be
            // compared with the others at all.
            if (!marks.hasR()) continue;

            if (marks.qOnset >= 0.0) {
                const double rel = marks.rel(marks.qOnset);
                g.perLeadQrsOnset[lead] = rel;
                if (rel < onset) { onset = rel; g.onsetLead = lead; }
            }
            if (marks.jPoint >= 0.0) {
                const double rel = marks.rel(marks.jPoint);
                g.perLeadQrsOffset[lead] = rel;
                if (rel > offset) { offset = rel; g.offsetLead = lead; }
            }
            if (marks.pOnset >= 0.0)  pOnset = std::min(pOnset, marks.rel(marks.pOnset));
            if (marks.tOffset >= 0.0) tOffset = std::max(tOffset, marks.rel(marks.tOffset));
        }

        const bool haveOnset = !g.perLeadQrsOnset.empty();
        const bool haveOffset = !g.perLeadQrsOffset.empty();
        if (haveOnset)  g.qrsOnset = onset;
        if (haveOffset) g.qrsOffset = offset;

        const double msPerSamp = (rateHz > 0.0) ? 1000.0 / rateHz : std::numeric_limits<double>::quiet_NaN();

        // QRS needs both edges. The two are located independently and in
        // different channels, so a garbled beat can produce an offset that
        // precedes the onset. That is a detection failure, not a
        // negative-width QRS: reported as NAN rather than as a negative number
        // that would quietly poison any average taken over bins.
        if (haveOnset && haveOffset && offset > onset) {
            g.qrsDuration_ms = (offset - onset) * msPerSamp;
            g.valid = true;
        }

        // PR and QT are both anchored on the GLOBAL QRS onset, so neither is
        // computable without it. Note QT runs onset -> T offset, matching
        // computeEcgFeatures' qt_ms (q_begin -> t_end), not R -> T offset.
        if (haveOnset && pOnset != std::numeric_limits<double>::max() && onset > pOnset)
            g.prInterval_ms = (onset - pOnset) * msPerSamp;
        if (haveOnset && tOffset != std::numeric_limits<double>::lowest() && tOffset > onset)
            g.qtInterval_ms = (tOffset - onset) * msPerSamp;

        return g;
    }

    /// Convenience overload: straight from a bin.
    inline GlobalIntervals computeGlobalIntervals(const TemplateBin& bin,
        AnchorType anchor, double rateHz, MarkerSource src = MarkerSource::USER) {
        return computeGlobalIntervals(allLeadMarkers(bin, anchor, src), rateHz);
    }

    /**
     * @brief Largest disagreement between any lead's QRS onset and the global
     *        onset, in ms. Expected to be small (a few ms up to a few tens);
     *        a large value flags a suspect fiducial rather than physiology.
     *
     * @return NAN when there is nothing to compare.
     */
    inline double onsetSpread_ms(const GlobalIntervals& g, double rateHz) {
        if (g.perLeadQrsOnset.empty() || std::isnan(g.qrsOnset) || rateHz <= 0.0)
            return std::numeric_limits<double>::quiet_NaN();
        double worst = 0.0;
        for (const auto& [lead, rel] : g.perLeadQrsOnset)
            worst = std::max(worst, rel - g.qrsOnset);
        return worst * (1000.0 / rateHz);
    }

}  // namespace global_intervals#pragma once
