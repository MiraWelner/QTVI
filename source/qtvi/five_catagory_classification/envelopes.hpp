#pragma once
//
// envelopes.hpp
//
// Spec 4.7: moving-average, power-spectral and wavelet envelopes of the QRS,
// P and T waves over a short (5-10 beat) and a long (30 beat) window. Short
// windows track PVC onset and conduction change; long windows track ischaemic
// ST drift. Both feed the dynamic master template selection of spec 9.5.
//
// STATUS: stubs.
//
// ---------------------------------------------------------------------------
// READ BEFORE IMPLEMENTING. The spec does not define what an "envelope" is
// numerically, and the three kinds x three segments x two window lengths comes
// to 18 series per lead per bin. Guessing wrong here is a lot of wasted work,
// so the open questions are written down rather than silently resolved:
//
//   * MOVING_AVERAGE -- is the envelope a per-sample min/max band across the
//     window's beats, mean +/- k*std, or an inter-quartile band? A band needs
//     two curves (lo/hi); a smoothed mean needs one. EnvelopeBand below carries
//     three (lo/mid/hi) so any of those fits, but only one is intended.
//
//   * POWER_SPECTRAL -- "power-spectral envelope of a P wave over 5 beats" is
//     the least defined of the three. A P wave is ~100 samples at 1000 Hz,
//     which is thin for a PSD. Unresolved: concatenate the window's segments
//     into one series and transform, or transform each and average the spectra?
//     And is the output a spectrum (indexed by frequency) or a scalar per beat
//     (e.g. band-power ratio) tracked over time? Those are different shapes.
//
//   * WAVELET -- no basis, no level count, no decomposition given. There is no
//     wavelet code anywhere in this repo yet, so this one is a new dependency
//     as well as a new algorithm.
//
//   * Alignment -- envelopes across beats only mean anything if the beats are
//     column-aligned first. alignment.hpp's ecg_beat_set is already that, so
//     these should consume aligned beats, never raw slices.
//
//   * Segment bounds -- P/QRS/T spans come from the FeatureMarks landmarks. On
//     a per-beat basis those move, so either the bounds are taken from the
//     reference template (fixed columns, simpler) or per beat (follows the
//     morphology, but then the window's members aren't sample-comparable).
//     Fixed columns from the reference is the assumption below.
// ---------------------------------------------------------------------------

#include <vector>
#include <cstdint>
#include <cmath>
#include <limits>
#include <algorithm>

namespace envelopes {

    // ---- What is being tracked --------------------------------------------
    enum class Segment : uint8_t { P, QRS, T };

    inline const char* segment_name(Segment s) {
        switch (s) {
        case Segment::P:   return "P";
        case Segment::QRS: return "QRS";
        default:           return "T";
        }
    }

    enum class Kind : uint8_t { MOVING_AVERAGE, POWER_SPECTRAL, WAVELET };

    inline const char* kind_name(Kind k) {
        switch (k) {
        case Kind::MOVING_AVERAGE: return "moving_average";
        case Kind::POWER_SPECTRAL: return "power_spectral";
        default:                   return "wavelet";
        }
    }

    // ---- Windows (spec 4.7) ------------------------------------------------
    struct WindowConfig {
        int short_beats = 8;    // spec says 5-10; 8 is the midpoint
        int long_beats = 30;   // spec fixes this one
        bool trailing = true; // window ends at the current beat (causal), vs centred
    };

    // Column span of one segment, in the aligned beats' shared sample space.
    struct SegmentBounds {
        int lo = -1, hi = -1;
        bool valid() const { return lo >= 0 && hi > lo; }
    };

    // ---- Output shapes ----------------------------------------------------

    // Per-sample band over a window of beats. Three curves so a min/max band,
    // a mean +/- k*std band, or an IQR band all fit -- pick ONE and document it
    // once the spec is pinned down.
    struct EnvelopeBand {
        std::vector<double> lo, mid, hi;
        int n_beats = 0;                  // beats actually contributing
        bool valid() const { return !mid.empty() && n_beats > 0; }
    };

    // One scalar per beat, tracked over the record -- the shape the spec's
    // "feeds template selection" language implies for the spectral/wavelet
    // kinds (a drift signal, not a spectrum per beat).
    struct EnvelopeTrace {
        std::vector<double> value;         // parallel to the beat series
        std::vector<char>   valid;         // 0 where the window was too short
    };

    // Everything for one segment at one window length.
    struct SegmentEnvelope {
        Segment       segment = Segment::QRS;
        int           window_beats = 0;
        EnvelopeBand  moving_average;      // Kind::MOVING_AVERAGE
        EnvelopeTrace spectral;            // Kind::POWER_SPECTRAL
        EnvelopeTrace wavelet;             // Kind::WAVELET
    };

    // The full 3-segment x 2-window set for one lead.
    struct EnvelopeSet {
        SegmentEnvelope short_p, short_qrs, short_t;
        SegmentEnvelope long_p, long_qrs, long_t;
    };

    // ---- Stubs ------------------------------------------------------------

    // Segment column spans, taken from the reference template's landmarks so
    // every beat in a window is sample-comparable.
    // TODO: P = [p_begin, p_end], QRS = [q_onset, s_end], T = [t_begin, t_end],
    //       all from FeatureMarks against the reference, clamped to the array.
    inline SegmentBounds bounds_for(Segment /*seg*/,
        const std::vector<double>& /*reference*/,
        int /*r_col*/, double /*fs*/)
    {
        return {};
    }

    // Per-sample band across the beats in one window. See the header note on
    // which of min/max, mean+/-k*std or IQR is meant.
    // TODO
    inline EnvelopeBand moving_average_band(
        const std::vector<std::vector<double>>& /*beats*/,
        const std::vector<char>& /*accepted*/,
        int /*t*/, int /*window_beats*/,
        SegmentBounds /*bounds*/,
        const WindowConfig & /*cfg*/ = {})
    {
        return {};
    }

    // Spectral measure for the window ending (or centred) at beat t.
    // TODO: define the output first -- band-power ratio is the most useful
    //       scalar for a drift signal, and is well conditioned even on a
    //       ~100-sample P segment, unlike a full PSD.
    inline double spectral_measure(
        const std::vector<std::vector<double>>& /*beats*/,
        const std::vector<char>& /*accepted*/,
        int /*t*/, int /*window_beats*/,
        SegmentBounds /*bounds*/, double /*fs*/,
        const WindowConfig & /*cfg*/ = {})
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // Wavelet measure for the same window.
    // TODO: needs a basis, a level count and a decomposition choice, none of
    //       which the spec gives, and no wavelet code exists in this repo yet.
    inline double wavelet_measure(
        const std::vector<std::vector<double>>& /*beats*/,
        const std::vector<char>& /*accepted*/,
        int /*t*/, int /*window_beats*/,
        SegmentBounds /*bounds*/, double /*fs*/,
        const WindowConfig & /*cfg*/ = {})
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // ---- Top level --------------------------------------------------------
    // Build the whole set for one lead's aligned beats.
    //
    // `accepted` gates membership: only category-1 beats should enter an
    // envelope, or a PVC widens the band it is supposed to be detected against.
    // The beats stay in the SERIES (so beat index and RR keep their meaning) --
    // they are just skipped as window members.
    // TODO
    inline EnvelopeSet build(const std::vector<std::vector<double>>& /*beats*/,
        const std::vector<char>& /*accepted*/,
        const std::vector<double>& /*reference*/,
        int /*r_col*/, double /*fs*/,
        const WindowConfig& cfg = {})
    {
        EnvelopeSet out;
        out.short_p.segment = Segment::P;    out.short_p.window_beats = cfg.short_beats;
        out.short_qrs.segment = Segment::QRS;  out.short_qrs.window_beats = cfg.short_beats;
        out.short_t.segment = Segment::T;    out.short_t.window_beats = cfg.short_beats;
        out.long_p.segment = Segment::P;    out.long_p.window_beats = cfg.long_beats;
        out.long_qrs.segment = Segment::QRS;  out.long_qrs.window_beats = cfg.long_beats;
        out.long_t.segment = Segment::T;    out.long_t.window_beats = cfg.long_beats;
        return out;
    }

}   // namespace envelopes