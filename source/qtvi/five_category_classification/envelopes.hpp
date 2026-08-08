#pragma once
//
// envelopes.hpp
//
// Spec 4.7: moving-average, power-spectral and wavelet envelopes of the QRS, P
// and T waves over a short (5-10 beat) and a long (30 beat) window. Short
// windows track PVC onset and conduction change; long windows track ischaemic ST
// drift. Both feed the dynamic master template selection of spec 9.5.
//
// ---------------------------------------------------------------------------
// THE SPEC DEFINES NONE OF THE THREE NUMERICALLY, so each is pinned down here.
// The choices are documented at the point they are made -- change them if the
// intent was different, but they are now at least explicit and testable:
//
//   MOVING_AVERAGE -- a per-sample IQR band (Q1 / median / Q3) across the
//     window's beats. Chosen over min/max because a single artifact beat sets
//     min/max forever, and over mean +/- k*sigma because the rest of this
//     codebase aggregates templates by median. Band WIDTH (Q3 - Q1) is the
//     scalar that tracks morphology instability.
//
//   POWER_SPECTRAL -- a scalar per beat, not a spectrum: the fraction of the
//     segment's power above a cutoff, via Goertzel bins. A ~100-sample P wave is
//     far too short for a meaningful PSD, and "feeds template selection" wants a
//     drift signal rather than a per-beat spectrum. Segments are NOT
//     concatenated: each beat's segment is transformed and the ratios averaged
//     over the window, so a single bad beat can't smear across the boundary.
//
//   WAVELET -- normalized Haar detail energy at a chosen level. Haar because it
//     needs no coefficient tables and no dependency, and this repo has no
//     wavelet code at all. Level 1 responds to sample-scale change, level 3 to
//     ~8-sample features; the default of 2 sits between.
//
// Alignment: envelopes across beats only mean anything if the beats are
// column-aligned first. alignment::ecg_beat_set already is, so these consume
// aligned beats, never raw slices.
//
// Segment bounds come from the REFERENCE template's landmarks, i.e. fixed
// columns for the whole window. Per-beat bounds would follow each beat's own
// morphology, but then the window's members are no longer sample-comparable and
// a band across them is meaningless.
// ---------------------------------------------------------------------------

#include <vector>
#include <cstdint>
#include <cmath>
#include <limits>
#include <algorithm>

#include "feature_marks.hpp"

namespace envelopes {

    inline constexpr double kPi = 3.14159265358979323846;

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

    struct WindowConfig {
        int    short_beats = 8;      // spec says 5-10
        int    long_beats = 30;     // spec fixes this
        bool   trailing = true;   // window ends at the current beat (causal)
        double spectral_cutoff_hz = 20.0;   // "high" band for the power ratio
        int    wavelet_level = 2;      // Haar detail level
        int    min_members = 3;      // fewer contributing beats -> invalid
    };

    struct SegmentBounds {
        int lo = -1, hi = -1;
        bool valid() const { return lo >= 0 && hi > lo; }
        int  width() const { return valid() ? hi - lo + 1 : 0; }
    };

    // Per-sample band over a window of beats. lo/mid/hi are Q1/median/Q3.
    struct EnvelopeBand {
        std::vector<double> lo, mid, hi;
        int  n_beats = 0;
        bool valid() const { return !mid.empty() && n_beats > 0; }
        // Mean band width -- the scalar that says "how unstable is this segment".
        double mean_width() const {
            double s = 0.0; int m = 0;
            for (size_t j = 0; j < mid.size(); ++j) {
                if (j < lo.size() && j < hi.size() && !std::isnan(lo[j]) && !std::isnan(hi[j])) {
                    s += hi[j] - lo[j]; ++m;
                }
            }
            return (m > 0) ? s / m : std::numeric_limits<double>::quiet_NaN();
        }
    };

    // One scalar per beat, tracked over the record.
    struct EnvelopeTrace {
        std::vector<double> value;
        std::vector<char>   valid;
        void resize(size_t n) {
            value.assign(n, std::numeric_limits<double>::quiet_NaN());
            valid.assign(n, 0);
        }
    };

    struct SegmentEnvelope {
        Segment       segment = Segment::QRS;
        int           window_beats = 0;
        SegmentBounds bounds;
        EnvelopeBand  band;          // at the last beat of the record
        EnvelopeTrace band_width;    // MOVING_AVERAGE, as a per-beat scalar
        EnvelopeTrace spectral;      // POWER_SPECTRAL
        EnvelopeTrace wavelet;       // WAVELET
    };

    struct EnvelopeSet {
        SegmentEnvelope short_p, short_qrs, short_t;
        SegmentEnvelope long_p, long_qrs, long_t;
    };

    // =======================================================================
    // Primitives
    // =======================================================================

    inline double quantile_of(std::vector<double> v, double q) {
        v.erase(std::remove_if(v.begin(), v.end(),
            [](double x) { return !std::isfinite(x); }), v.end());
        if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
        std::sort(v.begin(), v.end());
        const double pos = q * (v.size() - 1);
        const size_t i = static_cast<size_t>(std::floor(pos));
        const size_t k = std::min(i + 1, v.size() - 1);
        const double f = pos - i;
        return v[i] * (1.0 - f) + v[k] * f;
    }

    // Single-frequency power, one DFT bin, O(n), no library.
    inline double goertzel_power(const std::vector<double>& v, int lo, int hi,
        double freq_hz, double fs)
    {
        lo = std::max(0, lo); hi = std::min(hi, static_cast<int>(v.size()) - 1);
        const int n = hi - lo + 1;
        if (n < 8 || fs <= 0.0 || freq_hz <= 0.0 || freq_hz >= 0.5 * fs) return 0.0;
        const double coeff = 2.0 * std::cos(2.0 * kPi * freq_hz / fs);
        double s1 = 0.0, s2 = 0.0;
        for (int i = lo; i <= hi; ++i) {
            const double x = std::isnan(v[i]) ? 0.0 : v[i];
            const double s0 = x + coeff * s1 - s2;
            s2 = s1; s1 = s0;
        }
        const double p = s1 * s1 + s2 * s2 - coeff * s1 * s2;
        return (p > 0.0) ? p / (static_cast<double>(n) * n) : 0.0;
    }

    // Fraction of the segment's power above `cutoff_hz`. Swept in bins of
    // fs/width, which is the natural resolution of a window this short.
    inline double high_band_fraction(const std::vector<double>& v, int lo, int hi,
        double fs, double cutoff_hz)
    {
        lo = std::max(0, lo); hi = std::min(hi, static_cast<int>(v.size()) - 1);
        const int n = hi - lo + 1;
        if (n < 16 || fs <= 0.0) return std::numeric_limits<double>::quiet_NaN();
        const double df = fs / n;                     // bin spacing
        const double nyq = 0.5 * fs;
        double total = 0.0, high = 0.0;
        for (double f = df; f < nyq; f += df) {
            const double p = goertzel_power(v, lo, hi, f, fs);
            total += p;
            if (f >= cutoff_hz) high += p;
        }
        return (total > 0.0) ? high / total : std::numeric_limits<double>::quiet_NaN();
    }

    // Normalized Haar detail energy at `level`. One level = pairwise differences
    // over sqrt(2); the approximation (pairwise sums) is carried forward and the
    // process repeated. Returned as a fraction of total energy so it is
    // amplitude-independent and comparable across leads.
    inline double haar_detail_fraction(const std::vector<double>& v, int lo, int hi,
        int level)
    {
        lo = std::max(0, lo); hi = std::min(hi, static_cast<int>(v.size()) - 1);
        if (hi - lo < (1 << std::max(1, level))) return std::numeric_limits<double>::quiet_NaN();
        std::vector<double> a;
        a.reserve(hi - lo + 1);
        for (int i = lo; i <= hi; ++i) a.push_back(std::isnan(v[i]) ? 0.0 : v[i]);

        double totalEnergy = 0.0;
        for (double x : a) totalEnergy += x * x;
        if (totalEnergy <= 0.0) return std::numeric_limits<double>::quiet_NaN();

        const double invSqrt2 = 1.0 / std::sqrt(2.0);
        double detailEnergy = 0.0;
        for (int L = 1; L <= std::max(1, level); ++L) {
            const size_t m = a.size() / 2;
            if (m < 2) break;
            std::vector<double> approx(m), detail(m);
            for (size_t k = 0; k < m; ++k) {
                approx[k] = (a[2 * k] + a[2 * k + 1]) * invSqrt2;
                detail[k] = (a[2 * k] - a[2 * k + 1]) * invSqrt2;
            }
            detailEnergy = 0.0;
            for (double d : detail) detailEnergy += d * d;
            a = std::move(approx);
        }
        return detailEnergy / totalEnergy;
    }

    // =======================================================================
    // Segment bounds from the reference's landmarks
    // =======================================================================
    inline SegmentBounds bounds_for(Segment seg, const std::vector<double>& reference,
        int r_col, double fs)
    {
        SegmentBounds b;
        const int n = static_cast<int>(reference.size());
        if (n < 16 || r_col < 0 || r_col >= n || fs <= 0.0) return b;
        auto cl = [&](double d) { return std::clamp(static_cast<int>(std::lround(d)), 0, n - 1); };

        const double pPk = FeatureMarks::detect_p_peak(reference, r_col, fs);
        const double pOn = FeatureMarks::compute_p_begin(reference, fs, r_col, pPk);
        const int    pEnd = FeatureMarks::detect_p_end(reference, r_col, fs, pPk);
        const double qOn = FeatureMarks::compute_q_onset(reference, fs, r_col);
        const double j = FeatureMarks::compute_j_point(reference, fs, r_col);
        const double tOn = FeatureMarks::compute_t_begin(reference, fs, r_col, j);
        const double tEnd = FeatureMarks::compute_t_end(reference, fs, r_col, tOn);

        switch (seg) {
        case Segment::P:
            if (pOn >= 0.0 && pEnd > pOn) { b.lo = cl(pOn); b.hi = cl(pEnd); }
            break;
        case Segment::QRS:
            if (qOn >= 0.0 && j > qOn) { b.lo = cl(qOn); b.hi = cl(j); }
            break;
        case Segment::T:
            if (tOn >= 0.0 && tEnd > tOn) { b.lo = cl(tOn); b.hi = cl(tEnd); }
            break;
        }
        if (b.lo >= 0 && b.hi <= b.lo) b = SegmentBounds{};
        return b;
    }

    // =======================================================================
    // Window membership
    // =======================================================================
    // Only accepted (category-1) beats contribute, or a PVC widens the very band
    // it is meant to be detected against. Rejected beats stay in the SERIES, so
    // beat index and RR keep their meaning -- they are just skipped as members.
    inline std::vector<int> window_members(const std::vector<char>& accepted,
        int t, int window_beats, const WindowConfig& cfg)
    {
        std::vector<int> out;
        const int n = static_cast<int>(accepted.size());
        if (t < 0 || t >= n || window_beats < 1) return out;
        int lo, hi;
        if (cfg.trailing) { hi = t; lo = t - window_beats + 1; }
        else { lo = t - window_beats / 2; hi = lo + window_beats - 1; }
        lo = std::max(0, lo); hi = std::min(n - 1, hi);
        for (int i = lo; i <= hi; ++i) if (accepted[i]) out.push_back(i);
        return out;
    }

    // =======================================================================
    // The three envelope kinds
    // =======================================================================

    // Per-sample Q1 / median / Q3 across the window's accepted beats.
    inline EnvelopeBand moving_average_band(
        const std::vector<std::vector<double>>& beats,
        const std::vector<char>& accepted,
        int t, int window_beats,
        SegmentBounds bounds,
        const WindowConfig& cfg = {})
    {
        EnvelopeBand out;
        if (!bounds.valid()) return out;
        const std::vector<int> mem = window_members(accepted, t, window_beats, cfg);
        if (static_cast<int>(mem.size()) < cfg.min_members) return out;

        const int w = bounds.width();
        out.lo.assign(w, std::numeric_limits<double>::quiet_NaN());
        out.mid.assign(w, std::numeric_limits<double>::quiet_NaN());
        out.hi.assign(w, std::numeric_limits<double>::quiet_NaN());
        out.n_beats = static_cast<int>(mem.size());

        std::vector<double> col; col.reserve(mem.size());
        for (int k = 0; k < w; ++k) {
            const int j = bounds.lo + k;
            col.clear();
            for (int i : mem)
                if (j < static_cast<int>(beats[i].size()) && !std::isnan(beats[i][j]))
                    col.push_back(beats[i][j]);
            if (col.size() < static_cast<size_t>(cfg.min_members)) continue;
            out.lo[k] = quantile_of(col, 0.25);
            out.mid[k] = quantile_of(col, 0.50);
            out.hi[k] = quantile_of(col, 0.75);
        }
        return out;
    }

    // High-band power fraction, averaged over the window's members. Each beat's
    // segment is transformed separately -- NOT concatenated -- so a single bad
    // beat cannot smear across a boundary.
    inline double spectral_measure(
        const std::vector<std::vector<double>>& beats,
        const std::vector<char>& accepted,
        int t, int window_beats,
        SegmentBounds bounds, double fs,
        const WindowConfig& cfg = {})
    {
        if (!bounds.valid() || fs <= 0.0) return std::numeric_limits<double>::quiet_NaN();
        const std::vector<int> mem = window_members(accepted, t, window_beats, cfg);
        if (static_cast<int>(mem.size()) < cfg.min_members)
            return std::numeric_limits<double>::quiet_NaN();
        double s = 0.0; int m = 0;
        for (int i : mem) {
            const double r = high_band_fraction(beats[i], bounds.lo, bounds.hi,
                fs, cfg.spectral_cutoff_hz);
            if (std::isfinite(r)) { s += r; ++m; }
        }
        return (m >= cfg.min_members) ? s / m : std::numeric_limits<double>::quiet_NaN();
    }

    // Haar detail-energy fraction, averaged over the window's members.
    inline double wavelet_measure(
        const std::vector<std::vector<double>>& beats,
        const std::vector<char>& accepted,
        int t, int window_beats,
        SegmentBounds bounds, double /*fs*/,
        const WindowConfig& cfg = {})
    {
        if (!bounds.valid()) return std::numeric_limits<double>::quiet_NaN();
        const std::vector<int> mem = window_members(accepted, t, window_beats, cfg);
        if (static_cast<int>(mem.size()) < cfg.min_members)
            return std::numeric_limits<double>::quiet_NaN();
        double s = 0.0; int m = 0;
        for (int i : mem) {
            const double r = haar_detail_fraction(beats[i], bounds.lo, bounds.hi,
                cfg.wavelet_level);
            if (std::isfinite(r)) { s += r; ++m; }
        }
        return (m >= cfg.min_members) ? s / m : std::numeric_limits<double>::quiet_NaN();
    }

    // =======================================================================
    // Top level
    // =======================================================================
    inline SegmentEnvelope build_one(const std::vector<std::vector<double>>& beats,
        const std::vector<char>& accepted,
        Segment seg, int window_beats,
        SegmentBounds bounds, double fs,
        const WindowConfig& cfg)
    {
        SegmentEnvelope e;
        e.segment = seg;
        e.window_beats = window_beats;
        e.bounds = bounds;
        const int n = static_cast<int>(beats.size());
        e.band_width.resize(n);
        e.spectral.resize(n);
        e.wavelet.resize(n);
        if (!bounds.valid()) return e;

        for (int t = 0; t < n; ++t) {
            const EnvelopeBand b = moving_average_band(beats, accepted, t, window_beats, bounds, cfg);
            if (b.valid()) {
                const double w = b.mean_width();
                if (std::isfinite(w)) { e.band_width.value[t] = w; e.band_width.valid[t] = 1; }
            }
            const double sp = spectral_measure(beats, accepted, t, window_beats, bounds, fs, cfg);
            if (std::isfinite(sp)) { e.spectral.value[t] = sp; e.spectral.valid[t] = 1; }
            const double wv = wavelet_measure(beats, accepted, t, window_beats, bounds, fs, cfg);
            if (std::isfinite(wv)) { e.wavelet.value[t] = wv; e.wavelet.valid[t] = 1; }
        }
        // Keep the band itself at the last beat, for inspection.
        if (n > 0) e.band = moving_average_band(beats, accepted, n - 1, window_beats, bounds, cfg);
        return e;
    }

    inline EnvelopeSet build(const std::vector<std::vector<double>>& beats,
        const std::vector<char>& accepted,
        const std::vector<double>& reference,
        int r_col, double fs,
        const WindowConfig& cfg = {})
    {
        const SegmentBounds bP = bounds_for(Segment::P, reference, r_col, fs);
        const SegmentBounds bQ = bounds_for(Segment::QRS, reference, r_col, fs);
        const SegmentBounds bT = bounds_for(Segment::T, reference, r_col, fs);

        EnvelopeSet s;
        s.short_p = build_one(beats, accepted, Segment::P, cfg.short_beats, bP, fs, cfg);
        s.short_qrs = build_one(beats, accepted, Segment::QRS, cfg.short_beats, bQ, fs, cfg);
        s.short_t = build_one(beats, accepted, Segment::T, cfg.short_beats, bT, fs, cfg);
        s.long_p = build_one(beats, accepted, Segment::P, cfg.long_beats, bP, fs, cfg);
        s.long_qrs = build_one(beats, accepted, Segment::QRS, cfg.long_beats, bQ, fs, cfg);
        s.long_t = build_one(beats, accepted, Segment::T, cfg.long_beats, bT, fs, cfg);
        return s;
    }

}   // namespace envelopes