#pragma once
/**
 * @file   envelope_compute.hpp
 * @brief  Section 4.7 dynamic envelopes, as a COMPUTATION. No file I/O.
 *
 *         4.7 asks for moving-average, power-spectral and wavelet envelopes of
 *         the P, QRS and T waves over a short (5 to 10 beat) and a long (30
 *         beat) window, feeding the dynamic master-template selection of 9.5.
 *         envelopes.hpp does the envelope mathematics; this header supplies
 *         the one thing it cannot: the per-beat scalar for each wave. Without
 *         that, envelopes.hpp has nothing to run on -- BeatWaveScalars is
 *         three doubles the caller must fill, and 4.7 never says what they
 *         are.
 *
 *         WHY THIS IS A SEPARATE HEADER. The computation used to live inside
 *         envelope_output.hpp next to the CSV writers, and that file includes
 *         five_category_output.hpp for its geometry helpers. Having
 *         five_category_output call the computation would have closed a
 *         circular include. So the computation moved here, takes its geometry
 *         as four plain numbers, and depends on nothing but envelopes.hpp and
 *         Task B's Segments. five_category_output includes this and runs it
 *         per bin; nothing is written anywhere.
 *
 *         SEGMENTS ARE PASSED IN, NOT RECOMPUTED. buildSegments runs the wave
 *         detectors over the whole slice and the caller has already done it
 *         once per beat. Taking the vector as an argument removes what used to
 *         be a third pass over every beat in the record.
 *
 *         THE SCALARS, and why these three:
 *           P   signed peak deviation from baseline -- amplitude is what drops
 *               out on a PAC or in AF, and the sign catches a retrograde P.
 *           QRS ABSOLUTE area. The complex is biphasic, so a signed area
 *               cancels and a broad PVC can integrate to the same number as a
 *               narrow sinus beat.
 *           T   SIGNED area. Ischaemic ST/T drift is a slow signed shift,
 *               which is what the long-window slope is there to catch; the
 *               modulus would fold an inverting T back onto an upright one.
 *
 *         WHAT CONSUMES IT. EnvelopeVerdict per beat per wave -- stable, step,
 *         drifting, alternating -- is what 9.5 is specified to read. 9.5 is
 *         not implemented, so the verdicts are computed and held on the bin
 *         result and nothing reads them yet. That is the spec's ordering, not
 *         an oversight: 4.7 produces, 9.5 consumes.
 *
 * @date   2026-08-26
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include "five_category_classification/envelopes.hpp"
#include "logging/sqi_ecg.hpp"          // Segments

namespace envelope_compute {

    inline constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

    // A beat slice spans 1.8 RR resampled to a fixed width, so a column is not
    // a fixed number of milliseconds -- it depends on the bin's own rate.
    // Passed in as plain numbers so this header needs no other module.
    struct Geometry {
        int    width = 0;
        int    rCol = -1;
        double msPerCol = kNaN;
        double sliceFs = kNaN;
    };

    struct ScalarConfig {
        // Isoelectric baseline: the median of the slice ahead of the P wave,
        // or of this leading fraction when no P segment was found. The slice
        // starts 0.3 RR before R, so its opening columns are pre-atrial.
        double baselineLeadFrac = 0.10;
        // Fallback T window, as a fraction of the bin's median RR after R.
        double tStartFracRr = 0.10;
        double tEndFracRr = 0.55;
        // Trailing windows by default -- a centred window cannot be evaluated
        // until half of it is in the future. Pass true for the offline variant.
        bool   centred = false;
        int    shortWindow = envelopes::kShortWindow;   // 8, midpoint of 5-10
        int    longWindow = envelopes::kLongWindow;    // 30
    };


    namespace detail {

        // Local, so this header needs nothing from another module.
        inline double medianOf(std::vector<double> v) {
            v.erase(std::remove_if(v.begin(), v.end(),
                [](double x) { return !std::isfinite(x); }), v.end());
            if (v.empty()) return kNaN;
            std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
            return v[v.size() / 2];
        }

        // Does Task B's Segments carry T bounds? Detected rather than assumed:
        // this header must compile against the Segments that exists, and a
        // fabricated T window is reported as a fallback rather than passed off
        // as a detection.
        template <class S, class = void> struct HasTBounds : std::false_type {};
        template <class S>
        struct HasTBounds<S, std::void_t<decltype(std::declval<const S&>().tLo),
            decltype(std::declval<const S&>().tHi)>>
            : std::true_type {};

        enum class TSource { NONE = 0, DETECTED = 1, GEOMETRIC = 2 };

        template <class S>
        inline std::pair<int, int> tWindow(const S& seg, const Geometry& geo,
            double medRrMs, const ScalarConfig& cfg,
            TSource& src)
        {
            src = TSource::NONE;
            int lo = -1, hi = -1;
            if constexpr (HasTBounds<S>::value) {
                const int a = static_cast<int>(seg.tLo);
                const int b = static_cast<int>(seg.tHi);
                if (a >= 0 && b > a) { lo = a; hi = b; src = TSource::DETECTED; }
            }
            if (lo < 0) {
                if (geo.rCol < 0 || !std::isfinite(geo.msPerCol)
                    || geo.msPerCol <= 0.0 || !(medRrMs > 0.0))
                    return { -1, -1 };
                lo = geo.rCol + static_cast<int>(
                    std::lround(cfg.tStartFracRr * medRrMs / geo.msPerCol));
                hi = geo.rCol + static_cast<int>(
                    std::lround(cfg.tEndFracRr * medRrMs / geo.msPerCol));
                // Never let the fallback window swallow the QRS tail.
                if (seg.qrsHi > lo) lo = seg.qrsHi;
                src = TSource::GEOMETRIC;
            }
            lo = std::max(0, lo);
            hi = std::min(geo.width, hi);
            if (hi - lo < 2) { src = TSource::NONE; return { -1, -1 }; }
            return { lo, hi };
        }

        inline double baselineOf(const std::vector<double>& beat, int pLo,
            const ScalarConfig& cfg)
        {
            const int n = static_cast<int>(beat.size());
            if (n == 0) return kNaN;
            int hi = (pLo > 2) ? pLo
                : static_cast<int>(std::lround(cfg.baselineLeadFrac * n));
            hi = std::clamp(hi, 1, n);
            return medianOf(std::vector<double>(beat.begin(), beat.begin() + hi));
        }

        // One beat's three scalars, plus what they were measured over.
        struct BeatScalarRow {
            envelopes::BeatWaveScalars s;
            double  baseline = kNaN;
            int     pLo = -1, pHi = -1, qrsLo = -1, qrsHi = -1, tLo = -1, tHi = -1;
            TSource tSource = TSource::NONE;
        };

        // The scalar per wave. The envelope machinery is indifferent to which
        // scalar it is; these three are chosen for what 4.7 asks the envelopes
        // to detect:
        //   P   -- signed peak deviation from baseline. Amplitude is what
        //          drops out on a PAC or in AF, and the sign catches a
        //          retrograde P.
        //   QRS -- ABSOLUTE area. The complex is biphasic, so a signed area
        //          cancels and a broad PVC can integrate to the same number as
        //          a narrow sinus beat. Absolute area separates them.
        //   T   -- SIGNED area. Ischaemic ST/T drift is a slow signed shift,
        //          which is exactly what the long-window slope is meant to
        //          pick up; taking the modulus here would fold an inverting T
        //          wave back onto an upright one.
        // Areas are in amplitude-units x ms, so they are comparable across
        // bins whose rates differ.
        inline BeatScalarRow scalarsOf(const std::vector<double>& beat,
            const Segments& seg, const Geometry& geo,
            double medRrMs, const ScalarConfig& cfg)
        {
            BeatScalarRow d;
            if (beat.empty()) return d;
            const int n = static_cast<int>(beat.size());

            d.baseline = baselineOf(beat, seg.pLo, cfg);
            const double b = std::isfinite(d.baseline) ? d.baseline : 0.0;
            const double msPerCol = geo.msPerCol;

            // ---- P: signed peak deviation --------------------------------
            if (seg.pLo >= 0 && seg.pHi > seg.pLo && seg.pHi <= n) {
                d.pLo = seg.pLo; d.pHi = seg.pHi;
                double best = 0.0; bool any = false;
                for (int c = seg.pLo; c < seg.pHi; ++c) {
                    if (!std::isfinite(beat[static_cast<std::size_t>(c)])) continue;
                    const double v = beat[static_cast<std::size_t>(c)] - b;
                    if (!any || std::fabs(v) > std::fabs(best)) { best = v; any = true; }
                }
                if (any) d.s.p = best;
            }

            // ---- QRS: absolute area --------------------------------------
            if (seg.qrsLo >= 0 && seg.qrsHi > seg.qrsLo && seg.qrsHi <= n
                && std::isfinite(msPerCol))
            {
                d.qrsLo = seg.qrsLo; d.qrsHi = seg.qrsHi;
                double sum = 0.0; int m = 0;
                for (int c = seg.qrsLo; c < seg.qrsHi; ++c)
                    if (std::isfinite(beat[static_cast<std::size_t>(c)])) {
                        sum += std::fabs(beat[static_cast<std::size_t>(c)] - b);
                        ++m;
                    }
                if (m > 0) d.s.qrs = sum * msPerCol;
            }

            // ---- T: signed area ------------------------------------------
            TSource src = TSource::NONE;
            const std::pair<int, int> tw = tWindow(seg, geo, medRrMs, cfg, src);
            if (tw.first >= 0 && std::isfinite(msPerCol)) {
                d.tLo = tw.first; d.tHi = tw.second; d.tSource = src;
                double sum = 0.0; int m = 0;
                for (int c = tw.first; c < tw.second; ++c)
                    if (std::isfinite(beat[static_cast<std::size_t>(c)])) {
                        sum += beat[static_cast<std::size_t>(c)] - b;
                        ++m;
                    }
                if (m > 0) d.s.t = sum * msPerCol;
            }
            return d;
        }

        inline const char* tSourceName(TSource s) {
            switch (s) {
            case TSource::DETECTED:  return "detected";
            case TSource::GEOMETRIC: return "geometric";
            default:                 return "none";
            }
        }

        inline double medianAbs(std::vector<double> v) {
            for (double& x : v) x = std::fabs(x);
            return medianOf(std::move(v));
        }

    } // namespace detail

    // ---------------------------------------------------------------------
    // per-bin result
    // ---------------------------------------------------------------------
    inline constexpr std::array<envelopes::Wave, 3> kWaves = {
        envelopes::Wave::P, envelopes::Wave::QRS, envelopes::Wave::T
    };

    struct BinEnvelopes {
        int    bin = 0;
        int    nBeats = 0;
        double medianRrMs = kNaN;
        double msPerCol = kNaN;
        std::vector<detail::BeatScalarRow>            scalars;
        envelopes::DynamicEnvelopes                   env;
        // Parallel to kWaves: the readouts 9.5 is specified to consume.
        std::array<std::vector<envelopes::EnvelopeVerdict>, 3> verdicts;

        // Counts, so a caller can see the envelopes ran without walking them.
        int nStep = 0, nDrifting = 0, nAlternating = 0, nStable = 0;
    };

    // `segOf` must be parallel to `beats` -- the caller's already-computed
    // boundaries. Pass an empty vector to run without P/T windows (the QRS
    // scalar still works from geo.rCol).
    inline BinEnvelopes runBin(const std::vector<std::vector<double>>& beats,
        const std::vector<Segments>& segOf,
        const Geometry& geo,
        double medianRrMs,
        int binIndex,
        const ScalarConfig& cfg = {},
        const envelopes::VerdictThresholds& th = {})
    {
        BinEnvelopes out;
        out.bin = binIndex;
        if (beats.empty()) return out;

        out.nBeats = static_cast<int>(beats.size());
        out.medianRrMs = medianRrMs;
        out.msPerCol = geo.msPerCol;

        // ---- the three per-beat series (4.7) ----------------------------
        static const Segments kNoSeg{};
        out.scalars.reserve(beats.size());
        std::vector<envelopes::BeatWaveScalars> series;
        series.reserve(beats.size());
        for (std::size_t i = 0; i < beats.size(); ++i) {
            const Segments& seg = (i < segOf.size()) ? segOf[i] : kNoSeg;
            detail::BeatScalarRow r =
                detail::scalarsOf(beats[i], seg, geo, medianRrMs, cfg);
            series.push_back(r.s);
            out.scalars.push_back(std::move(r));
        }

        // ---- three envelope kinds x two windows x three waves -----------
        // moving-average, power-spectral and wavelet, short (8, the midpoint
        // of the specified 5-to-10) and long (30).
        out.env = envelopes::buildDynamicEnvelopes(series, cfg.shortWindow,
            cfg.longWindow, cfg.centred);

        // ---- the readouts Section 9.5 consumes --------------------------
        for (std::size_t w = 0; w < kWaves.size(); ++w) {
            out.verdicts[w] = envelopes::verdicts(out.env.of(kWaves[w]), th);
            for (const envelopes::EnvelopeVerdict& v : out.verdicts[w]) {
                if (v.step)        ++out.nStep;
                if (v.drifting)    ++out.nDrifting;
                if (v.alternating) ++out.nAlternating;
                if (v.stable)      ++out.nStable;
            }
        }
        return out;
    }

} // namespace envelope_compute