#pragma once
/**
 * @file   envelope_output.hpp
 * @brief  Writes the Section 4.7 dynamic envelopes as CSV. The missing half of
 *         Task C's output.
 *
 *         WHY THIS FILE EXISTS. five_category_output.hpp says "Sections 4.5 to
 *         4.7" in its header and emits six CSVs, but every one of them is 4.5
 *         or 4.6: categories, per-bin percentages, the template bank, NSVT
 *         runs, the reference average, the substituted beats. Nothing in it
 *         touches envelopes.hpp. So 4.7 -- moving-average, power-spectral and
 *         wavelet envelopes of the P, QRS and T waves over a short (5 to 10
 *         beat) and a long (30 beat) window -- was being computed by nobody
 *         and written nowhere. envelopes.hpp is a pure compute header: it
 *         returns DynamicEnvelopes and EnvelopeVerdict and opens no file. This
 *         is its writer.
 *
 *         WHAT IT WRITES, per record stem, into cfg.five_category_output:
 *
 *           {stem}_envelope_scalars.csv   the per-beat scalar each envelope
 *                                         tracks (P amplitude, QRS area, T
 *                                         area) plus the baseline and the
 *                                         window each was measured over. This
 *                                         is the audit trail: an envelope is
 *                                         only as good as the series under it.
 *           {stem}_envelopes.csv          one row per beat PER WAVE: the
 *                                         moving-average corridor, slope and
 *                                         z, the three spectral band
 *                                         fractions and the alternans SNR, and
 *                                         the per-level wavelet energies --
 *                                         each at both window lengths.
 *           {stem}_envelope_verdicts.csv  one row per beat per wave: stable /
 *                                         step / drifting / alternating and
 *                                         the six numbers behind them. This is
 *                                         what Section 9.5's master-template
 *                                         selection reads.
 *           {stem}_envelope_summary.csv   one row per bin per wave: the share
 *                                         of beats stable, stepping, drifting
 *                                         and alternating.
 *           {stem}_envelope_config.csv    one row: the windows, the wavelet
 *                                         level count and every verdict
 *                                         threshold actually used. 4.7 gives
 *                                         two window lengths and no
 *                                         thresholds, so the numbers are
 *                                         implementation and the CSV has to
 *                                         say which ones it ran with.
 *
 *         JOINS TO THE 4.5/4.6 OUTPUT. (bin, beat_in_bin) and the running
 *         `beat` column are generated the same way five_category_output.hpp
 *         generates them: bins in order, bad_segment bins contributing no
 *         rows. So {stem}_envelopes.csv joins straight onto
 *         {stem}_beats.csv and no category or template column is duplicated
 *         here.
 *
 *         GEOMETRY IS SHARED, NOT REIMPLEMENTED. SliceGeometry, geometryOf
 *         and medianOf are thin aliases over five_category_output.hpp's rather than being
 *         copied. A beat slice spans 1.8 RR resampled to a fixed width, so
 *         ms-per-column depends on the bin's own rate; if the two writers
 *         computed it separately they could disagree and the join would be
 *         between two different time axes.
 *
 *         ENVELOPES RUN ON EVERY BEAT, NOT THE CLEAN POOL. 4.7 wants the short
 *         window to capture PVC onset. A PVC excluded before the envelope is
 *         built is exactly the beat whose onset was supposed to show up in it,
 *         so nothing is filtered here. Use the join to _beats.csv if you want
 *         the category-1 subset.
 *
 *         THE T WINDOW. Task B's Segments is used when it carries T bounds
 *         (tLo/tHi, detected at compile time). When it does not, the window
 *         falls back to a rate-proportional one: 0.10 to 0.55 RR after R,
 *         floored at the QRS offset. That fallback is an implementation
 *         choice and is recorded as such -- t_window_source in the scalars CSV
 *         says which of the two produced each row, so a T-wave drift reading
 *         is never mistaken for a detected T wave when it was a geometric
 *         guess.
 *
 *         NON-FINITE VALUES ARE WRITTEN AS EMPTY FIELDS, not as "nan" and not
 *         as 0. An unmeasurable P wave is a gap; envelopes.hpp deliberately
 *         carries it through as NaN rather than interpolating it away, and a
 *         zero here would put a fabricated amplitude into whatever reads the
 *         CSV. (Worth backporting to five_category_output.hpp, which streams
 *         raw doubles -- on MSVC that emits -nan(ind), which most CSV readers
 *         will not parse.)
 *
 *         CALL SITE. Beside the existing Task C call in main.cpp, which
 *         already has both inputs in scope:
 *
 *           five_category_output::runAll(job->cfg, job->beats, job->tmpl,
 *                                       job->stem);
 *           envelope_output::runAll   (job->cfg, job->beats, job->tmpl,
 *                                       job->stem);
 *
 *         Same reasoning as the comment above that call applies: it must NOT
 *         go inside the needsFinalize block, or records that already have wave
 *         markings silently produce nothing.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <ostream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "config_file_handling/config_entry.hpp"
#include "template_generation/template_io.hpp"
#include "five_category_classification/envelopes.hpp"
#include "five_category_classification/five_category_output.hpp"
#include "logging/sqi_ecg.hpp"

namespace envelope_output {

    // Borrowed from five_category_output rather than reimplemented, so the
    // two writers cannot disagree about a bin's ms-per-column and their CSVs
    // stay joinable.
    //
    // ALIASED AND QUALIFIED, NOT using-DECLARED. A namespace-scope
    // `using five_category_output::SliceGeometry;` is legal C++ but MSVC
    // rejected the later `SliceGeometry geo = ...` with C2146 (missing ';'
    // before identifier) -- it did not treat the name as a type at that point.
    // A type alias and a qualified call are equivalent and unambiguous.
    using SliceGeometry = five_category_output::SliceGeometry;
    inline constexpr double kNaN = five_category_output::kNaN;

    inline SliceGeometry geometryOf(const template_io::ChannelMethodTemplate& t,
        double medianRrMs)
    {
        return five_category_output::geometryOf(t, medianRrMs);
    }

    inline double medianOf(std::vector<double> v) {
        return five_category_output::medianOf(std::move(v));
    }

    // ---------------------------------------------------------------------
    // configuration
    // ---------------------------------------------------------------------
    // 4.7 specifies the two window lengths and nothing else. Everything below
    // is implementation and is written into {stem}_envelope_config.csv so the
    // output states the parameters it was produced under.
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

    // ---------------------------------------------------------------------
    // CSV field helpers
    // ---------------------------------------------------------------------
    inline void putNum(std::ostream& os, double v) {
        if (std::isfinite(v)) os << v;      // else: empty field
    }
    inline void putInt(std::ostream& os, int v, bool valid = true) {
        if (valid) os << v;
    }

    namespace detail {

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
        inline std::pair<int, int> tWindow(const S& seg, const SliceGeometry& geo,
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
            const Segments& seg, const SliceGeometry& geo,
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
        // parallel to kWaves
        std::array<std::vector<envelopes::EnvelopeVerdict>, 3> verdicts;
    };

    inline BinEnvelopes runBin(const std::vector<std::vector<double>>& beats,
        const std::vector<double>& rrMs,
        const template_io::ChannelMethodTemplate& tmplRaw,
        int binIndex,
        const ScalarConfig& cfg,
        const envelopes::VerdictThresholds& th)
    {
        BinEnvelopes out;
        out.bin = binIndex;
        if (beats.empty()) return out;

        out.nBeats = static_cast<int>(beats.size());
        out.medianRrMs = medianOf(rrMs);

        // Geometry from the bin's own template when there is one, else from
        // the first beat -- the width and R column are all that is needed and
        // the slices in a bin share both.
        SliceGeometry geo = geometryOf(tmplRaw, out.medianRrMs);
        if (geo.width == 0 && !beats[0].empty()) {
            geo.width = static_cast<int>(beats[0].size());
            if (out.medianRrMs > 0.0) {
                geo.msPerCol = 1.8 * out.medianRrMs / geo.width;
                geo.sliceFs = 1000.0 / geo.msPerCol;
            }
        }
        out.msPerCol = geo.msPerCol;

        // ---- the three per-beat series (4.7) ----------------------------
        out.scalars.reserve(beats.size());
        std::vector<envelopes::BeatWaveScalars> series;
        series.reserve(beats.size());
        for (const std::vector<double>& beat : beats) {
            const Segments seg = buildSegments(beat, geo.rCol, geo.sliceFs);
            detail::BeatScalarRow r =
                detail::scalarsOf(beat, seg, geo, out.medianRrMs, cfg);
            series.push_back(r.s);
            out.scalars.push_back(std::move(r));
        }

        // ---- three envelope kinds x two windows x three waves -----------
        out.env = envelopes::buildDynamicEnvelopes(series, cfg.shortWindow,
            cfg.longWindow, cfg.centred);

        // ---- the readouts Section 9.5 consumes --------------------------
        for (std::size_t w = 0; w < kWaves.size(); ++w)
            out.verdicts[w] = envelopes::verdicts(out.env.of(kWaves[w]), th);

        return out;
    }

    // ---------------------------------------------------------------------
    // the entry point
    // ---------------------------------------------------------------------
    //
    // `rrPerBin` is optional and is used only for the rate-proportional
    // geometry (ms-per-column, the fallback T window). The envelopes
    // themselves are indexed by beat, not by time, so they do not depend on
    // it. Pass real intervals when the caller has them.
    inline void runAll(const template_io::BeatsFile& beats,
        const template_io::TemplateFile& tmpl,
        double ecgRate,
        const std::string& dir,
        const std::string& stem,
        const std::vector<std::vector<double>>& rrPerBin = {},
        const ScalarConfig& cfg = {},
        const envelopes::VerdictThresholds& th = {})
    {
        if (dir.empty() || stem.empty()) {
            std::fprintf(stderr, "envelope_output: empty dir or stem; "
                "nothing written\n");
            return;
        }
        {   // An ofstream onto a missing directory fails silently, which looks
            // identical to the feature being switched off.
            std::error_code ec;
            std::filesystem::create_directories(dir, ec);
            if (ec) {
                std::fprintf(stderr, "envelope_output: cannot create %s: %s\n",
                    dir.c_str(), ec.message().c_str());
                return;
            }
        }
        // pickBeatSource locates the channel that actually has beats and
        // returns that lead's templates with them.
        const five_category_output::BeatSource src =
            five_category_output::pickBeatSource(beats);
        if (!src.beats) {
            std::fprintf(stderr, "envelope_output: %s has no beats in "
                "per_channel_beats[CH1/CH2/CH3]; nothing "
                "written\n", stem.c_str());
            return;
        }
        const std::size_t nBins = std::min(src.beats->size(),
            tmpl.bins.size());
        if (nBins == 0) {
            std::fprintf(stderr, "envelope_output: no bins for %s\n",
                stem.c_str());
            return;
        }

        std::vector<BinEnvelopes> results;
        results.reserve(nBins);
        for (std::size_t b = 0; b < nBins; ++b) {
            // Bad bins contribute no rows -- same rule as
            // five_category_output, so the running `beat` index agrees.
            if (b < beats.bad_segment.size() && beats.bad_segment[b]) {
                results.push_back(BinEnvelopes{});
                continue;
            }
            const std::vector<std::vector<double>>& bb = (*src.beats)[b];
            const double sliceRr = (!bb.empty() && !bb[0].empty() && ecgRate > 0.0)
                ? (bb[0].size() / 1.8) * (1000.0 / ecgRate) : 800.0;
            const std::vector<double> rr =
                (b < rrPerBin.size() && rrPerBin[b].size() == bb.size())
                ? rrPerBin[b]
                : std::vector<double>(bb.size(), sliceRr);
            results.push_back(runBin(bb, rr, tmpl.bins[b].*src.raw,
                static_cast<int>(b), cfg, th));
        }

        const std::string base = dir + "/" + stem;
        const int levels = envelopes::kWaveletLevels;
        int nRows = 0, nBeatsTotal = 0, nStep = 0, nAlt = 0, nDrift = 0;

        // ---- 1. the scalar series ---------------------------------------
        {
            std::ofstream f(base + "_envelope_scalars.csv");
            f << std::setprecision(9);
            f << "bin,beat_in_bin,beat,baseline,p_amp,qrs_area,t_area,"
                "ms_per_col,median_rr_ms,p_lo,p_hi,qrs_lo,qrs_hi,t_lo,t_hi,"
                "t_window_source\n";
            int running = 0;
            for (const BinEnvelopes& r : results) {
                for (std::size_t i = 0; i < r.scalars.size(); ++i) {
                    const detail::BeatScalarRow& s = r.scalars[i];
                    f << r.bin << ',' << i << ',' << running++ << ',';
                    putNum(f, s.baseline);   f << ',';
                    putNum(f, s.s.p);        f << ',';
                    putNum(f, s.s.qrs);      f << ',';
                    putNum(f, s.s.t);        f << ',';
                    putNum(f, r.msPerCol);   f << ',';
                    putNum(f, r.medianRrMs); f << ',';
                    putInt(f, s.pLo, s.pLo >= 0);     f << ',';
                    putInt(f, s.pHi, s.pLo >= 0);     f << ',';
                    putInt(f, s.qrsLo, s.qrsLo >= 0); f << ',';
                    putInt(f, s.qrsHi, s.qrsLo >= 0); f << ',';
                    putInt(f, s.tLo, s.tLo >= 0);     f << ',';
                    putInt(f, s.tHi, s.tLo >= 0);     f << ','
                        << detail::tSourceName(s.tSource) << '\n';
                    ++nBeatsTotal;
                }
            }
        }

        // ---- 2. the envelopes themselves --------------------------------
        {
            std::ofstream f(base + "_envelopes.csv");
            f << std::setprecision(9);
            f << "bin,beat_in_bin,beat,wave,scalar,"
                "ma_s_mean,ma_s_sd,ma_s_lo,ma_s_hi,ma_s_z,ma_s_slope,"
                "ma_s_trend,ma_s_n,"
                "ma_l_mean,ma_l_sd,ma_l_lo,ma_l_hi,ma_l_z,ma_l_slope,"
                "ma_l_trend,ma_l_n,"
                "ps_s_total,ps_s_drift_frac,ps_s_mid_frac,ps_s_alt_frac,"
                "ps_s_alt_snr,ps_s_dom_cpb,"
                "ps_l_total,ps_l_drift_frac,ps_l_mid_frac,ps_l_alt_frac,"
                "ps_l_alt_snr,ps_l_dom_cpb";
            for (int l = 1; l <= levels; ++l) f << ",wv_s_l" << l;
            f << ",wv_s_dom_level";
            for (int l = 1; l <= levels; ++l) f << ",wv_l_l" << l;
            f << ",wv_l_dom_level\n";

            int running = 0;
            for (const BinEnvelopes& r : results) {
                for (std::size_t i = 0; i < r.scalars.size(); ++i) {
                    const int beatIdx = running++;
                    for (std::size_t wi = 0; wi < kWaves.size(); ++wi) {
                        const envelopes::Wave w = kWaves[wi];
                        const envelopes::WaveEnvelopes& e = r.env.of(w);
                        f << r.bin << ',' << i << ',' << beatIdx << ','
                            << envelopes::waveName(w) << ',';
                        putNum(f, r.scalars[i].s.get(w)); f << ',';

                        const auto ma = [&](const envelopes::MovingAverageEnvelope& m) {
                            const bool ok = i < m.mean.size();
                            putNum(f, ok ? m.mean[i] : kNaN);  f << ',';
                            putNum(f, ok ? m.sd[i] : kNaN);    f << ',';
                            putNum(f, ok ? m.lo[i] : kNaN);    f << ',';
                            putNum(f, ok ? m.hi[i] : kNaN);    f << ',';
                            putNum(f, ok ? m.z[i] : kNaN);     f << ',';
                            putNum(f, ok ? m.slope[i] : kNaN); f << ',';
                            putNum(f, ok ? m.trend[i] : kNaN); f << ',';
                            f << (ok ? m.nContrib[i] : 0) << ',';
                            };
                        ma(e.maShort);
                        ma(e.maLong);

                        const auto ps = [&](const envelopes::SpectralEnvelope& s,
                            bool trailingComma) {
                                const envelopes::SpectralPoint p =
                                    (i < s.pt.size()) ? s.pt[i] : envelopes::SpectralPoint{};
                                putNum(f, p.total);          f << ',';
                                putNum(f, p.driftFrac);      f << ',';
                                putNum(f, p.midFrac);        f << ',';
                                putNum(f, p.alternansFrac);  f << ',';
                                putNum(f, p.alternansSnr);   f << ',';
                                putNum(f, p.dominantCpb);
                                if (trailingComma) f << ',';
                            };
                        ps(e.psShort, true);
                        ps(e.psLong, false);

                        const auto wv = [&](const envelopes::WaveletEnvelope& v) {
                            for (int l = 0; l < levels; ++l) {
                                f << ',';
                                const bool ok =
                                    l < static_cast<int>(v.energy.size())
                                    && i < v.energy[static_cast<std::size_t>(l)].size();
                                putNum(f, ok
                                    ? v.energy[static_cast<std::size_t>(l)][i] : kNaN);
                            }
                            f << ',';
                            putNum(f, (i < v.dominantLevel.size())
                                ? v.dominantLevel[i] : kNaN);
                            };
                        wv(e.wvShort);
                        wv(e.wvLong);
                        f << '\n';
                        ++nRows;
                    }
                }
            }
        }

        // ---- 3. the verdicts Section 9.5 reads --------------------------
        {
            std::ofstream f(base + "_envelope_verdicts.csv");
            f << std::setprecision(9);
            f << "bin,beat_in_bin,beat,wave,stable,step,drifting,alternating,"
                "z_short,z_long,alt_frac_short,alt_snr_short,drift_frac_long,"
                "trend_long\n";
            int running = 0;
            for (const BinEnvelopes& r : results) {
                for (std::size_t i = 0; i < r.scalars.size(); ++i) {
                    const int beatIdx = running++;
                    for (std::size_t wi = 0; wi < kWaves.size(); ++wi) {
                        if (i >= r.verdicts[wi].size()) continue;
                        const envelopes::EnvelopeVerdict& v = r.verdicts[wi][i];
                        f << r.bin << ',' << i << ',' << beatIdx << ','
                            << envelopes::waveName(kWaves[wi]) << ','
                            << (v.stable ? 1 : 0) << ',' << (v.step ? 1 : 0) << ','
                            << (v.drifting ? 1 : 0) << ','
                            << (v.alternating ? 1 : 0) << ',';
                        putNum(f, v.zShort);    f << ',';
                        putNum(f, v.zLong);     f << ',';
                        putNum(f, v.altFrac);   f << ',';
                        putNum(f, v.altSnr);    f << ',';
                        putNum(f, v.driftFrac); f << ',';
                        putNum(f, v.trendLong); f << '\n';
                        if (v.step)        ++nStep;
                        if (v.alternating) ++nAlt;
                        if (v.drifting)    ++nDrift;
                    }
                }
            }
        }

        // ---- 4. per-bin per-wave summary --------------------------------
        {
            std::ofstream f(base + "_envelope_summary.csv");
            f << std::setprecision(9);
            f << "bin,wave,n_beats,n_scorable,pct_stable,pct_step,"
                "pct_drifting,pct_alternating,median_abs_z_short,"
                "median_trend_long,median_alt_snr_short,short_window,"
                "long_window\n";
            for (const BinEnvelopes& r : results) {
                if (r.nBeats == 0) continue;
                for (std::size_t wi = 0; wi < kWaves.size(); ++wi) {
                    const std::vector<envelopes::EnvelopeVerdict>& v = r.verdicts[wi];
                    // Scorable, not total: an unmeasurable beat is not a stable
                    // one, and dividing by the total would report a record with
                    // no detectable P wave as 100 percent P-unstable.
                    int scorable = 0, stable = 0, step = 0, drift = 0, alt = 0;
                    std::vector<double> zs, tr, sn;
                    for (const envelopes::EnvelopeVerdict& q : v) {
                        if (!std::isfinite(q.zShort)) continue;
                        ++scorable;
                        stable += q.stable ? 1 : 0;
                        step += q.step ? 1 : 0;
                        drift += q.drifting ? 1 : 0;
                        alt += q.alternating ? 1 : 0;
                        zs.push_back(q.zShort);
                        tr.push_back(q.trendLong);
                        sn.push_back(q.altSnr);
                    }
                    const double d = scorable > 0 ? 100.0 / scorable : kNaN;
                    f << r.bin << ',' << envelopes::waveName(kWaves[wi]) << ','
                        << r.nBeats << ',' << scorable << ',';
                    putNum(f, scorable ? stable * d : kNaN); f << ',';
                    putNum(f, scorable ? step * d : kNaN); f << ',';
                    putNum(f, scorable ? drift * d : kNaN); f << ',';
                    putNum(f, scorable ? alt * d : kNaN); f << ',';
                    putNum(f, detail::medianAbs(zs)); f << ',';
                    putNum(f, medianOf(tr));          f << ',';
                    putNum(f, medianOf(sn));          f << ','
                        << cfg.shortWindow << ',' << cfg.longWindow << '\n';
                }
            }
        }

        // ---- 5. the parameters this run used ----------------------------
        {
            std::ofstream f(base + "_envelope_config.csv");
            f << std::setprecision(9);
            f << "short_window,long_window,wavelet_levels,centred,"
                "min_beats_for_spectrum,max_attainable_z_short,z_step,"
                "alt_frac,alt_snr,trend_sd,drift_frac,baseline_lead_frac,"
                "t_start_frac_rr,t_end_frac_rr\n";
            f << cfg.shortWindow << ',' << cfg.longWindow << ',' << levels << ','
                << (cfg.centred ? 1 : 0) << ',' << envelopes::kMinForSpectrum << ',';
            // The z ceiling for the chosen short window. A step threshold above
            // it can never fire, which is why zStep is not a textbook 3-sigma.
            putNum(f, envelopes::maxAttainableZ(cfg.shortWindow)); f << ',';
            putNum(f, th.zStep);     f << ',';
            putNum(f, th.altFrac);   f << ',';
            putNum(f, th.altSnr);    f << ',';
            putNum(f, th.trendSd);   f << ',';
            putNum(f, th.driftFrac); f << ',';
            putNum(f, cfg.baselineLeadFrac); f << ',';
            putNum(f, cfg.tStartFracRr);     f << ',';
            putNum(f, cfg.tEndFracRr);       f << '\n';
        }

        std::fprintf(stderr, "  envelope output: %d beats, %zu bins, %d rows "
            "(%d step, %d alternans, %d drift) -> %s_envelope*.csv\n",
            nBeatsTotal, nBins, nRows, nStep, nAlt, nDrift, base.c_str());
    }

    inline void runAll(const config_entry& cfg,
        const template_io::BeatsFile& beats,
        const template_io::TemplateFile& tmpl,
        const std::string& stem,
        const std::vector<std::vector<double>>& rrPerBin = {},
        const ScalarConfig& scfg = {},
        const envelopes::VerdictThresholds& th = {})
    {
        runAll(beats, tmpl, cfg.ecg_upsample_rate, cfg.five_category_output,
            stem, rrPerBin, scfg, th);
    }

} // namespace envelope_output