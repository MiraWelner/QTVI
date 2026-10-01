/**
 * @file   channel_offset.hpp
 * @brief  One-time per-file measurement of the fixed hardware delay between
 *         the ECG acquisition path and EACH acquisition path it's paired
 *         against -- PPG on its own, and ABP/ART/ART_PULM as one shared
 *         arterial group (ABP is the representative channel measured; the
 *         resulting lag is applied to all three, since they share one
 *         acquisition path and hardware delay). PPG and the arterial group
 *         are measured and applied INDEPENDENTLY -- their lags are not
 *         assumed equal.
 *
 *         Cross-correlates the R-peak event train against a foot-event
 *         train (PPG valleys, or ABP's self-detected pulse locations) at
 *         1 ms steps over a lag range of 0 to 500 ms. The lag at maximum
 *         correlation is the delay. Runs once per file during the initial
 *         load, before any template slicing.
 *
 *         MEASURED, NEVER APPLIED TO THE SIGNALS. This used to shift the pulse
 *         and arterial signals earlier by the measured lag before anything was
 *         sliced, which silently subtracted a hardware delay from every
 *         R-to-pulse time downstream -- a transit time measured on a shifted
 *         signal is not a transit time. The signals now stay as recorded. A
 *         confident lag is used in exactly two places, both of which change
 *         where we LOOK and neither of which changes a sample:
 *
 *           the pulse window   each beat's systolic-peak search runs over
 *                              [R + lag, next R + lag] instead of [R, next R]
 *                              (alignment::pulseLandmarks), so a long delay
 *                              cannot hand a beat the previous beat's pulse;
 *           the display        the viewer draws the pulse panels shifted
 *                              earlier by the lag, so review sees physiology
 *                              rather than wiring.
 *
 *         The lag, the method and the confidence are stored per record
 *         (RecordLag) as a trailer on <stem>_templates.bin and on
 *         <stem>_template_markings.bin, so every reported time can be
 *         corrected by whoever reads it.
 *
 *         Includes NO project header: bin members are reached through a
 *         template parameter. This matters in this tree because
 *         template_structs.hpp declares its members with an unqualified
 *         `vector` supplied by a forced include / PCH, so pulling it in
 *         earlier than it used to appear makes struct members fail to
 *         declare and surfaces as errors in unrelated files.
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

namespace channel_offset {

    constexpr int kMaxLagMs = 500;   // search 0..500 ms, per spec
    constexpr int kTolMs = 25;    // a foot counts as aligned to an R peak
    // if it lands within this many ms of it
    constexpr double kMinPeakRatio = 2.0;   // best lag must beat the runner-up
    // by this much, else ambiguous
    constexpr double first_n_secs_to_measure_on = 3600.0 * 8;

    struct Result {
        int    lag_ms = 0;
        double ratio = 0.0;       // best correlation / runner-up
        bool   ambiguous = true;  // true => not used for the window or display; review
        size_t n_r_peaks = 0;
        size_t n_feet = 0;
        uint64_t best_score = 0;
        double analyzed_sec = 0.0;
    };

    // Log destination. Same dir/stem pattern as ecg_move_log and
    // premark_beats. Empty dir/stem => no log.
    inline std::string g_dir, g_stem;
    inline void set(const std::string& dir, const std::string& stem) {
        g_dir = dir; g_stem = stem;
    }

    // -------------------------------------------------------------------
    // Slide a foot train against the R-peak train, one ms at a time, and
    // keep the lag that lines them up best.
    //
    // `feetPerBin[i]` are the foot-event sample indices (at `footRate`) for
    // bins[i] -- positionally aligned with `bins`, NOT read off a bin field
    // by name. This is what lets the SAME correlation logic serve both
    // groups: the caller builds feetPerBin from bin.ppgMinAmps for the PPG
    // measurement, or from a fresh detection on the ABP signal for the
    // arterial group (ABP has no precomputed foot-index field the way PPG
    // does -- see create_arterial_templates.hpp's
    // pulse_matched_filter::derivativePulseLocations for the detector used
    // to build that on the fly).
    //
    // Bins must expose: bad_segment, ch1.raw, ecgSignal.
    // Templated so the bin type need not be complete in this header.
    // -------------------------------------------------------------------
    template <class Bins>
    inline Result measure(const Bins& bins, double ecgRate, double footRate,
        const std::vector<std::vector<std::size_t>>& feetPerBin)
    {
        Result out;
        if (ecgRate <= 0.0 || footRate <= 0.0) return out;

        // Both trains as millisecond timestamps on one shared clock.
        std::vector<int> rMs, fMs;

        std::size_t i = 0;
        for (const auto& bin : bins) {
            const std::vector<std::size_t> emptyFeet;
            const std::vector<std::size_t>& feet =
                (i < feetPerBin.size()) ? feetPerBin[i] : emptyFeet;
            ++i;

            const double budget = first_n_secs_to_measure_on - out.analyzed_sec;
            if (budget <= 0.0) break;
            if (bin.bad_segment) continue;
            if (bin.ch1.raw.size() < 2 || feet.size() < 2) continue;

            const double binSec = static_cast<double>(bin.ecgSignal.size()) / ecgRate;
            const double useSec = std::min(binSec, budget);
            if (useSec <= 0.0) continue;

            // Offset each bin onto the shared clock so bins don't overlap.
            const double t0 = out.analyzed_sec;
            out.analyzed_sec += useSec;

            for (size_t r : bin.ch1.raw) {
                const double t = static_cast<double>(r) / ecgRate;
                if (t >= useSec) break;                 // ascending; prefix only
                rMs.push_back(static_cast<int>(std::lround((t0 + t) * 1000.0)));
            }
            for (size_t f : feet) {
                const double t = static_cast<double>(f) / footRate;
                if (t >= useSec) break;
                fMs.push_back(static_cast<int>(std::lround((t0 + t) * 1000.0)));
            }
        }

        out.n_r_peaks = rMs.size();
        out.n_feet = fMs.size();
        if (rMs.size() < 2 || fMs.size() < 2) return out;

        // Membership stamp: rHit[t] is true if an R peak sits within kTolMs
        // of millisecond t. Built once, so scoring a lag is one lookup per
        // foot instead of a search.
        const int span = rMs.back() + kTolMs + 1;
        std::vector<uint8_t> rHit(static_cast<size_t>(span), 0);
        for (int t : rMs) {
            const int lo = std::max(0, t - kTolMs);
            const int hi = std::min(span - 1, t + kTolMs);
            for (int i = lo; i <= hi; ++i) rHit[static_cast<size_t>(i)] = 1;
        }

        // Correlation at every lag: shift the whole foot train back by lag,
        // count how many feet then land on an R peak.
        std::vector<uint64_t> score(static_cast<size_t>(kMaxLagMs) + 1, 0);
        for (int lag = 0; lag <= kMaxLagMs; ++lag) {
            uint64_t s = 0;
            for (int f : fMs) {
                const int t = f - lag;
                if (t >= 0 && t < span && rHit[static_cast<size_t>(t)]) ++s;
            }
            score[static_cast<size_t>(lag)] = s;
        }

        // Best lag, and the best rival outside its immediate neighbourhood
        // (the peak spreads over ~kTolMs, which is not a competing answer).
        int best = 0;
        for (int lag = 1; lag <= kMaxLagMs; ++lag)
            if (score[static_cast<size_t>(lag)] > score[static_cast<size_t>(best)]) best = lag;

        uint64_t rival = 0;
        for (int lag = 0; lag <= kMaxLagMs; ++lag)
            if (std::abs(lag - best) > 2 * kTolMs)
                rival = std::max(rival, score[static_cast<size_t>(lag)]);

        out.lag_ms = best;
        out.best_score = score[static_cast<size_t>(best)];
        if (out.best_score == 0) return out;
        out.ratio = (rival > 0) ? static_cast<double>(out.best_score) / static_cast<double>(rival)
            : static_cast<double>(out.best_score);
        out.ambiguous = (out.ratio < kMinPeakRatio);
        return out;
    }

    // -------------------------------------------------------------------
    // Build a truncated COPY of the annealed segments covering only the
    // first kMaxAnalysisSec, for the probe pass that measures the lag.
    // Without this, measuring would cost a full extra pass over the whole
    // recording; with it the probe costs 5 minutes of processing.
    //
    // Segments must expose: ppg_signal, ecg_signal_1, ecg_signal_2,
    // ecg_signal_3.
    // -------------------------------------------------------------------
    template <class Segments>
    inline Segments make_probe(const Segments& segs, double ecgRate, double ppgRate)
    {
        Segments probe;
        if (ecgRate <= 0.0 || ppgRate <= 0.0) return probe;
        double covered = 0.0;

        auto cut = [](std::vector<double>& sig, double sec, double rate) {
            const size_t keep = static_cast<size_t>(std::llround(sec * rate));
            if (sig.size() > keep) sig.resize(keep);
            };

        for (const auto& seg : segs) {
            const double budget = first_n_secs_to_measure_on - covered;
            if (budget <= 0.0) break;
            const double segSec = static_cast<double>(seg.ecg_signal_1.size()) / ecgRate;
            const double useSec = std::min(segSec, budget);
            if (useSec <= 0.0) continue;

            probe.push_back(seg);                 // copy, then truncate the copy
            auto& p = probe.back();
            cut(p.ppg_signal, useSec, ppgRate);
            cut(p.ecg_signal_1, useSec, ecgRate);
            cut(p.ecg_signal_2, useSec, ecgRate);
            cut(p.ecg_signal_3, useSec, ecgRate);
            covered += useSec;
        }
        return probe;
    }

    // =======================================================================
    // THE RECORD'S LAGS, AS STORED
    // =======================================================================
    //
    // `apply`, which shifted the signals, is gone -- see the header. What a
    // measurement produces now is this: per group, the lag, how it was found,
    // and how sure the method was, written beside the record's templates and
    // markings so the number travels with every output it bears on.
    enum class LagMethod : uint8_t {
        NONE = 0,                    // not measured (not a CHAOS record, or no data)
        RPEAK_FOOT_TRAIN_XCORR = 1,  // measure(): R train vs foot train, 1 ms steps
    };
    inline const char* methodName(uint8_t m) {
        switch (static_cast<LagMethod>(m)) {
        case LagMethod::RPEAK_FOOT_TRAIN_XCORR: return "rpeak_foot_train_xcorr";
        default:                                return "none";
        }
    }

    struct GroupLag {
        int32_t  lag_ms = 0;
        // Confidence: the best lag's score over the best rival's (Result::
        // ratio). confident = it cleared kMinPeakRatio.
        double   confidence = 0.0;
        uint8_t  method = static_cast<uint8_t>(LagMethod::NONE);
        uint8_t  confident = 0;
        uint32_t n_r_peaks = 0;
        uint32_t n_feet = 0;
        double   analyzed_sec = 0.0;

        bool measured() const { return method != static_cast<uint8_t>(LagMethod::NONE); }
        /// The lag the pulse window and the display may use: the measured
        /// one when the method was sure of it, else 0 (look where R says).
        double usableLagMs() const {
            return (measured() && confident) ? static_cast<double>(lag_ms) : 0.0;
        }
    };

    /// PPG on its own; ABP / ART / ART_PULM as one arterial group.
    struct RecordLag {
        GroupLag ppg;
        GroupLag arterial;
    };

    inline GroupLag fromResult(const Result& r, bool measured) {
        GroupLag g;
        if (!measured) return g;
        g.lag_ms = r.lag_ms;
        g.confidence = r.ratio;
        g.method = static_cast<uint8_t>(LagMethod::RPEAK_FOOT_TRAIN_XCORR);
        g.confident = r.ambiguous ? 0 : 1;
        g.n_r_peaks = static_cast<uint32_t>(r.n_r_peaks);
        g.n_feet = static_cast<uint32_t>(r.n_feet);
        g.analyzed_sec = r.analyzed_sec;
        return g;
    }

    // ---- THE TRAILER --------------------------------------------------------
    //
    // APPENDED TO THE END of a file whose own format is untouched, and found
    // from the end: [payload][uint32 payload bytes][8-byte magic]. Both C++
    // readers (templates_io::readBin, readTemplateMarkingsBin) stop after the
    // blocks they know and never look past them, so a file with a trailer reads
    // exactly as before and a file without one simply has no lag. A reader that
    // insists on EOF after its last block (read_morphology_bin.py, if it does)
    // needs to stop at the trailer instead.
    //
    // Fields written one by one, little-endian as the rest of these files are,
    // so the layout is the code below and not a compiler's padding.
    inline constexpr char kTrailerMagic[8] = { 'C','H','L','A','G','R','C','1' };

    namespace detail {
        inline void putGroup(std::string& b, const GroupLag& g) {
            auto put = [&b](const void* p, std::size_t n) {
                b.append(static_cast<const char*>(p), n); };
            put(&g.lag_ms, 4); put(&g.confidence, 8);
            put(&g.method, 1); put(&g.confident, 1);
            put(&g.n_r_peaks, 4); put(&g.n_feet, 4); put(&g.analyzed_sec, 8);
        }
        inline bool getGroup(const std::string& b, std::size_t& at, GroupLag& g) {
            auto get = [&](void* p, std::size_t n) {
                if (at + n > b.size()) return false;
                std::copy(b.data() + at, b.data() + at + n, static_cast<char*>(p));
                at += n; return true; };
            return get(&g.lag_ms, 4) && get(&g.confidence, 8)
                && get(&g.method, 1) && get(&g.confident, 1)
                && get(&g.n_r_peaks, 4) && get(&g.n_feet, 4)
                && get(&g.analyzed_sec, 8);
        }
    }

    inline bool appendTrailer(const std::string& path, const RecordLag& lag) {
        std::string payload;
        detail::putGroup(payload, lag.ppg);
        detail::putGroup(payload, lag.arterial);
        std::ofstream f(path, std::ios::binary | std::ios::app);
        if (!f) return false;
        const uint32_t n = static_cast<uint32_t>(payload.size());
        f.write(payload.data(), static_cast<std::streamsize>(payload.size()));
        f.write(reinterpret_cast<const char*>(&n), 4);
        f.write(kTrailerMagic, sizeof kTrailerMagic);
        return f.good();
    }

    /// False, and `lag` untouched, when the file has no trailer.
    inline bool readTrailer(const std::string& path, RecordLag& lag) {
        std::ifstream f(path, std::ios::binary | std::ios::ate);
        if (!f) return false;
        const std::streamoff end = f.tellg();
        if (end < 12) return false;
        char magic[8] = {};
        uint32_t n = 0;
        f.seekg(end - 12);
        if (!f.read(reinterpret_cast<char*>(&n), 4) || !f.read(magic, 8)) return false;
        if (!std::equal(magic, magic + 8, kTrailerMagic)) return false;
        if (static_cast<std::streamoff>(n) > end - 12) return false;
        std::string payload(n, '\0');
        f.seekg(end - 12 - static_cast<std::streamoff>(n));
        if (n && !f.read(payload.data(), n)) return false;
        RecordLag out;
        std::size_t at = 0;
        if (!detail::getGroup(payload, at, out.ppg)
            || !detail::getGroup(payload, at, out.arterial)) return false;
        lag = out;
        return true;
    }

    // `group` labels which measurement this row is (e.g. "PPG" or
    // "ARTERIAL"); `append` lets the caller log both groups to the same
    // per-file CSV (first call append=false to truncate+header, second
    // call append=true).
    inline void write_log(const Result& r, const std::string& group, bool append) {
        if (g_dir.empty() || g_stem.empty()) return;
        const std::string path = g_dir + "/" + g_stem + "_channel_offset.csv";
        std::ofstream f(path, append ? std::ios::app : std::ios::trunc);
        if (!f) {
            std::fprintf(stderr, "[channel_offset] cannot open %s\n", path.c_str());
            return;
        }
        // `applied` KEEPS ITS MEANING: were the signals shifted by lag_ms. It
        // is 0 on every row this build writes, because nothing is shifted any
        // more. A log from an earlier build with applied = 1 marks a record
        // whose outputs were measured on shifted signals -- the ones
        // restore_chaos_lag.py adds the lag back to. used_for_window_display
        // is what a confident lag is used for now.
        if (!append)
            f << "stem,group,lag_ms,applied,needs_manual_review,ratio,best_score,"
            "n_r_peaks,n_feet,analyzed_sec,method,used_for_window_display\n";
        f << g_stem << ',' << group << ',' << r.lag_ms << ",0,"
            << (r.ambiguous ? 1 : 0) << ',' << r.ratio << ',' << r.best_score << ','
            << r.n_r_peaks << ',' << r.n_feet << ',' << r.analyzed_sec << ','
            << methodName(static_cast<uint8_t>(LagMethod::RPEAK_FOOT_TRAIN_XCORR)) << ','
            << (r.ambiguous ? 0 : 1) << '\n';
    }

}   // namespace channel_offset#pragma once
