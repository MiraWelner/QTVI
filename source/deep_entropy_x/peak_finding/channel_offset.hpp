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
        bool   ambiguous = true;  // true => do NOT shift, flag for review
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

    // Shift PPG earlier by its OWN measured lag (ppgResult), and shift
    // ABP/ART/ART_PULM (slots 33/34/35 in the generic pass-through,
    // all_upsampled -- none of the three has a dedicated field the way PPG
    // does) earlier by a SEPARATE, independently-measured arterial-group
    // lag (arterialResult). PPG and the arterial group are NOT assumed to
    // share a lag -- they have different acquisition paths and different
    // hardware delays -- so each is measured on its own train and applied
    // with its own Result, not reused from the other. Within the arterial
    // group, all three channels share ONE lag (arterialResult, typically
    // measured on ABP -- see the caller) since they share one acquisition
    // path. Each channel's shift amount is still recomputed at THAT
    // channel's own sample rate, not reused as a sample count, since the
    // channels are not necessarily sampled at the same rate. A rate <= 0.0
    // for a given channel means it's absent for this dataset; that channel
    // alone is skipped, the others still shift. No-op per group when that
    // group's Result is ambiguous.
    //
    // This runs on the ANNEALED SEGMENTS, before create_ecg_ppg_pairs_raw --
    // i.e. before SegmentPPG, before R-pairing, before templating. Everything
    // derived from any of these channels is then computed from the shifted
    // signal, so no index arrays need renumbering and nothing can fall out
    // of sync.
    //
    // Segments must expose: ppg_signal, all_upsampled (for the ABP/ART/
    // ART_PULM slots).
    template <class Segments>
    inline void apply(Segments& segs,
        const Result& ppgResult, double ppgRate,
        const Result& arterialResult,
        double abpRate = 0.0, double artRate = 0.0, double artPulmRate = 0.0)
    {
        if (!ppgResult.ambiguous && ppgResult.lag_ms > 0 && ppgRate > 0.0) {
            const size_t shift = static_cast<size_t>(
                std::llround(ppgResult.lag_ms * 0.001 * ppgRate));
            if (shift != 0) {
                for (auto& seg : segs) {
                    if (seg.ppg_signal.size() <= shift) continue;
                    seg.ppg_signal.erase(seg.ppg_signal.begin(),
                        seg.ppg_signal.begin() + static_cast<std::ptrdiff_t>(shift));
                }
            }
        }

        if (!arterialResult.ambiguous && arterialResult.lag_ms > 0) {
            // Slots match file_to_bin / gui_handler: CH_ABP=33, CH_ART=34,
            // CH_ART_PULM=35.
            auto shiftSlot = [&](size_t slot, double rate) {
                if (rate <= 0.0) return;
                const size_t shift = static_cast<size_t>(
                    std::llround(arterialResult.lag_ms * 0.001 * rate));
                if (shift == 0) return;
                for (auto& seg : segs) {
                    if (seg.all_upsampled.size() <= slot) continue;
                    auto& sig = seg.all_upsampled[slot];
                    if (sig.size() <= shift) continue;
                    sig.erase(sig.begin(), sig.begin() + static_cast<std::ptrdiff_t>(shift));
                }
                };
            shiftSlot(33, abpRate);
            shiftSlot(34, artRate);
            shiftSlot(35, artPulmRate);
        }
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
        if (!append)
            f << "stem,group,lag_ms,applied,needs_manual_review,ratio,best_score,"
            "n_r_peaks,n_feet,analyzed_sec\n";
        f << g_stem << ',' << group << ',' << r.lag_ms << ',' << (r.ambiguous ? 0 : 1) << ','
            << (r.ambiguous ? 1 : 0) << ',' << r.ratio << ',' << r.best_score << ','
            << r.n_r_peaks << ',' << r.n_feet << ',' << r.analyzed_sec << '\n';
    }

}   // namespace channel_offset