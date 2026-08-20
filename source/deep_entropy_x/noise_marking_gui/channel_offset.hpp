/**
 * @file   channel_offset.hpp
 * @brief  One-time per-file measurement of the fixed hardware delay between
 *         the ECG and PPG acquisition paths.
 *
 *         Runs once per file during the initial load (prepareViewerJob),
 *         before any template slicing. Builds two event trains -- R-peak
 *         times from the ECG and foot times from the PPG -- pools them
 *         across every bin in the file, and finds the single time shift
 *         that best aligns the second onto the first.
 *
 *         The two trains are in DIFFERENT sample units (rates.ecg vs
 *         rates.ppg), so everything here is in seconds.
 *
 *         Rather than literally sliding the foot train 500 times, this
 *         computes each foot's delay from its preceding R peak and
 *         histograms those delays. The histogram peak is the same
 *         maximum a slide-and-score cross-correlation would find, at one
 *         pass instead of 500.
 */
#pragma once

#include "template_generation/find_foot_pulseox.hpp"
#include "template_generation/pulse_matched_filter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

namespace channel_offset {

    // Search window and resolution. kMaxLagSec must stay well BELOW one
    // RR interval: if the window were wider than a heartbeat, a candidate
    // lag could align each foot with the NEXT beat's R peak and score just
    // as well, with no way to tell the two apart.
    constexpr double kMaxLagSec = 0.500;
    constexpr double kBinSec = 0.005;   // histogram resolution (5 ms)

    // How much of the recording to measure on. The offset is a per-file
    // constant, so the whole file adds cost without adding information --
    // 5 minutes at 60 bpm is ~300 beats, well past kMinFeet. Bins are
    // consumed in order until this budget is spent; the last bin used is
    // truncated rather than skipped.
    constexpr double kMaxAnalysisSec = 300.0;

    // A lag is accepted only if its histogram bin holds at least this
    // fraction of all matched feet, AND beats the best rival bin (one
    // outside the peak's immediate neighbourhood) by this ratio. Together
    // these are the "is there a clear peak" test: a flat plateau or two
    // similar humps fails both.
    constexpr double kMinPeakFraction = 0.30;
    constexpr double kMinPeakRatio = 2.0;
    constexpr int    kMinFeet = 100;

    struct Result {
        double lag_sec = 0.0;    // measured delay, PPG late by this much
        double confidence = 0.0; // peak count / best rival count
        double share = 0.0;      // peak count / total matched feet
        bool   ambiguous = true; // true => do NOT shift, flag for review
        size_t n_r_peaks = 0;
        size_t n_feet = 0;
        size_t n_matched = 0;
        double analyzed_sec = 0.0;   // how much of the file was measured on
    };

    // ---------------------------------------------------------------------
    // Foot census on a CONTINUOUS pulse signal, independent of the ECG.
    //
    // find_foot_pulseox takes a matrix of pre-sliced beats, so the slicing
    // has to come from the PPG's own derivative -- not from R peaks, which
    // would make the two trains dependent and the measured lag meaningless.
    // This is the same sequence create_arterial_templates.hpp uses:
    // derivative-max upstroke census -> apex walk -> [prevPeak, thisPeak]
    // segments -> find_foot_pulseox -> absolute indices.
    // ---------------------------------------------------------------------
    inline std::vector<int> detect_feet(const std::vector<double>& signal,
        double channelRate)
    {
        std::vector<int> footAbs;
        if (signal.empty() || channelRate <= 0.0) return footAbs;
        const int n = static_cast<int>(signal.size());

        const int minSep = std::max(1, static_cast<int>(std::llround(0.25 * channelRate)));
        const std::vector<int> upstrokes =
            pulse_matched_filter::derivativePulseLocations(signal, minSep);
        if (upstrokes.size() < 2) return footAbs;

        // Apex walk: from each upstroke forward to the local maximum,
        // bounded by the next upstroke so we cannot cross into the next beat.
        std::vector<int> peaks;
        peaks.reserve(upstrokes.size());
        for (size_t k = 0; k < upstrokes.size(); ++k) {
            const int start = upstrokes[k];
            const int hardEnd = (k + 1 < upstrokes.size()) ? upstrokes[k + 1] : n;
            const int end = std::min(hardEnd, start + minSep);
            int pk = start;
            double pkVal = -std::numeric_limits<double>::infinity();
            for (int i = start; i < end && i < n; ++i) {
                const double v = signal[i];
                if (std::isnan(v)) continue;
                if (v > pkVal) { pkVal = v; pk = i; }
            }
            if (std::isfinite(pkVal)) peaks.push_back(pk);
        }
        if (peaks.size() < 2) return footAbs;

        const size_t nBeats = peaks.size() - 1;
        std::vector<std::vector<double>> segments(nBeats);
        for (size_t k = 0; k < nBeats; ++k) {
            const int a = peaks[k], b = peaks[k + 1];
            segments[k].assign(signal.begin() + a, signal.begin() + b + 1);
        }
        const FootResult feet = find_foot_pulseox(segments);

        footAbs.resize(nBeats);
        for (size_t k = 0; k < nBeats; ++k)
            footAbs[k] = peaks[k] + static_cast<int>(feet.idx[k]);
        return footAbs;
    }

    // ---------------------------------------------------------------------
    // Pool the delay histogram across every bin in the file and pick the
    // peak. One call per file.
    //
    // Templated on the container rather than taking
    // std::vector<output_binfile_data> directly: that type is defined in a
    // header this one does not include, and a range-for needs the COMPLETE
    // type (a vector can be declared with an incomplete element type, but
    // not iterated -- MSVC reports that as C3536 on the loop). As a template
    // the requirement moves to the instantiation point in post_process.hpp,
    // where every relevant header is already in scope.
    //
    // Bins must expose: bad_segment, ch1.raw, ppgSignal, ecgSignal.
    // ---------------------------------------------------------------------
    template <class Bins>
    inline Result measure(const Bins& bins, double ecgRate, double ppgRate)
    {
        Result out;
        if (ecgRate <= 0.0 || ppgRate <= 0.0) return out;

        const int nBins = static_cast<int>(std::ceil(kMaxLagSec / kBinSec));
        std::vector<uint64_t> hist(static_cast<size_t>(nBins), 0);

        for (const auto& bin : bins) {
            const double budget = kMaxAnalysisSec - out.analyzed_sec;
            if (budget <= 0.0) break;                  // 5 minutes already covered

            if (bin.bad_segment) continue;
            const std::vector<size_t>& rp = bin.ch1.raw;
            if (rp.size() < 2 || bin.ppgSignal.empty()) continue;

            // Bin duration from the ECG, which drives the shared time window.
            const double binSec = static_cast<double>(bin.ecgSignal.size()) / ecgRate;
            const double useSec = std::min(binSec, budget);
            if (useSec <= 0.0) continue;

            // PPG truncated to the budget, then its own foot census.
            const size_t ppgKeep = std::min(bin.ppgSignal.size(),
                static_cast<size_t>(std::llround(useSec * ppgRate)));
            if (ppgKeep < 4) continue;
            const std::vector<double> ppgHead(bin.ppgSignal.begin(),
                bin.ppgSignal.begin() + static_cast<std::ptrdiff_t>(ppgKeep));

            const std::vector<int> feet = detect_feet(ppgHead, ppgRate);
            if (feet.empty()) { out.analyzed_sec += useSec; continue; }

            // R peaks truncated to the same window. ch1.raw is ascending, so
            // this is a prefix.
            const size_t rKeep = static_cast<size_t>(
                std::upper_bound(rp.begin(), rp.end(),
                    static_cast<size_t>(std::llround(useSec * ecgRate))) - rp.begin());
            if (rKeep < 2) { out.analyzed_sec += useSec; continue; }

            out.analyzed_sec += useSec;
            out.n_r_peaks += rKeep;
            out.n_feet += feet.size();

            // R times in seconds, ascending.
            std::vector<double> rSec(rKeep);
            for (size_t i = 0; i < rKeep; ++i)
                rSec[i] = static_cast<double>(rp[i]) / ecgRate;

            // For each foot: delay from the last R peak at or before it.
            for (int f : feet) {
                const double fSec = static_cast<double>(f) / ppgRate;
                const auto it = std::upper_bound(rSec.begin(), rSec.end(), fSec);
                if (it == rSec.begin()) continue;          // no R before this foot
                const double delay = fSec - *(it - 1);
                if (delay < 0.0 || delay >= kMaxLagSec) continue;
                const int b = static_cast<int>(delay / kBinSec);
                if (b >= 0 && b < nBins) { ++hist[static_cast<size_t>(b)]; ++out.n_matched; }
            }
        }

        if (out.n_matched < static_cast<size_t>(kMinFeet)) return out;   // stays ambiguous

        // Peak bin.
        size_t peakIdx = 0;
        uint64_t peakCount = 0;
        for (size_t i = 0; i < hist.size(); ++i)
            if (hist[i] > peakCount) { peakCount = hist[i]; peakIdx = i; }
        if (peakCount == 0) return out;

        // Best rival, excluding the peak and its two immediate neighbours
        // (a real peak spreads a little across adjacent bins; that spread is
        // not a competing answer).
        uint64_t rival = 0;
        for (size_t i = 0; i < hist.size(); ++i) {
            const size_t d = (i > peakIdx) ? (i - peakIdx) : (peakIdx - i);
            if (d <= 2) continue;
            rival = std::max(rival, hist[i]);
        }

        // Bin-centre, then a count-weighted centroid over the peak and its
        // two neighbours for sub-bin resolution.
        double wsum = 0.0, csum = 0.0;
        for (size_t i = (peakIdx >= 2 ? peakIdx - 2 : 0);
            i <= std::min(peakIdx + 2, hist.size() - 1); ++i) {
            const double centre = (static_cast<double>(i) + 0.5) * kBinSec;
            wsum += static_cast<double>(hist[i]);
            csum += static_cast<double>(hist[i]) * centre;
        }
        out.lag_sec = (wsum > 0.0) ? (csum / wsum)
            : ((static_cast<double>(peakIdx) + 0.5) * kBinSec);

        out.share = static_cast<double>(peakCount) / static_cast<double>(out.n_matched);
        out.confidence = (rival > 0)
            ? static_cast<double>(peakCount) / static_cast<double>(rival)
            : static_cast<double>(peakCount);

        out.ambiguous = (out.share < kMinPeakFraction)
            || (out.confidence < kMinPeakRatio);
        return out;
    }

    // ---------------------------------------------------------------------
    // Apply: drop the leading `lag` worth of PPG samples in every bin, so
    // every PPG timestamp effectively moves earlier by lag_sec. No-op when
    // the measurement was ambiguous.
    // ---------------------------------------------------------------------
    template <class Bins>
    inline void apply(Bins& bins, const Result& r, double ppgRate)
    {
        if (r.ambiguous || r.lag_sec <= 0.0 || ppgRate <= 0.0) return;
        const size_t shift = static_cast<size_t>(std::llround(r.lag_sec * ppgRate));
        if (shift == 0) return;
        for (auto& bin : bins) {
            if (bin.ppgSignal.size() <= shift) continue;
            bin.ppgSignal.erase(bin.ppgSignal.begin(),
                bin.ppgSignal.begin() + static_cast<std::ptrdiff_t>(shift));
        }
    }

    // ---------------------------------------------------------------------
    // Log destination. Same dir/stem pattern as ecg_move_log and
    // premark_beats: set once from post_process before the measurement,
    // read-only afterwards, single-threaded writer. Empty dir/stem => no log.
    //
    // A logfile rather than a .bin header field: the header is fixed-size
    // (FILE_HEADER_SIZE static_assert) and version gated, so adding a field
    // there means a BIN_HEADER_VERSION bump and regenerating every bin.
    // ---------------------------------------------------------------------
    inline std::string g_dir, g_stem;
    inline void set(const std::string& dir, const std::string& stem) {
        g_dir = dir; g_stem = stem;
    }

    // One row per file. Carries the stem in column 1 so a batch of these
    // concatenates into a single review table.
    inline void write_log(const Result& r) {
        if (g_dir.empty() || g_stem.empty()) return;
        const std::string path = g_dir + "/" + g_stem + "_channel_offset.csv";
        std::ofstream f(path, std::ios::trunc);
        if (!f) {
            std::fprintf(stderr, "[channel_offset] cannot open %s for writing\n",
                path.c_str());
            return;
        }
        f << "stem,lag_ms,applied,needs_manual_review,peak_share,peak_ratio,"
            "n_r_peaks,n_ppg_feet,n_matched,analyzed_sec\n";
        f << g_stem << ','
            << (r.lag_sec * 1000.0) << ','
            << (r.ambiguous ? 0 : 1) << ','
            << (r.ambiguous ? 1 : 0) << ','
            << r.share << ','
            << r.confidence << ','
            << r.n_r_peaks << ','
            << r.n_feet << ','
            << r.n_matched << ','
            << r.analyzed_sec << '\n';
    }

}   // namespace channel_offset