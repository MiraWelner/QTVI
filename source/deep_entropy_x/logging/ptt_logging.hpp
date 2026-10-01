#pragma once
/**
 * @file   ptt_log.hpp
 * @brief  <stem>_ptt.csv: pulse transit time for every ECG/PPG beat pair.
 *
 *         ONE ROW PER PULSE BEAT THE SLICER COULD MEASURE -- every R-pair whose
 *         pulse slice had a systolic upstroke peak, a foot and a 50% crossing
 *         (alignment::extract_ppg_beats_and_align). That is BEFORE the pulse QC
 *         correlation filter, so beats the filter rejected are logged too, with
 *         kept = 0. A pair whose pulse landmarks could not be found has no row:
 *         there is nothing to measure a transit time to.
 *
 *         COLUMNS
 *           stem, bin      the record and bin
 *           slice          R-pair ordinal (rPeaks[slice], rPeaks[slice+1]) the
 *                          beat was sliced from -- the join key shared with the
 *                          ECG beats (ecg_beat_set::slice_index)
 *           kept           1 if the pulse survived the QC filter into the bin's
 *                          template cohort, 0 if it was rejected
 *           r_s            the R peak, seconds from the start of the bin
 *           ptt_peak       ms from that R peak to this pulse's systolic peak
 *           ptt_t50        ms from that R peak to the pulse's 50% point: where
 *                          the upstroke first crosses halfway from its foot to
 *                          its peak (the same up50 the build time-aligns on)
 *
 *         THE TIMES ARE ABSOLUTE, NOT ALIGNED-FRAME COLUMNS. R comes from the
 *         ECG detection at the ECG rate (bin.ch1.raw, the same peaks every
 *         slicer uses); the pulse landmarks are positions in the pulse signal at
 *         the pulse rate. Converting each to seconds on its own clock keeps the
 *         R exact, rather than rounding it onto the pulse grid the way the
 *         slicer's R column is. ptt_t50 is sub-sample; ptt_peak is resolved to
 *         one pulse sample (the peak detector returns a whole column).
 *
 *         Destination set once per record from analysis_job::prepare. Empty
 *         dir/stem => no log. write() is called single-threaded.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace ptt_log {

    // One measured pulse beat, in seconds from the start of its bin.
    struct Beat {
        uint32_t slice = 0;
        bool     kept = false;
        double   r_peak_distance_from_binstart_in_s = std::numeric_limits<double>::quiet_NaN();
        double   peak_s = std::numeric_limits<double>::quiet_NaN();
        double   t50_s = std::numeric_limits<double>::quiet_NaN();
        double   foot_s = std::numeric_limits<double>::quiet_NaN();
        double   maxup_s = std::numeric_limits<double>::quiet_NaN();
        double   notch_s = std::numeric_limits<double>::quiet_NaN();
    };

    inline std::string g_dir;
    inline std::string g_stem;

    // Every pulse channel's beats, stashed by the build (PPG in
    // GenerateTemplatesFast, the arterial channels in
    // buildTemplatesAndBeatsFast) and written together by writeTransit().
    inline std::vector<std::pair<std::string, std::vector<std::vector<Beat>>>> g_transit;
    inline void clearTransit() { g_transit.clear(); }
    inline void stashTransit(const std::string& channel,
        const std::vector<std::vector<Beat>>& perBin) {
        for (auto& e : g_transit) if (e.first == channel) { e.second = perBin; return; }
        g_transit.emplace_back(channel, perBin);
    }
    inline void set(const std::string& dir, const std::string& stem) {
        g_dir = dir; g_stem = stem;
        clearTransit();
    }

    // perBin[bin] = that bin's measured beats, in slice order.
    inline bool write(const std::vector<std::vector<Beat>>& perBin) {
        if (g_dir.empty() || g_stem.empty()) return false;
        std::ofstream f(g_dir + "/" + g_stem + "_ptt.csv", std::ios::trunc);
        if (!f) return false;
        f << "stem,bin,slice,kept,r_peak_distance_from_binstart_in_s,ptt_peak,ptt_t50\n";
        for (std::size_t b = 0; b < perBin.size(); ++b)
            for (const Beat& bt : perBin[b])
                f << g_stem << ',' << b << ',' << bt.slice << ','
                << (bt.kept ? 1 : 0) << ',' << bt.r_peak_distance_from_binstart_in_s << ','
                << (bt.peak_s - bt.r_peak_distance_from_binstart_in_s) * 1000.0 << ','
                << (bt.t50_s - bt.r_peak_distance_from_binstart_in_s) * 1000.0 << '\n';
        return f.good();
    }

    // =======================================================================
    // TRANSIT TIMES, EVERY PULSE CHANNEL: <stem>_transit_beats.csv and
    // <stem>_transit_bins.csv
    // =======================================================================
    //
    // R peak to four landmarks -- foot, maximum upslope, systolic peak,
    // dicrotic notch -- plus the 50% point, on PPG and on each arterial
    // channel, every one defined by alignment::pulseLandmarks so the same
    // column name means the same measurement on every channel. All in ms.
    //
    // RAW TIMES, INCLUDING ANY HARDWARE DELAY. The signals are no longer
    // shifted (channel_offset.hpp), so these are R-to-landmark as recorded. The
    // record's measured lag, method and confidence are in
    // <stem>_channel_offset.csv and on the templates and markings files; to get
    // a physiological transit time, subtract the group's lag where it is
    // confident.
    //
    // _transit_beats.csv: the beat series. One row per measured beat per
    // channel, in slice order; kept = 1 when the beat is in the bin's template
    // cohort (passed the pulse QC).
    //
    // _transit_bins.csv: per bin, per channel, per landmark -- the median and
    // the spread (IQR, and p25 / p75 so it is not only a width) over the KEPT
    // beats, with n. Kept beats only, because those are the beats the bin's
    // template, and every other per-bin pulse number, describes. A landmark a
    // beat lacks (a notch) is simply not in that landmark's n.
    namespace detail {
        inline double quantile(std::vector<double> v, double q) {
            if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
            std::sort(v.begin(), v.end());
            const double pos = q * static_cast<double>(v.size() - 1);
            const std::size_t i = static_cast<std::size_t>(std::floor(pos));
            const double f = pos - static_cast<double>(i);
            return (i + 1 < v.size()) ? v[i] + f * (v[i + 1] - v[i]) : v[i];
        }
        inline double ms(double tS, double rS) {
            return (std::isfinite(tS) && std::isfinite(rS))
                ? (tS - rS) * 1000.0 : std::numeric_limits<double>::quiet_NaN();
        }
    }

    inline bool writeTransit() {
        if (g_dir.empty() || g_stem.empty()) return false;
        std::ofstream fb(g_dir + "/" + g_stem + "_transit_beats.csv", std::ios::trunc);
        std::ofstream fs(g_dir + "/" + g_stem + "_transit_bins.csv", std::ios::trunc);
        if (!fb || !fs) return false;

        fb << "stem,channel,bin,slice,kept,r_s,foot_ms,max_upslope_ms,"
            "peak_ms,notch_ms,t50_ms\n";
        static const char* kNames[5] = { "foot", "max_upslope", "peak", "notch", "t50" };
        fs << "stem,channel,bin,n_beats,n_kept";
        for (const char* nm : kNames)
            fs << ',' << nm << "_median_ms," << nm << "_iqr_ms,"
            << nm << "_p25_ms," << nm << "_p75_ms," << nm << "_n";
        fs << '\n';

        for (const auto& [chan, perBin] : g_transit) {
            for (std::size_t b = 0; b < perBin.size(); ++b) {
                std::vector<double> col[5];
                std::size_t nKept = 0;
                for (const Beat& bt : perBin[b]) {
                    const double v[5] = {
                        detail::ms(bt.foot_s, bt.r_peak_distance_from_binstart_in_s), detail::ms(bt.maxup_s, bt.r_peak_distance_from_binstart_in_s),
                        detail::ms(bt.peak_s, bt.r_peak_distance_from_binstart_in_s), detail::ms(bt.notch_s, bt.r_peak_distance_from_binstart_in_s),
                        detail::ms(bt.t50_s, bt.r_peak_distance_from_binstart_in_s) };
                    fb << g_stem << ',' << chan << ',' << b << ',' << bt.slice << ','
                        << (bt.kept ? 1 : 0) << ',' << bt.r_peak_distance_from_binstart_in_s;
                    for (const double x : v) {
                        fb << ',';
                        if (std::isfinite(x)) fb << x;   // empty, not "nan", when absent
                    }
                    fb << '\n';
                    if (!bt.kept) continue;
                    ++nKept;    
                    for (int k = 0; k < 5; ++k) if (std::isfinite(v[k])) col[k].push_back(v[k]);
                }
                if (perBin[b].empty()) continue;
                fs << g_stem << ',' << chan << ',' << b << ',' << perBin[b].size()
                    << ',' << nKept;
                for (int k = 0; k < 5; ++k) {
                    const double p25 = detail::quantile(col[k], 0.25);
                    const double med = detail::quantile(col[k], 0.50);
                    const double p75 = detail::quantile(col[k], 0.75);
                    auto cell = [&fs](double x) { fs << ','; if (std::isfinite(x)) fs << x; };
                    cell(med); cell(p75 - p25); cell(p25); cell(p75);
                    fs << ',' << col[k].size();
                }
                fs << '\n';
            }
        }
        return fb.good() && fs.good();
    }

}  // namespace ptt_log