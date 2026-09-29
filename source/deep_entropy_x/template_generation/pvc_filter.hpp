#pragma once
/**
 * @file   pvc_filter.hpp
 * @brief  Prematurity filter and 5-of-8 voting, run across ALL
 *         detected beats as a standalone step.
 */

#include <algorithm>
#include <cstdint>
#include <vector>

#include "template_bank.hpp"   // tbank::PvcFilter

namespace pvc_filter {

    // ---------------------------------------------------------------------
    // Section 4.6, VERBATIM from the specification. Do not edit these two to
    // fix a convention mismatch -- adapt the input instead, via
    // toPrecedingIntervals() below.
    // ---------------------------------------------------------------------

    inline bool isPremature(const std::vector<double>& rr, int t) { // rr in ms
        if (t < 10) return false;
        std::vector<double> w(rr.begin() + t - 10, rr.begin() + t);
        std::sort(w.begin(), w.end());
        double med = w[w.size() / 2];
        return rr[t] < 0.80 * med;   // RR(t) < 0.80 * median
    }

    // 5-of-8 voting: flag beat t if >=5 of the surrounding 8 beats are flagged
    inline bool voteFlag(const std::vector<char>& flag, int t) {
        int lo = std::max(0, t - 4), hi = std::min((int)flag.size(), t + 4), c = 0;
        for (int i = lo; i < hi; ++i) c += flag[i];
        return c >= 5;
    }
    inline std::vector<double> toPrecedingIntervals(
        const std::vector<double>& rr_after)
    {
        const size_t n = rr_after.size();
        std::vector<double> rr_pre(n, std::numeric_limits<double>::quiet_NaN());
        for (size_t i = 1; i < n; ++i) rr_pre[i] = rr_after[i - 1];
        return rr_pre;
    }

    // Section 4.6 constants, for reporting only -- the values live inside the
    // verbatim functions above. Named here so counts and comments can refer to
    // them without a second definition drifting from the first.
    inline constexpr double kPrematureRatio = 0.80;
    inline constexpr int    kTrailingWindow = 10;
    inline constexpr int    kVoteWindow = 8;
    inline constexpr int    kVoteRequired = 5;
    inline constexpr size_t kMinBeatsForFilter = 12;

    struct FilterResult {
        // Parallel to the detected beat list, never compacted.
        std::vector<tbank::PvcFilter> verdict;

        // The ratio the test compared, per beat: the continuous form of a binary
        // decision. A run of beats sitting at 0.82 is a different situation from
        // a run at 0.40, and only the ratio shows it. NaN where no trailing
        // median was available.
        std::vector<double> ratio;

        uint32_t n_premature = 0;
        uint32_t n_vote_only = 0;   // voted but not itself premature
        uint32_t n_no_median = 0;   // t < 10, no trailing window
        uint32_t n_vote_blind = 0;   // vote window clamped below 5 beats
    };

    // rr_after: alignment.hpp's rr_lens convention, in milliseconds.
    inline FilterResult runFilter(const std::vector<double>& rr_after) {
        const size_t n = rr_after.size();
        FilterResult out;
        out.verdict.assign(n, tbank::PvcFilter::NONE);
        out.ratio.assign(n, std::numeric_limits<double>::quiet_NaN());
        if (n == 0) return out;
        if (n < kMinBeatsForFilter) {
            // alignment.hpp requires nb >= 12 before running either test.
            // Matched rather than relaxed: a filter that fires on short records
            // here and not there would make the two disagree exactly where a bin
            // is too short to trust anyway.
            out.n_no_median = static_cast<uint32_t>(n);
            return out;
        }

        const std::vector<double> rr = toPrecedingIntervals(rr_after);

        // ---- pass A: the spec's isPremature(), per beat -------------------
        std::vector<char> flag(n, 0);
        for (size_t t = 0; t < n; ++t) {
            if (t < static_cast<size_t>(kTrailingWindow)) { ++out.n_no_median; continue; }
            if (std::isnan(rr[t])) continue;

            // Ratio recorded over the same window the spec's function uses, so
            // the two cannot disagree about which beats were in the median.
            std::vector<double> w(rr.begin() + t - kTrailingWindow, rr.begin() + t);
            w.erase(std::remove_if(w.begin(), w.end(),
                [](double v) { return std::isnan(v); }), w.end());
            if (!w.empty()) {
                std::sort(w.begin(), w.end());
                const double med = w[w.size() / 2];
                if (med > 0.0) out.ratio[t] = rr[t] / med;
            }

            if (isPremature(rr, static_cast<int>(t))) flag[t] = 1;
        }

        // ---- pass B: the spec's voteFlag(), over pass A's flags -----------
        // A SECOND PASS over the completed array, never interleaved: a vote that
        // could see votes would propagate a flag along the whole record from a
        // single dense neighbourhood.
        for (size_t t = 0; t < n; ++t) {
            const long lo = std::max<long>(0, static_cast<long>(t) - 4);
            const long hi = std::min<long>(static_cast<long>(n),
                static_cast<long>(t) + 4);
            // Near the record edges the spec's own clamping leaves fewer than
            // kVoteRequired beats in the window, so the vote cannot fire there
            // however ectopic the neighbourhood is. Preserved rather than
            // "fixed" -- alignment.hpp clamps identically, and a silently
            // different edge rule between the two would be worse than a known
            // blind spot -- but counted, so the blind region shows up in the
            // archive instead of being discovered later.
            if (hi - lo < kVoteRequired) ++out.n_vote_blind;

            if (flag[t]) {
                out.verdict[t] = tbank::PvcFilter::PREMATURE;
                ++out.n_premature;
            }
            else if (voteFlag(flag, static_cast<int>(t))) {
                // Premature wins over voted: direct evidence over inferred.
                // Same precedence as create_ecg_templates.hpp lines 170-172.
                out.verdict[t] = tbank::PvcFilter::VOTE;
                ++out.n_vote_only;
            }
        }

        return out;
    }


}  // namespace pvc_filter