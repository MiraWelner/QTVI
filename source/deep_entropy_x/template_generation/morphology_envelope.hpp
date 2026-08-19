/**
 * @file   morphology_envelope.hpp
 * @brief  Per-sample morphology envelope and band-match scoring.
 *         Spec Section 4.7.1.
 *
 *         Three corrections against the literal transcription:
 *
 *          1. PERCENTILES ARE INTERPOLATED. The truncated form
 *             col[(int)(0.025*n)] floors to index 0 for n < 40 and
 *             col[(int)(0.975*n)] reaches n-1 for n <= 40, so any bin with
 *             20-39 contributing beats produced a min/max corridor rather
 *             than a 95% corridor -- every beat then scored ~100% and the
 *             clean-pool gate passed the whole set, ectopics included.
 *             quantile_sorted() uses the same linear interpolation as
 *             envelopes::quantile_of, so the two envelope families now agree.
 *
 *          2. SMALL-N IS REPORTED, NOT HIDDEN. A 2.5/97.5 corridor estimated
 *             from 20 beats sits on the extremes of a 20-sample distribution
 *             and is unstable however it is indexed. tight() says whether the
 *             corridor had enough contributors to be treated as authoritative
 *             (kTightMinBeats). This is a report, not a policy: callers decide
 *             whether to raise their own minimum or down-weight the bin.
 *
 *          3. ST AND T ARE SEPARATED. The transcription scored
 *             pct_ST = [qrsEnd, tEnd) -- J-point to T-end, i.e. ST *and* T
 *             together -- and pct_T = [tEnd, W), the isoelectric tail after
 *             the T wave. The T-wave column measured flat baseline. scoreBeat
 *             now takes tBegin, so ST = [qrsEnd, tBegin) and T = [tBegin,
 *             tEnd), with the tail kept separately as pct_tail rather than
 *             discarded (it is a useful wander/noise indicator, just not a T
 *             wave).
 *
 *             NOTE: this assumes the spec's fourth scoreBeat parameter, named
 *             tEnd, is the end of the T wave. If it was meant as T-begin the
 *             original transcription was faithful and this change belongs in
 *             the document instead -- worth confirming against the 4.7.1 text.
 *
 *          4. An empty band returns NaN, not 0.0. Zero is indistinguishable
 *             from "every sample out of band", so a bin with an undetected
 *             landmark used to read as a totally failing beat.
 *
 * @date   2026-08-19
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

struct MorphologyEnvelope {
    std::vector<double> lo, hi;       // 2.5th/97.5th percentile at each sample
    std::vector<double> mean, sd;
    std::vector<int>    nContrib;     // non-NaN beats per column
    int nBeats = 0;                   // contributing beat count

    // Below this the percentile estimate rests on too few order statistics to
    // be treated as a tight corridor.
    static constexpr int kTightMinBeats = 40;
    bool tight() const { return nBeats >= kTightMinBeats; }
    bool empty() const { return lo.empty(); }
};

// Linear-interpolated quantile of an ALREADY SORTED, finite-valued vector.
// Same formula as envelopes::quantile_of so the two envelope families report
// comparable numbers.
inline double quantile_sorted(const std::vector<double>& v, double q) {
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    if (v.size() == 1) return v[0];
    const double pos = q * static_cast<double>(v.size() - 1);
    const size_t i = static_cast<size_t>(std::floor(pos));
    const size_t k = std::min(i + 1, v.size() - 1);
    const double f = pos - static_cast<double>(i);
    return v[i] * (1.0 - f) + v[k] * f;
}

// Same value as quantile_sorted, without sorting the whole column. Only two
// order statistics are ever needed per column (the 2.5th and 97.5th
// percentile), so a full O(n log n) sort per column was doing far more work
// than the result required -- and buildEnvelope does this W times per bin,
// twice per bin, for every bin and channel.
//
// nth_element partitions in O(n): after it, v[i] holds the i-th order
// statistic and everything below index i is <= it. The interpolation partner
// v[i+1] is then the SMALLEST element of the upper partition, which is
// min_element over [i+1, end) -- also O(n) and restricted to the part of the
// range that nth_element already isolated.
//
// Destructive: reorders `v`. Callers must take the low percentile before the
// high one, since each call re-partitions. buildEnvelope does.
//
// Same pattern as the column medians in create_arterial_templates.hpp.
inline double quantile_nth(std::vector<double>& v, double q) {
    const size_t n = v.size();
    if (n == 0) return std::numeric_limits<double>::quiet_NaN();
    if (n == 1) return v[0];
    const double pos = q * static_cast<double>(n - 1);
    const size_t i = static_cast<size_t>(std::floor(pos));
    const double f = pos - static_cast<double>(i);

    std::nth_element(v.begin(), v.begin() + i, v.end());
    const double vi = v[i];
    if (f == 0.0 || i + 1 >= n) return vi;
    const double vk = *std::min_element(v.begin() + i + 1, v.end());
    return vi * (1.0 - f) + vk * f;
}

inline MorphologyEnvelope buildEnvelope(const std::vector<std::vector<double>>& cleanBeats,
    int W)
{
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    MorphologyEnvelope env;
    env.nBeats = static_cast<int>(cleanBeats.size());
    if (W <= 0) return env;

    // NaN-filled, not zero-filled: a column with no contributors must be
    // distinguishable from a column whose corridor genuinely straddles zero.
    env.lo.assign(W, NaN);
    env.hi.assign(W, NaN);
    env.mean.assign(W, NaN);
    env.sd.assign(W, NaN);
    env.nContrib.assign(W, 0);

    std::vector<double> col;
    col.reserve(cleanBeats.size());
    for (int c = 0; c < W; ++c) {
        col.clear();
        for (const auto& beat : cleanBeats)
            if (c < static_cast<int>(beat.size()) && !std::isnan(beat[c]))
                col.push_back(beat[c]);
        env.nContrib[c] = static_cast<int>(col.size());
        if (col.empty()) continue;
        // Low percentile first: quantile_nth re-partitions `col` on each call,
        // and taking 0.025 first leaves the upper region intact for 0.975.
        env.lo[c] = quantile_nth(col, 0.025);
        env.hi[c] = quantile_nth(col, 0.975);
        double m = 0.0;
        for (double v : col) m += v;
        m /= static_cast<double>(col.size());
        env.mean[c] = m;
        double s = 0.0;
        for (double v : col) s += (v - m) * (v - m);
        env.sd[c] = std::sqrt(s / std::max(1, static_cast<int>(col.size()) - 1));
    }
    return env;
}

struct BandMatchResult {
    // Percentage of scorable samples inside the corridor. NaN means the band
    // had no scorable samples (landmark missing, or no corridor there) -- a
    // different statement from 0.0, "every sample fell outside".
    double pct_overall = std::numeric_limits<double>::quiet_NaN();
    double pct_P = std::numeric_limits<double>::quiet_NaN();
    double pct_QRS = std::numeric_limits<double>::quiet_NaN();
    double pct_ST = std::numeric_limits<double>::quiet_NaN();
    double pct_T = std::numeric_limits<double>::quiet_NaN();
    double pct_tail = std::numeric_limits<double>::quiet_NaN();  // after T-end
};

// Bands: P [0,pEnd), QRS [qrsStart,qrsEnd), ST [qrsEnd,tBegin),
//        T [tBegin,tEnd), tail [tEnd,W), overall [0,W).
//
// The P band starts at column 0, so it carries whatever pre-P baseline the
// slice holds and is diluted the same way pct_overall is -- see the dilution
// note in premark_beats.hpp. It is left as the spec has it; if P-onset ever
// becomes available, narrow it.
// Two-band variant: overall and QRS only, for callers that read just those.
// The pass-1 clean-pool gate in premark_beats.hpp is one, and it runs over
// every beat in every bin, so the four bands it never reads are worth not
// computing. Fields it does not fill stay NaN.
inline BandMatchResult scoreBeatBands(const std::vector<double>& beat,
    const MorphologyEnvelope& env, int qrsStart, int qrsEnd)
{
    const int W = static_cast<int>(std::min(beat.size(), env.lo.size()));
    auto pctInRange = [&](int lo, int hi) -> double {
        lo = std::max(0, lo);
        hi = std::min(hi, W);
        int in = 0, total = 0;
        for (int i = lo; i < hi; ++i) {
            if (std::isnan(beat[i]) || std::isnan(env.lo[i]) || std::isnan(env.hi[i]))
                continue;
            ++total;
            if (beat[i] >= env.lo[i] && beat[i] <= env.hi[i]) ++in;
        }
        return total > 0 ? 100.0 * in / total
            : std::numeric_limits<double>::quiet_NaN();
        };
    BandMatchResult r;
    r.pct_overall = pctInRange(0, W);
    r.pct_QRS = pctInRange(qrsStart, qrsEnd);
    return r;
}

inline BandMatchResult scoreBeat(const std::vector<double>& beat,
    const MorphologyEnvelope& env,
    int pEnd, int qrsStart, int qrsEnd, int tBegin, int tEnd)
{
    const int W = static_cast<int>(std::min(beat.size(), env.lo.size()));

    auto pctInRange = [&](int lo, int hi) -> double {
        lo = std::max(0, lo);
        hi = std::min(hi, W);
        int in = 0, total = 0;
        for (int i = lo; i < hi; ++i) {
            if (std::isnan(beat[i]) || std::isnan(env.lo[i]) || std::isnan(env.hi[i]))
                continue;                                  // nothing to score here
            ++total;
            if (beat[i] >= env.lo[i] && beat[i] <= env.hi[i]) ++in;
        }
        return total > 0 ? 100.0 * in / total
            : std::numeric_limits<double>::quiet_NaN();
        };

    BandMatchResult r;
    r.pct_overall = pctInRange(0, W);
    r.pct_P = pctInRange(0, pEnd);
    r.pct_QRS = pctInRange(qrsStart, qrsEnd);
    r.pct_ST = pctInRange(qrsEnd, tBegin);
    r.pct_T = pctInRange(tBegin, tEnd);
    r.pct_tail = pctInRange(tEnd, W);
    return r;
}