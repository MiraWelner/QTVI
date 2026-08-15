/**
 * @file   morphology_envelope.hpp
 * @brief  Per-sample morphology envelope and band-match scoring.
 *         Spec Section 4.7.1, transcribed as written.
 *
 * @date   2026-08-14
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

struct MorphologyEnvelope {
    std::vector<double> lo, hi;       // 2.5th/97.5th percentile at each sample
    std::vector<double> mean, sd;
    int nBeats = 0;                   // contributing beat count
};

inline MorphologyEnvelope buildEnvelope(const std::vector<std::vector<double>>& cleanBeats,
    int W)
{
    MorphologyEnvelope env; env.nBeats = cleanBeats.size();
    env.lo.resize(W); env.hi.resize(W); env.mean.resize(W); env.sd.resize(W);
    for (int c = 0; c < W; ++c) {
        std::vector<double> col;
        for (auto& beat : cleanBeats)
            if (c < (int)beat.size() && !std::isnan(beat[c])) col.push_back(beat[c]);
        if (col.empty()) continue;
        std::sort(col.begin(), col.end());
        env.lo[c] = col[(int)(0.025 * col.size())];
        env.hi[c] = col[(int)(0.975 * col.size())];
        double m = 0; for (double v : col) m += v; m /= col.size();
        env.mean[c] = m;
        double s = 0; for (double v : col) s += (v - m) * (v - m);
        env.sd[c] = std::sqrt(s / std::max(1, (int)col.size() - 1));
    }
    return env;
}

struct BandMatchResult {
    double pct_overall;
    double pct_P, pct_QRS, pct_ST, pct_T;
};

inline BandMatchResult scoreBeat(const std::vector<double>& beat,
    const MorphologyEnvelope& env,
    int pEnd, int qrsStart, int qrsEnd, int tEnd)
{
    auto pctInRange = [&](int lo, int hi) {
        int in = 0, total = 0;
        for (int i = lo; i < hi && i < (int)beat.size(); ++i) {
            if (beat[i] >= env.lo[i] && beat[i] <= env.hi[i]) ++in;
            ++total;
        }
        return total > 0 ? 100.0 * in / total : 0.0;
        };
    int W = beat.size();
    return { pctInRange(0, W), pctInRange(0, pEnd), pctInRange(qrsStart, qrsEnd),
             pctInRange(qrsEnd, tEnd), pctInRange(tEnd, W) };
}