// ============================================================================
// File: JoinedRR.hpp
// Ensemble R-R detection using multiple weighted algorithms (header-only)
// ============================================================================
#pragma once

#include "rpeakdetect.hpp"
#include "pan_tompkin.hpp"
#include "ecgLms.hpp"
#include "RRSimpleSquared.hpp"
#include "FilterUtils.hpp"

/**
 * @brief  Result from JoinedRR: the accepted R-peak indices plus
 *         the detrended (pre-bandpass) signal that was fed into the
 *         detection algorithms.
 */
struct JoinedRRResult {
    std::vector<std::size_t> peaks;         ///< Accepted R-peak sample indices
};

namespace joinedrr_detail {

    constexpr int diff_range = 8; // was 2 for 256hz, scaled to 1000hz: 2/256*1000 ~ 8

    // Refinement helper: snap detected R-wave locations to the local extremum
    inline vector<size_t> RPeakfromRWave(const vector<double>& ecg,
        const vector<size_t>& rWaveIdx, double fs, bool inverted)
    {
        if (rWaveIdx.empty()) return {};
        if (rWaveIdx.size() == 1) return rWaveIdx;

        vector<double> d;
        for (size_t i = 1; i < rWaveIdx.size(); ++i)
            d.push_back((double)rWaveIdx[i] - rWaveIdx[i - 1]);
        double half_win = median(d) / 6.0;

        vector<size_t> refined = rWaveIdx;
        for (size_t i = 0; i < rWaveIdx.size(); ++i) {
            int start = std::max(0, (int)std::round((double)rWaveIdx[i] - half_win));
            int end = std::min((int)ecg.size() - 1, (int)std::round((double)rWaveIdx[i] + half_win));

            double best_v = inverted ? 1e30 : -1e30;
            size_t best_i = rWaveIdx[i];
            for (int j = start; j <= end; ++j) {
                bool better = inverted ? (ecg[j] < best_v) : (ecg[j] > best_v);
                if (better) { best_v = ecg[j]; best_i = (size_t)j; }
            }
            refined[i] = best_i;
        }
        return refined;
    }

} // namespace joinedrr_detail


/**
 * @brief  Full interface that also returns the pre-bandpass (detrended) signal.
 *
 * @param[in] inverted  True if this channel's lead is inverted (from the
 *                      GUI's "Inverted Lead?" checkbox, persisted through
 *                      the annealed .bin). Determines whether the true
 *                      R-peak is a local max or local min at every
 *                      polarity-sensitive decision point below.
 */
inline JoinedRRResult JoinedRR_full(const vector<double>& ecgSeg, double ecgRate, const std::string fileID, bool inverted) {
    using namespace joinedrr_detail;

    JoinedRRResult result;
    if (ecgSeg.empty() || std_dev(ecgSeg) == 0) return result;

    // 0. Detrend
    vector<double> processedEcg = ecgSeg;
    double mu = 0.0;
    if (!processedEcg.empty()) {
        for (double v : processedEcg) mu += v;
        mu /= processedEcg.size();
        for (auto& val : processedEcg) val -= mu;
    }

    // 1. Run ensemble (6 algorithms)
    vector<vector<size_t>> output(6);
    vector<double> weights = { 0.75, 0.25, 0.25, 1.25, 1.5, 0.75 };

    // Algorithms 0-2 share threshold-independent preprocessing. Original code
    // ran rpeakdetect() three times -- bandpass+filter+medfilt1 once each
    // -- when only the threshold differs. We do the prep once and apply
    // three thresholds. Output is bit-identical.
    auto rpd_prep = rpeakdetect_prep(processedEcg, ecgRate, 0, fileID);
    output[0] = rpeakdetect_apply(rpd_prep, 0.2, inverted).r_peak_index;
    output[1] = rpeakdetect_apply(rpd_prep, 0.1, inverted).r_peak_index;
    output[2] = rpeakdetect_apply(rpd_prep, 0.4, inverted).r_peak_index;

    PanTompkinResult pt_res = pan_tompkin(ecgSeg, ecgRate, 0, fileID);
    output[3].assign(pt_res.qrs_i_raw.begin(), pt_res.qrs_i_raw.end());

    vector<double> b = { 5.0 };
    vector<double> a = { 12.0 };
    output[4] = ecgLms(processedEcg, (int)ecgRate, b, a, 0, fileID);

    // Algorithm 6 depends on median distances from the first 5
    vector<double> dists;
    for (int i = 0; i < 5; ++i) {
        if (output[i].size() > 1) {
            vector<double> d = diff(vector<double>(output[i].begin(), output[i].end()));
            dists.push_back(median(d));
        }
    }
    double m_dist = (!dists.empty()) ? median(dists) : (ecgRate * 0.6);
    output[5] = RRsimpleSquared(ecgSeg, m_dist / 2.0).first;

    // 2. Refine algorithms 4, 5, 6 (indices 3-5, matching MATLAB)
    for (size_t r = 3; r < 6; ++r)
        output[r] = RPeakfromRWave(ecgSeg, output[r], ecgRate, inverted);

    // 3. Build weighted detection list
    struct DetWithWeight { size_t pos; double weight; };
    vector<DetWithWeight> all_weighted;
    for (int i = 0; i < 6; ++i)
        for (size_t p : output[i])
            all_weighted.push_back({ p, weights[i] });

    std::sort(all_weighted.begin(), all_weighted.end(),
        [](const auto& a, const auto& b) { return a.pos < b.pos; });

    if (all_weighted.empty()) return result;

    // 4. Get unique positions
    vector<size_t> uniq;
    for (const auto& dw : all_weighted)
        if (uniq.empty() || uniq.back() != dw.pos)
            uniq.push_back(dw.pos);

    // 5. Merge adjacent positions within diff_range.
    for (size_t i = 0; i + 1 < uniq.size(); ++i) {
        if ((double)(uniq[i + 1] - uniq[i]) <= diff_range) {
            size_t winner, loser;
            bool first_wins = inverted
                ? (ecgSeg[uniq[i]] < ecgSeg[uniq[i + 1]])
                : (ecgSeg[uniq[i]] > ecgSeg[uniq[i + 1]]);
            if (first_wins) {
                winner = uniq[i]; loser = uniq[i + 1];
            }
            else {
                winner = uniq[i + 1]; loser = uniq[i];
            }
            for (auto& dw : all_weighted)
                if (dw.pos == loser) dw.pos = winner;
        }
    }

    // 6. Sum weights per unique position
    std::map<size_t, double> weighted_peaks;
    for (const auto& dw : all_weighted)
        weighted_peaks[dw.pos] += dw.weight;

    // 7. Accept peaks with weight >= 2.4
    vector<size_t> candidates;
    for (const auto& [pos, weight] : weighted_peaks)
        if (weight >= 2.399)
            candidates.push_back(pos);

    std::sort(candidates.begin(), candidates.end());
    result.peaks = candidates;
    return result;
}