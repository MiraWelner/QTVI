#include "JoinedRR.h"
#include "rpeakdetect.h"
#include "pan_tompkin.h"
#include "ecgLms.h"
#include "RRsimpleSquared.h"
#include "StatsUtils.h"
#include "FilterUtils.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>

// Refinement helper (ensures algorithms snap to the local maximum)
vector<size_t> RPeakfromRWave(const vector<double>& ecg, const vector<size_t>& rWaveIdx, double fs) {
    if (rWaveIdx.empty()) return {};
    if (rWaveIdx.size() == 1) return rWaveIdx;

    vector<double> d;
    for (size_t i = 1; i < rWaveIdx.size(); ++i) d.push_back((double)rWaveIdx[i] - rWaveIdx[i - 1]);
    double half_win = median(d) / 6.0;

    vector<size_t> refined = rWaveIdx;
    for (size_t i = 0; i < rWaveIdx.size(); ++i) {
        int start = std::max(0, (int)std::round((double)rWaveIdx[i] - half_win));
        int end = std::min((int)ecg.size() - 1, (int)std::round((double)rWaveIdx[i] + half_win));

        double max_v = -1e30;
        size_t max_i = rWaveIdx[i];
        for (int j = start; j <= end; ++j) {
            if (ecg[j] > max_v) { max_v = ecg[j]; max_i = (size_t)j; }
        }
        refined[i] = max_i;
    }
    return refined;
}

vector<size_t> JoinedRR(const vector<double>& ecgSeg, double ecgSamplingRate, double diff_range, const std::string fileID) {
    if (ecgSeg.empty() || std_dev(ecgSeg) == 0) return {};

    // 0. Proper Detrending
    vector<double> processedEcg = ecgSeg;
    double mu = 0.0;
    if (!processedEcg.empty()) {
        for (double v : processedEcg) mu += v;
        mu /= processedEcg.size();
        for (auto& val : processedEcg) val -= mu;
    }


    // 1. Run Ensemble
    vector<vector<size_t>> output(6);
    vector<double> weights = { 0.75, 0.25, 0.25, 1.25, 1.5, 0.75 };

    // Pass the detrended signal to all algorithms
    output[0] = rpeakdetect(processedEcg, ecgSamplingRate, 0.2, 0, fileID).R_index;
    output[1] = rpeakdetect(processedEcg, ecgSamplingRate, 0.1, 0, fileID).R_index;
    output[2] = rpeakdetect(processedEcg, ecgSamplingRate, 0.4, 0, fileID).R_index;

    PanTompkinResult pt_res = pan_tompkin(processedEcg, ecgSamplingRate, 0);
    output[3].assign(pt_res.qrs_i_raw.begin(), pt_res.qrs_i_raw.end());

    vector<double> b = { 5.0 };
    vector<double> a = { 12.0 };
    output[4] = ecgLms(processedEcg, (int)ecgSamplingRate, b, a, 0);

    vector<double> dists;
    for (int i = 0; i < 5; ++i) {
        if (output[i].size() > 1) {
            vector<double> d = diff(vector<double>(output[i].begin(), output[i].end()));
            dists.push_back(median(d));
        }
    }
    double m_dist = (!dists.empty()) ? median(dists) : (ecgSamplingRate * 0.6);
    output[5] = RRsimpleSquared(ecgSeg, m_dist / 2.0).first;

    // 2. Refine ONLY algorithms 4, 5, 6 (matching MATLAB)
    for (size_t r = 3; r < 6; ++r) {
        output[r] = RPeakfromRWave(ecgSeg, output[r], ecgSamplingRate);
    }

    // 3. Build weighted detection list (MATLAB: sortedList)
    struct DetWithWeight { size_t pos; double weight; };
    vector<DetWithWeight> all_weighted;
    for (int i = 0; i < 6; ++i) {
        for (size_t p : output[i]) {
            all_weighted.push_back({ p, weights[i] });
        }
    }
    std::sort(all_weighted.begin(), all_weighted.end(), [](const auto& a, const auto& b) {
        return a.pos < b.pos;
    });

    if (all_weighted.empty()) return {};

    // 4. Get unique positions
    vector<size_t> uniq;
    for (const auto& dw : all_weighted) {
        if (uniq.empty() || uniq.back() != dw.pos) {
            uniq.push_back(dw.pos);
        }
    }


    // Inside JoinedRR.cpp, right before the merger (Step 5)
    std::ofstream alg_dbg(fileID + "_alg_indices.csv", std::ios::app);
    alg_dbg << "--- BIN_START ---\n";
    for (int i = 0; i < 6; ++i) {
        alg_dbg << "Algorithm_" << i << "_Weight_" << weights[i] << ":";
        for (size_t p : output[i]) alg_dbg << p << ",";
        alg_dbg << "\n";
    }
    alg_dbg.close();


    // 5. Merge adjacent unique positions (FIXED: Replicates MATLAB's in-place logic)
 // This ensures that clusters (e.g. indices 10, 11, 12) all merge into a single peak.
    for (size_t i = 0; i + 1 < uniq.size(); ++i) {
        if ((double)(uniq[i + 1] - uniq[i]) <= diff_range) {
            size_t winner, loser;
            if (ecgSeg[uniq[i]] > ecgSeg[uniq[i + 1]]) {
                winner = uniq[i];
                loser = uniq[i + 1];
            }
            else {
                winner = uniq[i + 1];
                loser = uniq[i];
            }

            // MATLAB: potentialPeaks(potentialPeaks(:,1)==loser, 1) = winner;
            // We must update EVERY instance in our weighted list immediately.
            for (auto& dw : all_weighted) {
                if (dw.pos == loser) dw.pos = winner;
            }

            // We must also update the 'uniq' list so the NEXT loop iteration 
            // compares the correct winner against the next peak.
            for (size_t k = 0; k < uniq.size(); ++k) {
                if (uniq[k] == loser) uniq[k] = winner;
            }
        }
    }

    // 6. Recalculate unique and sum weights
    std::map<size_t, double> weighted_peaks;
    for (const auto& dw : all_weighted) {
        weighted_peaks[dw.pos] += dw.weight;
    }

    // 7. Filter by weight threshold (Floating point safe: 2.39)
    vector<size_t> candidates;
    for (const auto& [pos, weight] : weighted_peaks) {
        if (weight >= 2.399) {
            candidates.push_back(pos);
        }
    }
    std::sort(candidates.begin(), candidates.end());

    // ============================================================================
    // DEBUG BLOCK: Final verification right before output
    // ============================================================================
    std::string dbg_filename = fileID + "_joined_debug.csv";
    std::ofstream dbg(dbg_filename, std::ios::app);
    if (dbg.is_open()) {
        // Each call represents one bin. This logs EVERY potential peak before thresholding.
        dbg << "--- NEW_BIN_START ---\n";
        dbg << "Position,SumWeight,Result\n";
        for (auto const& [pos, weight] : weighted_peaks) {
            dbg << pos << "," << std::fixed << std::setprecision(4) << weight << ","
                << (weight >= 2.399 ? "ACCEPTED" : "REJECTED") << "\n";
        }
        dbg << "TOTAL_CANDIDATES_IN_BIN," << candidates.size() << "\n";
        dbg.close();
    }
    // ============================================================================

    return candidates;
}
