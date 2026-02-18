#include "JoinedRR.h"
#include "rpeakdetect.h"
#include "pan_tompkin.h"
#include "ecgLms.h"
#include "RRsimpleSquared.h"
#include "StatsUtils.h"
#include "FilterUtils.h"
#include <algorithm>
#include <vector>
#include <map> // Required for weight merging

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

vector<size_t> JoinedRR(const vector<double>& ecgSeg, double ecgSamplingRate, double diff_range) {
    if (std_dev(ecgSeg) == 0) return {};

    // 1. Run Ensemble
    vector<vector<size_t>> output(6);
    vector<double> weights = { 0.75, 0.25, 0.25, 1.25, 1.5, 0.75 };

    output[0] = rpeakdetect(ecgSeg, ecgSamplingRate).R_index;
    output[1] = rpeakdetect(ecgSeg, ecgSamplingRate, 0.1).R_index;
    output[2] = rpeakdetect(ecgSeg, ecgSamplingRate, 0.4).R_index;
    output[3] = pan_tompkin(ecgSeg, ecgSamplingRate).qrs_i_raw;

    vector<double> b, a;
    butter(3, { 5.0 * 2.0 / ecgSamplingRate, 12.0 * 2.0 / ecgSamplingRate }, b, a);
    output[4] = ecgLms(ecgSeg, (int)ecgSamplingRate, b, a, 0);

    vector<double> dists;
    for (int i = 0; i < 5; ++i) {
        if (output[i].size() > 1) {
            vector<double> d = diff(vector<double>(output[i].begin(), output[i].end()));
            dists.push_back(median(d));
        }
    }
    double m_dist = (!dists.empty()) ? median(dists) : (ecgSamplingRate * 0.6);
    output[5] = RRsimpleSquared(ecgSeg, m_dist / 2.0).first;

    // 2. Refine ALL algorithms to ensure they are looking at the same peak sample
    for (size_t r = 0; r < 6; ++r) {
        output[r] = RPeakfromRWave(ecgSeg, output[r], ecgSamplingRate);
    }

    // 3. Flatten detections with Algorithm ID
    struct PeakDet { size_t pos; double weight; int algo_idx; };
    vector<PeakDet> all_dets;
    for (int i = 0; i < 6; ++i) {
        for (size_t p : output[i]) all_dets.push_back({ p, weights[i], i });
    }
    std::sort(all_dets.begin(), all_dets.end(), [](const PeakDet& a, const PeakDet& b) {
        return a.pos < b.pos;
        });

    // 4. Chained Cluster Merging
    vector<size_t> candidates;
    if (all_dets.empty()) return candidates;

    vector<PeakDet> cluster;
    for (size_t i = 0; i < all_dets.size(); ++i) {
        if (cluster.empty() || (double)all_dets[i].pos - (double)cluster.back().pos <= diff_range) {
            cluster.push_back(all_dets[i]);
        }
        else {
            // Process the finished cluster
            std::map<int, double> algo_weights;
            size_t best_pos = cluster[0].pos;
            double max_v = -1e30;
            for (auto& pd : cluster) {
                if (pd.weight > algo_weights[pd.algo_idx]) algo_weights[pd.algo_idx] = pd.weight;
                if (ecgSeg[pd.pos] > max_v) { max_v = ecgSeg[pd.pos]; best_pos = pd.pos; }
            }
            double total_w = 0;
            for (auto const& [id, w] : algo_weights) total_w += w;

            // USE 2.0: Higher threshold to prevent T-waves/S-waves from being counted
            if (total_w >= 2.0) candidates.push_back(best_pos);

            cluster.clear();
            cluster.push_back(all_dets[i]);
        }
    }
    // Handle last cluster
    if (!cluster.empty()) {
        std::map<int, double> algo_weights;
        size_t best_pos = cluster[0].pos;
        double max_v = -1e30;
        for (auto& pd : cluster) {
            if (pd.weight > algo_weights[pd.algo_idx]) algo_weights[pd.algo_idx] = pd.weight;
            if (ecgSeg[pd.pos] > max_v) { max_v = ecgSeg[pd.pos]; best_pos = pd.pos; }
        }
        double total_w = 0;
        for (auto const& [id, w] : algo_weights) total_w += w;
        if (total_w >= 2.0) candidates.push_back(best_pos);
    }

    // 5. --- REFRACTORY PERIOD FILTER (250ms) ---
    vector<size_t> final_rr;
    double refractory_samples = 0.25 * ecgSamplingRate; // 250ms
    for (size_t p : candidates) {
        if (final_rr.empty() || (p - final_rr.back()) > refractory_samples) {
            final_rr.push_back(p);
        }
        else {
            // If two beats are too close, keep the one with higher ECG voltage
            if (ecgSeg[p] > ecgSeg[final_rr.back()]) {
                final_rr.back() = p;
            }
        }
    }

    return final_rr;
}
