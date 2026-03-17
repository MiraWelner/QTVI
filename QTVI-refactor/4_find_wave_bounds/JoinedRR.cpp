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
#include <iostream>

// ============================================================================
// DEBUG CONFIG — set to the bin index you want to inspect, -1 to disable
// ============================================================================
static const int DEBUG_BIN = 1;      // which bin to dump (0-based)
static int g_bin_counter = -1;       // incremented each call

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

// Helper: dump peaks near a time range for debugging
static void dump_algo_peaks(const std::string& name, const vector<size_t>& peaks,
    double fs, double t_start, double t_end) {
    std::cout << "  " << name << " (" << peaks.size() << " total): peaks in ["
        << t_start << "s, " << t_end << "s]: ";
    int count = 0;
    for (size_t p : peaks) {
        double t = (double)p / fs;
        if (t >= t_start && t <= t_end) {
            std::cout << p << "(" << std::fixed << std::setprecision(2) << t << "s) ";
            count++;
        }
    }
    if (count == 0) std::cout << "NONE";
    std::cout << std::endl;
}

vector<size_t> JoinedRR(const vector<double>& ecgSeg, double ecgSamplingRate, double diff_range, const std::string fileID) {
    if (ecgSeg.empty() || std_dev(ecgSeg) == 0) return {};

    g_bin_counter++;

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

    output[0] = rpeakdetect(processedEcg, ecgSamplingRate, 0.2, 0, fileID).R_index;
    output[1] = rpeakdetect(processedEcg, ecgSamplingRate, 0.1, 0, fileID).R_index;
    output[2] = rpeakdetect(processedEcg, ecgSamplingRate, 0.4, 0, fileID).R_index;

    PanTompkinResult pt_res = pan_tompkin(ecgSeg, ecgSamplingRate, 0, fileID);
    output[3].assign(pt_res.qrs_i_raw.begin(), pt_res.qrs_i_raw.end());

    vector<double> b = { 5.0 };
    vector<double> a = { 12.0 };
    output[4] = ecgLms(processedEcg, (int)ecgSamplingRate, b, a, 0, fileID);

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

    // 3. Build weighted detection list
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

    // 5. Merge adjacent unique positions
    {
        std::ofstream dbg(fileID + "_boundary_debug.csv", std::ios::app);
        if (dbg.is_open() && uniq.size() >= 2) {
            dbg << "  PRE_MERGE first gaps: ";
            for (size_t j = 0; j + 1 < std::min((size_t)8, uniq.size()); ++j) {
                dbg << uniq[j] << "->" << uniq[j + 1] << "(d=" << (uniq[j + 1] - uniq[j]) << ") ";
            }
            dbg << "\n  PRE_MERGE last gaps: ";
            size_t start = (uniq.size() > 8) ? uniq.size() - 8 : 0;
            for (size_t j = start; j + 1 < uniq.size(); ++j) {
                dbg << uniq[j] << "->" << uniq[j + 1] << "(d=" << (uniq[j + 1] - uniq[j]) << ") ";
            }
            dbg << "\n";
        }
    }
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
            for (auto& dw : all_weighted) {
                if (dw.pos == loser) dw.pos = winner;
            }
        }
    }

    // 6. Recalculate unique and sum weights
    std::map<size_t, double> weighted_peaks;
    for (const auto& dw : all_weighted) {
        weighted_peaks[dw.pos] += dw.weight;
    }

    // 7. Filter by weight threshold
    vector<size_t> candidates;
    for (const auto& [pos, weight] : weighted_peaks) {
        if (weight >= 2.399) {
            candidates.push_back(pos);
        }
    }
    std::sort(candidates.begin(), candidates.end());

    // ========================================================================
        // DEBUG: Boundary peak analysis
        // ========================================================================
    {
        std::string dbg_fname = fileID + "_boundary_debug.csv";
        std::ofstream dbg(dbg_fname, std::ios::app);
        if (dbg.is_open()) {
            dbg << "--- BIN " << g_bin_counter
                << " | N_samples=" << ecgSeg.size()
                << " | fs=" << ecgSamplingRate
                << " ---\n";

            std::string alg_names[] = {
                "rpeakdetect_0.2", "rpeakdetect_0.1", "rpeakdetect_0.4",
                "pan_tompkin", "ecgLms", "RRsimpleSquared"
            };

            for (int i = 0; i < 6; ++i) {
                dbg << "  " << alg_names[i]
                    << " (w=" << std::fixed << std::setprecision(2) << weights[i]
                    << ", count=" << output[i].size() << "): ";

                if (output[i].empty()) {
                    dbg << "NONE";
                }
                else {
                    size_t print_head = std::min((size_t)3, output[i].size());
                    for (size_t j = 0; j < print_head; ++j) {
                        dbg << output[i][j] << " ";
                    }
                    if (output[i].size() > 6) dbg << "... ";
                    if (output[i].size() > 3) {
                        size_t start_tail = output[i].size() - 3;
                        if (start_tail < print_head) start_tail = print_head;
                        for (size_t j = start_tail; j < output[i].size(); ++j) {
                            dbg << output[i][j] << " ";
                        }
                    }
                }
                dbg << "\n";
            }

            dbg << "  ACCEPTED (" << candidates.size() << "): ";
            if (!candidates.empty()) {
                size_t print_head = std::min((size_t)3, candidates.size());
                for (size_t j = 0; j < print_head; ++j) {
                    dbg << candidates[j] << " ";
                }
                if (candidates.size() > 6) dbg << "... ";
                if (candidates.size() > 3) {
                    size_t start_tail = candidates.size() - 3;
                    if (start_tail < print_head) start_tail = print_head;
                    for (size_t j = start_tail; j < candidates.size(); ++j) {
                        dbg << candidates[j] << " ";
                    }
                }
            }
            dbg << "\n";

            dbg << "  REJECTED_NEAR_BOUNDARIES: ";
            bool any_rejected = false;
            for (const auto& [pos, weight] : weighted_peaks) {
                if (weight < 2.399 && (pos < 500 || pos + 500 > ecgSeg.size())) {
                    dbg << pos << "(w=" << std::fixed << std::setprecision(2)
                        << weight << ") ";
                    any_rejected = true;
                }
            }
            if (!any_rejected) dbg << "NONE";
            dbg << "\n\n";
            dbg.flush();
        }
    }

    return candidates;
}