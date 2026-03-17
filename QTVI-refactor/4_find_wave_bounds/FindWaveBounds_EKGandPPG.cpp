// ============================================================================
// File: FindWaveBounds_EKGandPPG.cpp
// ============================================================================
#include "FindWaveBounds_EKGandPPG.h"
#include "SegmentPPG.h"
#include "JoinedRR.h"
#include "pairRtoPPGBeat.h"
#include "StatsUtils.h"

std::vector<WaveData> FindWaveBounds_EKGandPPG(const std::vector<AnnealedSegment>& annealedSegments,
    int dbg_plot, bool use_R_algorithms, std::string fileID) {
    std::vector<WaveData> data(annealedSegments.size());
    for (std::size_t i = 0; i < annealedSegments.size(); ++i) {
        const auto& seg = annealedSegments[i];
        auto& d = data[i];
        d.index = static_cast<int>(i);
        d.ecgSamplingRate = seg.ecgSampleRate;
        d.ppgSamplingRate = seg.ppgSampleRate;
        d.ppgSignal = seg.po;
        d.ecgSignal = seg.ecg;
        d.ecgSignal2 = seg.ecg2;
        d.ecgSignal3 = seg.ecg3;
        d.ppg_bin_indexs = seg.ppg_bin_indexs;
        d.ecg_bin_indexs = seg.ecg_bin_indexs;

        bool hasPPG = !seg.po.empty();
        bool hasEcg2 = !seg.ecg2.empty();
        bool hasEcg3 = !seg.ecg3.empty();
        bool rIsNoise = false;

        // 1. Segment PPG only if we have PPG data
        if (hasPPG) {
            try {
                SegmentPPGResult ppgResult = SegmentPPG(seg.po, seg.ppgSampleRate);
                d.ppgMinAmps = ppgResult.ppgMinAmps;
                d.ppgMaxAmps = ppgResult.maxAmps;
                d.bad_segment = false;
            }
            catch (...) {
                d.ppgMinAmps.clear();
                d.ppgMaxAmps.clear();
                d.bad_segment = true;
            }
        }
        else {
            d.ppgMinAmps.clear();
            d.ppgMaxAmps.clear();
            d.bad_segment = false;
        }

        // 2. Noise Logic (channel 1)
        if (d.ecgRIndex.empty() && use_R_algorithms) {
            if (std_dev(seg.ecg) == 0) {
                rIsNoise = true;
            }
        }

        // ================================================================
        // 3. Detection Logic — Channel 1 (original behavior)
        // ================================================================
        if (!rIsNoise && d.ecgRIndex.empty() && use_R_algorithms) {
            try {
                d.ecgRIndex = JoinedRR(seg.ecg, seg.ecgSampleRate, 2.0, fileID);
                if (hasPPG) {
                    double r_count = static_cast<double>(d.ecgRIndex.size());
                    double p_count = static_cast<double>(d.ppgMinAmps.size());
                    if (r_count < p_count / 2.0 || p_count * 1.5 < r_count) {
                        rIsNoise = true;
                    }
                }
            }
            catch (...) {
                rIsNoise = true;
            }
        }

        // ================================================================
        // 3b. Detection Logic — Channel 2 (independent, no PPG validation)
        // ================================================================
        if (hasEcg2 && use_R_algorithms) {
            if (std_dev(seg.ecg2) != 0) {
                try {
                    d.ecgRIndex2 = JoinedRR(seg.ecg2, seg.ecgSampleRate, 2.0, fileID);
                }
                catch (...) {
                    d.ecgRIndex2.clear();
                }
            }
        }

        // ================================================================
        // 3c. Detection Logic — Channel 3 (independent, no PPG validation)
        // ================================================================
        if (hasEcg3 && use_R_algorithms) {
            if (std_dev(seg.ecg3) != 0) {
                try {
                    d.ecgRIndex3 = JoinedRR(seg.ecg3, seg.ecgSampleRate, 2.0, fileID);
                }
                catch (...) {
                    d.ecgRIndex3.clear();
                }
            }
        }

        // ================================================================
        // 4. Pairing and Clearing (channel 1 only, same as before)
        // ================================================================
        if (!rIsNoise && !d.ecgRIndex.empty() && hasPPG) {
            try {
                d.pairs = pairRtoPPGBeat(seg.ecg, seg.po, seg.ecgSampleRate,
                    seg.ppgSampleRate, d.ecgRIndex, d.ppgMinAmps);
            }
            catch (...) {
                if (!d.bad_segment) {
                    d.ecgRIndex.clear();
                    d.pairs.clear();
                    for (size_t k = 0; k < d.ppgMinAmps.size(); ++k) {
                        d.pairs.push_back({ (double)d.ppgMinAmps[k], -1.0 });
                    }
                }
                else {
                    d.pairs.clear();
                }
            }
        }
        else if (!rIsNoise && !d.ecgRIndex.empty() && !hasPPG) {
            // ECG-only: keep R-peaks, no pairing needed
            d.pairs.clear();
        }
        else {
            if (!d.bad_segment && hasPPG) {
                d.ecgRIndex.clear();
                d.pairs.clear();
                for (size_t k = 0; k < d.ppgMinAmps.size(); ++k) {
                    d.pairs.push_back({ (double)d.ppgMinAmps[k], -1.0 });
                }
            }
            else {
                d.ecgRIndex.clear();
                d.pairs.clear();
            }
        }
    }
    return data;
}