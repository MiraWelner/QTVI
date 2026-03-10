// ============================================================================
// File: FindWaveBounds_EKGandPPG.cpp (MATCHES MATLAB FLOW EXACTLY)
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
        d.ppg_bin_indexs = seg.ppg_bin_indexs;
        d.ecg_bin_indexs = seg.ecg_bin_indexs;

        // 1. Load existing R-peaks if available
        d.ecgRIndex = seg.r_peaks;
        bool rIsNoise = false;

        // 2. Segment PPG (Matches MATLAB try-catch)
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

        // 3. Noise Logic (Line 31-38 in MATLAB)
        if (d.ecgRIndex.empty() && use_R_algorithms) {
            if (std_dev(seg.ecg) == 0) {
                rIsNoise = true;
            }
        }

        // 4. Detection Logic (Line 41-49 in MATLAB)
        // CRITICAL: Only run JoinedRR and Heuristic IF ecgRIndex is still empty
        if (!rIsNoise && d.ecgRIndex.empty() && use_R_algorithms) {
            try {
                d.ecgRIndex = JoinedRR(seg.ecg, seg.ecgSampleRate, 2.0, fileID);

                // Heuristic Check: Only triggers during new detection
                double r_count = static_cast<double>(d.ecgRIndex.size());
                double p_count = static_cast<double>(d.ppgMinAmps.size());
                if (r_count < p_count / 2.0 || p_count * 1.5 < r_count) {
                    rIsNoise = true;
                }
            }
            catch (...) {
                rIsNoise = true;
            }
        }

        // 5. Pairing and Clearing (Line 54-73 in MATLAB)
        if (!rIsNoise && !d.ecgRIndex.empty()) {
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
        else {
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

    return data;
}
