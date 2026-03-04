#include "FindWaveBounds_EKGandPPG.h"
#include "SegmentPPG.h"
#include "JoinedRR.h"
#include "pairRtoPPGBeat.h"
#include "StatsUtils.h"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

std::vector<WaveData> FindWaveBounds_EKGandPPG(const std::vector<AnnealedSegment>& annealedSegments,
    int dbg_plot,
    bool use_R_algorithms,
    std::string fileID) {

    if (annealedSegments.empty()) return {};

    std::vector<WaveData> data(annealedSegments.size());

    // Loop through each 30-second bin individually (Matlab-style)
    for (std::size_t i = 0; i < annealedSegments.size(); ++i) {
        const auto& seg = annealedSegments[i];

        // --- 1. Basic Metadata & Waveform Storage ---
        data[i].index = static_cast<int>(i);
        data[i].ecgSamplingRate = seg.ecgSampleRate;
        data[i].ppgSamplingRate = seg.ppgSampleRate;

        // These assignments ensure saveWaveData finds the actual signal arrays
        data[i].ppgSignal = seg.po;
        data[i].ecgSignal = seg.ecg;

        data[i].ppg_bin_indexs = seg.ppg_bin_indexs;
        data[i].ecg_bin_indexs = seg.ecg_bin_indexs;

        bool rIsNoise = false;

        // --- 2. Segment PPG (Dips and Peaks) ---
        try {
            // SegmentPPG is called once per 30s bin
            SegmentPPGResult ppgResult = SegmentPPG(seg.po, seg.ppgSampleRate);
            data[i].ppgMinAmps = ppgResult.ppgMinAmps;
            data[i].ppgMaxAmps = ppgResult.maxAmps;
            data[i].bad_segment = false;
        }
        catch (...) {
            data[i].bad_segment = true;
            data[i].ppgMinAmps.clear();
            data[i].ppgMaxAmps.clear();
        }

        // --- 3. Noise Check (Matches Matlab logic) ---
        if (use_R_algorithms && std_dev(seg.ecg) == 0) {
            rIsNoise = true;
        }

        // --- 4. R-Peak Detection (JoinedRR) ---
        if (!rIsNoise && use_R_algorithms) {
            try {
                // JoinedRR runs on a single 30s bin for maximum speed
                data[i].ecgRIndex = JoinedRR(seg.ecg, seg.ecgSampleRate, 20.0, fileID);

                // Noise heuristic: R-peaks should roughly correlate with PPG heart rate
                double r_count = static_cast<double>(data[i].ecgRIndex.size());
                double p_count = static_cast<double>(data[i].ppgMinAmps.size());

                if (r_count < p_count / 2.0 || p_count * 1.5 < r_count) {
                    rIsNoise = true;
                }
            }
            catch (...) {
                rIsNoise = true;
            }
        }

        // --- 5. Pairing R-peaks to PPG Beats ---
        if (!rIsNoise && !data[i].ecgRIndex.empty()) {
            try {
                data[i].pairs = pairRtoPPGBeat(seg.ecg, seg.po,
                    seg.ecgSampleRate, seg.ppgSampleRate,
                    data[i].ecgRIndex, data[i].ppgMinAmps);
            }
            catch (...) {
                // Fallback for failed pairing: match PPG to a dummy ECG index (-1)
                if (!data[i].bad_segment) {
                    data[i].ecgRIndex.clear();
                    data[i].pairs.assign(data[i].ppgMinAmps.size(), { 0.0, -1.0, 0.0, 0.0 });
                    for (std::size_t j = 0; j < data[i].ppgMinAmps.size(); ++j) {
                        data[i].pairs[j][0] = static_cast<double>(data[i].ppgMinAmps[j]);
                    }
                }
            }
        }
        else {
            // Case: Noisy ECG or Empty ECG, but valid PPG exists
            if (!data[i].bad_segment && !data[i].ppgMinAmps.empty()) {
                data[i].ecgRIndex.clear();
                data[i].pairs.assign(data[i].ppgMinAmps.size(), { 0.0, -1.0, 0.0, 0.0 });
                for (std::size_t j = 0; j < data[i].ppgMinAmps.size(); ++j) {
                    data[i].pairs[j][0] = static_cast<double>(data[i].ppgMinAmps[j]);
                }
            }
        }
    }

    return data;
}
