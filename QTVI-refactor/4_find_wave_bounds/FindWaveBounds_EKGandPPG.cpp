// ============================================================================
// File: FindWaveBounds_EKGandPPG.cpp
// ============================================================================
#include "FindWaveBounds_EKGandPPG.h"
#include "SegmentPPG.h"
#include "JoinedRR.h"
#include "pairRtoPPGBeat.h"
#include "StatsUtils.h"
#include <iostream>

vector<WaveData> FindWaveBounds_EKGandPPG(const vector<AnnealedSegment>& annealedSegments,
    int dbg_plot,
    bool use_R_algorithms) {
    vector<WaveData> data(annealedSegments.size());

    for (size_t i = 0; i < annealedSegments.size(); ++i) {
        double ppgSamplingRate = annealedSegments[i].ppgSampleRate;
        double ecgSamplingRate = annealedSegments[i].ecgSampleRate;
        vector<double> ecgSeg = annealedSegments[i].ecg;
        bool rIsNoise = false;

        vector<size_t> ecgRIndex;
        /*
        try {
            ecgRIndex = annealedSegments[i].r_peaks;
        }
        catch (...) {
            ecgRIndex.clear();
        }
        */

        vector<double> ppgSeg = annealedSegments[i].po;

        vector<size_t> ppgMinAmps, ppgMaxAmps;
        try {
            auto ppgResult = SegmentPPG(ppgSeg, ppgSamplingRate);
            ppgMinAmps = ppgResult.ppgMinAmps;
            ppgMaxAmps = ppgResult.maxAmps;
            data[i].bad_segment = false;
        }
        catch (...) {
            ppgMinAmps.clear();
            ppgMaxAmps.clear();
            data[i].bad_segment = true;
        }
        if (ecgRIndex.empty() && use_R_algorithms) {
            if (std_dev(ecgSeg) == 0) {
                rIsNoise = true;
            }
        }

        if (!rIsNoise && ecgRIndex.empty() && use_R_algorithms) {
            try {
                ecgRIndex = JoinedRR(ecgSeg, ecgSamplingRate, 2.0);
                if (ecgRIndex.size() < ppgMinAmps.size() / 2 ||
                    ppgMinAmps.size() * 1.5 < ecgRIndex.size()) {
                    rIsNoise = true;
                }
            }
            catch (...) {
                rIsNoise = true;
            }
        }

        vector<vector<double>> pairs;
        if (!rIsNoise && !ecgRIndex.empty()) {
            try {
                pairs = pairRtoPPGBeat(ecgSeg, ppgSeg, ecgSamplingRate, ppgSamplingRate,
                    ecgRIndex, ppgMinAmps);
            }
            catch (...) {
                if (!data[i].bad_segment) {
                    //ecgRIndex.clear();
                    pairs.resize(ppgMinAmps.size(), vector<double>(2));
                    for (size_t j = 0; j < ppgMinAmps.size(); ++j) {
                        pairs[j][0] = ppgMinAmps[j];
                        pairs[j][1] = -1;
                    }
                }
                else {
                    pairs.clear();
                }
            }
        }
        else {
            if (!data[i].bad_segment) {
               // ecgRIndex.clear();
                pairs.resize(ppgMinAmps.size(), vector<double>(2));
                for (size_t j = 0; j < ppgMinAmps.size(); ++j) {
                    pairs[j][0] = ppgMinAmps[j];
                    pairs[j][1] = -1;
                }
            }
            else {
                pairs.clear();
            }
        }
        // Populate data structure
        data[i].ecgSeg = ecgSeg;
        data[i].ppgSeg = ppgSeg;
        data[i].ecgRIndex = ecgRIndex;
        data[i].ppgMinAmps = ppgMinAmps;
        data[i].ppgMaxAmps = ppgMaxAmps;
        data[i].pairs = pairs;
        data[i].index = i;
        data[i].ecgSamplingRate = ecgSamplingRate;
        data[i].ppgSamplingRate = ppgSamplingRate;
        data[i].ppg_bin_indexs = annealedSegments[i].ppg_bin_indexs;
        data[i].ecg_bin_indexs = annealedSegments[i].ecg_bin_indexs;

    }

    return data;
}
