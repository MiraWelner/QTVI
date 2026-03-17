#ifndef FINDWAVEBOUNDS_H
#define FINDWAVEBOUNDS_H

#include <vector>
#include <string>
#include <cstddef>

struct AnnealedSegment {
    std::vector<double> po;
    std::vector<double> ecg;
    std::vector<double> ecg2;
    std::vector<double> ecg3;
    std::vector<double> sleepStages;
    std::vector<std::pair<uint64_t, uint64_t>> ppg_bin_indexs;
    std::vector<std::pair<uint64_t, uint64_t>> ecg_bin_indexs;
    double ecgSampleRate;
    double ppgSampleRate;
};

struct AnnealedData {
    std::vector<AnnealedSegment> bins;
};

struct WaveData {
    std::vector<std::vector<double>> pairs;
    std::vector<double> ecgPeaks;
    bool bad_segment;

    // Channel 1 (always used)
    std::vector<double> ecgSignal;
    std::vector<std::size_t> ecgRIndex;

    // Channel 2 (populated only for multi-channel inputs like Bittium/CHAOS)
    std::vector<double> ecgSignal2;
    std::vector<std::size_t> ecgRIndex2;

    // Channel 3 (populated only for multi-channel inputs like Bittium/CHAOS)
    std::vector<double> ecgSignal3;
    std::vector<std::size_t> ecgRIndex3;

    std::vector<std::size_t> ppgMinAmps;
    std::vector<std::size_t> ppgMaxAmps;
    std::size_t index;
    double ecgSamplingRate;
    double ppgSamplingRate;
    std::vector<std::pair<uint64_t, uint64_t>> ppg_bin_indexs;
    std::vector<std::pair<uint64_t, uint64_t>> ecg_bin_indexs;
    std::vector<double> ppgSignal;
};

AnnealedData readCppBin(const std::string& path, double ecgFs, double ppgFs);
void saveWaveData(const std::string& path, const std::vector<WaveData>& results);

#endif