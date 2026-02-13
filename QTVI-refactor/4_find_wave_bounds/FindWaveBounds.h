#ifndef FINDWAVEBOUNDS_H
#define FINDWAVEBOUNDS_H

#include <vector>
#include <string>
#include <cstddef> // 1. Added for std::size_t

struct AnnealedSegment {
    std::vector<double> ecg;
    std::vector<double> po;
    double ecgSampleRate;
    double ppgSampleRate;
    std::vector<std::size_t> r_peaks;        // 2. Added std:: prefix
    std::vector<std::size_t> ppg_bin_indexs;
    std::vector<std::size_t> ecg_bin_indexs;
};

struct AnnealedData {
    std::vector<AnnealedSegment> bins;
};

struct WaveData {
    std::vector<std::vector<double>> pairs;
    std::vector<double> ecgPeaks;
    bool bad_segment;
    std::vector<double> ecgSeg;
    std::vector<double> ppgSeg;
    std::vector<std::size_t> ecgRIndex;
    std::vector<std::size_t> ppgMinAmps;
    std::vector<std::size_t> ppgMaxAmps;
    std::size_t index;
    double ecgSamplingRate;
    double ppgSamplingRate;
    std::vector<std::size_t> ppg_bin_indexs;
    std::vector<std::size_t> ecg_bin_indexs;
};

AnnealedData readCppBin(const std::string& path, double ecgFs, double ppgFs);
void saveWaveData(const std::string& path, const std::vector<WaveData>& results);

#endif
