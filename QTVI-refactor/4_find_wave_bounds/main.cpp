// ============================================================================
// File: main.cpp
// Entry point: reads annealed binary segments, detects R-peaks and PPG
// pulses, pairs them, and writes results to binary output.
// ============================================================================
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <filesystem>
#include <algorithm>
#include <stdexcept>
#include <cstdint>
#include <cmath>

#include "FindWaveBounds.h"
#include "FindWaveBounds_EKGandPPG.h"

namespace fs = std::filesystem;

struct ConfigSettings {
    std::string dataType;
    std::string annealedPath;
    std::string wavePath;
    double ecgFs;
    double ppgFs;
};

// ============================================================================
// Binary I/O
// ============================================================================

AnnealedData readCppBin(const std::string& path) {
    AnnealedData data;
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("Could not open: " + path);

    uint64_t numBins = 0;
    file.read(reinterpret_cast<char*>(&numBins), 8);

    double filePpgSR, fileEcgSR, scoringEpoch;
    file.read(reinterpret_cast<char*>(&filePpgSR), 8);
    file.read(reinterpret_cast<char*>(&fileEcgSR), 8);
    file.read(reinterpret_cast<char*>(&scoringEpoch), 8);

    data.bins.resize(numBins);

    for (uint64_t i = 0; i < numBins; ++i) {
        auto& bin = data.bins[i];

        // PPG index pairs
        uint64_t nPpgPairs;
        if (!file.read(reinterpret_cast<char*>(&nPpgPairs), 8)) break;
        bin.ppg_bin_indexs.resize(nPpgPairs);
        for (uint64_t j = 0; j < nPpgPairs; ++j) {
            file.read(reinterpret_cast<char*>(&bin.ppg_bin_indexs[j].first), 8);
            file.read(reinterpret_cast<char*>(&bin.ppg_bin_indexs[j].second), 8);
        }

        // ECG index pairs
        uint64_t nEcgPairs;
        if (!file.read(reinterpret_cast<char*>(&nEcgPairs), 8)) break;
        bin.ecg_bin_indexs.resize(nEcgPairs);
        for (uint64_t j = 0; j < nEcgPairs; ++j) {
            file.read(reinterpret_cast<char*>(&bin.ecg_bin_indexs[j].first), 8);
            file.read(reinterpret_cast<char*>(&bin.ecg_bin_indexs[j].second), 8);
        }

        // Helper: read a double array preceded by its uint64 count
        auto readDoubleArray = [&](std::vector<double>& vec) -> bool {
            uint64_t sz;
            if (!file.read(reinterpret_cast<char*>(&sz), 8)) return false;
            vec.resize(sz);
            if (sz > 0) file.read(reinterpret_cast<char*>(vec.data()), sz * 8);
            return true;
            };

        if (!readDoubleArray(bin.po)) break;
        if (!readDoubleArray(bin.ecg)) break;
        if (!readDoubleArray(bin.ecg2)) break;
        if (!readDoubleArray(bin.ecg3)) break;
        if (!readDoubleArray(bin.sleepStages)) break;

        bin.ecgSampleRate = fileEcgSR;
        bin.ppgSampleRate = filePpgSR;
    }
    return data;
}

void saveWaveData(const std::string& path, const std::vector<WaveData>& results) {
    std::ofstream file(path, std::ios::binary);
    if (!file.is_open()) return;

    uint64_t numBins = results.size();
    file.write(reinterpret_cast<const char*>(&numBins), 8);

    for (const auto& bin : results) {
        file.write(reinterpret_cast<const char*>(&bin.ecgSamplingRate), 8);
        file.write(reinterpret_cast<const char*>(&bin.ppgSamplingRate), 8);

        // Write a size_t index array, converting to 1-based for MATLAB compatibility
        auto writeIdx = [&](const std::vector<std::size_t>& v) {
            uint64_t sz = v.size();
            file.write(reinterpret_cast<const char*>(&sz), 8);
            for (std::size_t val : v) {
                uint64_t out = static_cast<uint64_t>(val) + 1;
                file.write(reinterpret_cast<const char*>(&out), 8);
            }
            };

        auto writeSignal = [&](const std::vector<double>& sig) {
            uint64_t sz = sig.size();
            file.write(reinterpret_cast<const char*>(&sz), 8);
            if (sz > 0) file.write(reinterpret_cast<const char*>(sig.data()), sz * 8);
            };

        writeIdx(bin.ecgRIndex);
        writeIdx(bin.ecgRIndex2);
        writeIdx(bin.ecgRIndex3);
        writeIdx(bin.ppgMaxAmps);
        writeIdx(bin.ppgMinAmps);

        writeSignal(bin.ppgSignal);
        writeSignal(bin.ecgSignal);
        writeSignal(bin.ecgSignal2);
        writeSignal(bin.ecgSignal3);

        // Pairs: [ppg_idx, ecg_idx], converted to 1-based (-1 for invalid)
        uint64_t numPairs = bin.pairs.size();
        file.write(reinterpret_cast<const char*>(&numPairs), 8);
        for (const auto& p : bin.pairs) {
            int64_t ppg_out = (std::isnan(p[0]) || p[0] < -0.1) ? -1 : static_cast<int64_t>(std::round(p[0])) + 1;
            int64_t ecg_out = (std::isnan(p[1]) || p[1] < -0.1) ? -1 : static_cast<int64_t>(std::round(p[1])) + 1;
            file.write(reinterpret_cast<char*>(&ppg_out), 8);
            file.write(reinterpret_cast<char*>(&ecg_out), 8);
        }
    }
}

// ============================================================================
// Config
// ============================================================================

ConfigSettings parseConfig(const std::string& configPath, int type_row) {
    std::ifstream file(configPath);
    if (!file.is_open()) throw std::runtime_error("Could not open config.csv");

    std::string line;
    std::getline(file, line); // Skip header

    int currentIdx = 0;
    while (std::getline(file, line)) {
        if (currentIdx++ == type_row - 1) {
            std::stringstream ss(line);
            std::vector<std::string> row;
            std::string val;
            while (std::getline(ss, val, ',')) row.push_back(val);

            return { row[0], row[6], row[7], std::stod(row[12]),
                   (row.size() > 13 && !row[13].empty()) ? std::stod(row[13]) : std::stod(row[12]) };
        }
    }
    throw std::runtime_error("Config type row not found");
}

// ============================================================================
// Main
// ============================================================================

int main() {
    try {
        int choice;
        std::cout << "Enter Data Type:\n1) MESA\n2) Bittium\n3) CHAOS" << std::endl;
        if (!(std::cin >> choice)) return 1;

        ConfigSettings cfg = parseConfig("config.csv", choice);

        if (!fs::exists(cfg.wavePath)) {
            fs::create_directories(cfg.wavePath);
        }

        for (const auto& entry : fs::directory_iterator(cfg.annealedPath)) {
            if (entry.path().extension() == ".bin") {
                std::string currentID = entry.path().stem().string();
                std::cout << "Processing: " << currentID << "..." << std::endl;

                AnnealedData annealedData = readCppBin(entry.path().string());

                std::vector<WaveData> results = FindWaveBounds_EKGandPPG(
                    annealedData.bins, 0, true, currentID);

                std::string outputPath = cfg.wavePath + "/" + currentID + "_wave_data.bin";
                saveWaveData(outputPath, results);
            }
        }
        std::cout << "Processing Complete." << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "CRITICAL ERROR: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}