#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <filesystem>
#include <algorithm>
#include <stdexcept>
#include <cstdint> // For uint64_t

// This include brings in the structures automatically
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

// --- Binary Loading Logic ---
AnnealedData readCppBin(const std::string& path, double ecgFs, double ppgFs) {
    AnnealedData data;
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("Could not open: " + path);

    uint64_t numBins = 0;
    if (!file.read(reinterpret_cast<char*>(&numBins), 8)) return data;
    data.bins.resize(numBins);

    for (uint64_t i = 0; i < numBins; ++i) {
        uint64_t pS, eS, sS;
        
        if (!file.read(reinterpret_cast<char*>(&pS), 8)) break;
        data.bins[i].po.resize(pS);
        file.read(reinterpret_cast<char*>(data.bins[i].po.data()), pS * 8);

        if (!file.read(reinterpret_cast<char*>(&eS), 8)) break;
        data.bins[i].ecg.resize(eS);
        file.read(reinterpret_cast<char*>(data.bins[i].ecg.data()), eS * 8);

        data.bins[i].ecgSampleRate = ecgFs;
        data.bins[i].ppgSampleRate = ppgFs;

        if (!file.read(reinterpret_cast<char*>(&sS), 8)) break;
        std::vector<double> sleep_data(sS);
        if (sS > 0) {
            file.read(reinterpret_cast<char*>(sleep_data.data()), sS * 8);
        }
    }
    return data;
}

// --- Binary Saving Logic ---
void saveWaveData(const std::string& path, const std::vector<WaveData>& results) {
    std::ofstream file(path, std::ios::binary);
    if (!file.is_open()) return;

    // 1. Write the number of segments (bins)
    uint64_t numBins = results.size();
    file.write(reinterpret_cast<const char*>(&numBins), 8);

    for (const auto& bin : results) {
        // 2. Save Sampling Rates (Matches Matlab metadata)
        file.write(reinterpret_cast<const char*>(&bin.ecgSamplingRate), 8);
        file.write(reinterpret_cast<const char*>(&bin.ppgSamplingRate), 8);

        // Helper: write indices with +1 offset for MATLAB compatibility
        auto writeIdx = [&](const std::vector<size_t>& v) {
            uint64_t sz = v.size();
            file.write(reinterpret_cast<const char*>(&sz), 8);
            for (size_t val : v) {
                uint64_t out = static_cast<uint64_t>(val) + 1; // Convert to 1-based
                file.write(reinterpret_cast<const char*>(&out), 8);
            }
            };

        // 3. Save Feature Indices
        writeIdx(bin.ecgRIndex);
        writeIdx(bin.ppgMaxAmps);
        writeIdx(bin.ppgMinAmps);

        // Helper: write raw signals
        auto writeSignal = [&](const std::vector<double>& sig) {
            uint64_t sz = sig.size();
            file.write(reinterpret_cast<const char*>(&sz), 8);
            if (sz > 0) {
                file.write(reinterpret_cast<const char*>(sig.data()), sz * 8);
            }
            };

        // 4. Save Raw Signals (Both PPG and ECG)
        writeSignal(bin.ppgSignal);
        writeSignal(bin.ecgSeg);    // ADDED: Saves the raw ECG signal waveform

        // 5. Save Pairs (Mapping)
        uint64_t numPairs = bin.pairs.size();
        file.write(reinterpret_cast<const char*>(&numPairs), 8);
        for (const auto& p : bin.pairs) {
            // Convert indices to 1-based, keep -1 for unmatched points
            int64_t ppg_out = (std::isnan(p[0]) || p[0] < -0.1) ? -1 : (int64_t)std::round(p[0]) + 1;
            int64_t ecg_out = (std::isnan(p[1]) || p[1] < -0.1) ? -1 : (int64_t)std::round(p[1]) + 1;
            file.write(reinterpret_cast<const char*>(&ppg_out), 8);
            file.write(reinterpret_cast<const char*>(&ecg_out), 8);
        }
    }
}


// --- Config Parser ---
ConfigSettings parseConfig(const std::string& configPath, const int type_row) {
    std::ifstream file(configPath);
    std::string line; std::getline(file, line);
    int currentIdx = 0; // Track the current data row index
    while (std::getline(file, line)) {
        // If the current index matches the requested type_row (0-indexed)
        if (currentIdx++ == type_row - 1) {
            std::stringstream ss(line);
            std::vector<std::string> row;
            std::string val;
            while (std::getline(ss, val, ','))
                row.push_back(val);

            return { row[0], row[6], row[7], std::stod(row[12]),
                   (row.size() > 13 && !row[13].empty()) ? std::stod(row[13]) : std::stod(row[12]) };
        }
    }

    throw std::runtime_error("Type not found");
}

int main() {
    try {
        int choice;
        std::cout << "Enter Data Type:\n1) MESA\n2) Bittium\n3) CHAOS" << std::endl;
        std::cin >> choice;
        ConfigSettings cfg = parseConfig("config.csv", choice);

        if (!fs::exists(cfg.wavePath)) fs::create_directories(cfg.wavePath);

        for (const auto& entry : fs::directory_iterator(cfg.annealedPath)) {
            if (entry.path().extension() == ".bin") {
                AnnealedData annealedData = readCppBin(entry.path().string(), cfg.ecgFs, cfg.ppgFs);

                // Process logic
                auto results = FindWaveBounds_EKGandPPG(annealedData.bins, 0, true, entry.path().stem().string());

                // Save the results, now including ppgMinAmps and pairs, to the binary file
                saveWaveData(cfg.wavePath + "/" + entry.path().stem().string() + "_wave_data.bin", results);
                std::cout << "Processed and saved binary for: " << entry.path().filename() << std::endl;
            }
        }
    }
    catch (const std::exception& e) { std::cerr << e.what() << std::endl; }
    return 0;
}
