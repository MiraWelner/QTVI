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
        if (sS > 0) file.seekg(sS * 8, std::ios_base::cur); // Skip broken sleep data
    }
    return data;
}

// --- Binary Saving Logic ---
void saveWaveData(const std::string& path, const std::vector<WaveData>& results) {
    std::ofstream file(path, std::ios::binary);
    if (!file.is_open()) return;

    uint64_t numBins = results.size();
    file.write(reinterpret_cast<char*>(&numBins), 8);

    for (const auto& bin : results) {
        // --- SECTION 1: ECG R-PEAKS ---
        uint64_t numRPeaks = bin.ecgRIndex.size();
        file.write(reinterpret_cast<char*>(&numRPeaks), 8);
        for (size_t rIdx : bin.ecgRIndex) {
            uint64_t val = static_cast<uint64_t>(rIdx);
            file.write(reinterpret_cast<char*>(&val), 8);
        }

        // --- SECTION 2: PPG SYSTOLIC PEAKS (ppgMaxAmps) ---
        uint64_t numMaxAmps = bin.ppgMaxAmps.size();
        file.write(reinterpret_cast<char*>(&numMaxAmps), 8);
        for (size_t pIdx : bin.ppgMaxAmps) {
            uint64_t val = static_cast<uint64_t>(pIdx);
            file.write(reinterpret_cast<char*>(&val), 8);
        }

        // --- SECTION 3: PPG MIN AMPS (ppgMinAmps) ---
        uint64_t numMinAmps = bin.ppgMinAmps.size();
        file.write(reinterpret_cast<char*>(&numMinAmps), 8);
        for (size_t mIdx : bin.ppgMinAmps) {
            uint64_t val = static_cast<uint64_t>(mIdx);
            file.write(reinterpret_cast<char*>(&val), 8);
        }

        // --- SECTION 4: PAIRS ---
        uint64_t numPairs = bin.pairs.size();
        file.write(reinterpret_cast<char*>(&numPairs), 8);
        for (const auto& pair : bin.pairs) {
            // Assuming each pair always has 2 elements
            if (pair.size() == 2) {
                uint64_t ppg_idx = static_cast<uint64_t>(pair[0]);
                uint64_t ecg_idx = static_cast<uint64_t>(pair[1]);
                file.write(reinterpret_cast<char*>(&ppg_idx), 8);
                file.write(reinterpret_cast<char*>(&ecg_idx), 8);
            }
            else {
                // If a pair is malformed, write two zeros as placeholders
                uint64_t zero = 0;
                file.write(reinterpret_cast<char*>(&zero), 8);
                file.write(reinterpret_cast<char*>(&zero), 8);
            }
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
                auto results = FindWaveBounds_EKGandPPG(annealedData.bins, cfg.ecgFs, true);

                // Save the results, now including ppgMinAmps and pairs, to the binary file
                saveWaveData(cfg.wavePath + "/" + entry.path().stem().string() + "_wave_data.bin", results);
                std::cout << "Processed and saved binary for: " << entry.path().filename() << std::endl;
            }
        }
    }
    catch (const std::exception& e) { std::cerr << e.what() << std::endl; }
    return 0;
}
